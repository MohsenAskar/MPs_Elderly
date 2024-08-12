**********************************************************************************
* 1st exclusion 
**********************************************************************************
// Drop non elderly 
use "\mod_file.dta", replace
drop kodenavn
drop if alder < 65
// Removing entries with miscoded ICD codes
drop if kodeverdi==" "
// to remove irrelevant ICD COdes
drop if strpos( kodeverdi,"X9")
drop if strpos( kodeverdi,"Y4")
drop if strpos( kodeverdi,"X8")
drop if strpos( kodeverdi,"V4")
drop if strpos( kodeverdi,"Z3")
drop if strpos( kodeverdi,"X0n")
drop if strpos( kodeverdi,"Z0")
drop if strpos( kodeverdi,"T4n")
drop if strpos( kodeverdi,"Z1")

// removing patients that had only 1 diagnosis multiple times
sort lopenr
by lopenr : gen n_diseases = _N
sort lopenr kodeverdi
by lopenr kodeverdi : gen n_diseases2 = _n
by lopenr: gen okay =1 if  n_diseases2 == n_diseases
tab okay
xfill okay , i(lopenr)
drop if okay == 1
drop n_diseases2 n_diseases okay
// Removing patients if they have only 1 entry
bysort lopenr:drop if _N==1
// removing duplicates in terms of all variables
duplicates drop
save "D:\data_after_all.dta" 
save "\data_after_all.dta", replace

**********************************************
*  Determine chronic ICD-10 code using CCIR
**********************************************
// Fixing CCIR list 
import delimited "\Chronic_Condition_Indicatior_Refined_CCIR\CCIR_v2023-1.csv"
drop in 1/2
rename v1 ICD10CM_CODE
rename v2 ICD10CM_DESCRIPTION
rename v3 CHRONIC_INDICATOR
drop in 1
gen test = ""
replace test = "Not chronic" if CHRONIC_INDICATOR == "0"
replace test = "Chronic" if CHRONIC_INDICATOR == "1"
replace test = "No determination" if CHRONIC_INDICATOR == "9"
replace ICD10CM_CODE = substr( ICD10CM_CODE , 1, length( ICD10CM_CODE ) - 1)
preserve
replace ICD10CM_CODE = substr( ICD10CM_CODE , 2, .)
rename test CHRONIC_INDICATOR2
drop CHRONIC_INDICATOR
rename CHRONIC_INDICATOR2 CHRONIC_INDICATOR
gen ICD_3=substr(ICD10CM_CODE ,1,3)
save "\Chronic_Condition_Indicatior_Refined_CCIR\CCIR_List_Stata.dta", replace
keep if CHRONIC_INDICATOR == "Chronic"
codebook ICD_3
save "\Chronic_Condition_Indicatior_Refined_CCIR\CCIR_List_3_Character_ICD.dta", replace


// merge with chronic conditions from CCIR list 
use "\Chronic_Condition_Indicatior_Refined_CCIR\CCIR_List_3_Character_ICD.dta", replace
keep ICD_3 CHRONIC_INDICATOR
duplicates drop 
save "\Chronic_Condition_Indicatior_Refined_CCIR\CCIR_List_3_Character_ICD_To_Merge.dta"

use "\Analysis\data_after_all.dta", replace
rename kodeverdi ICD_3
merge m:1 ICD_3 using "\Chronic_Condition_Indicatior_Refined_CCIR\CCIR_List_3_Character_ICD_To_Merge.dta"
keep if _merge ==3
drop _merge
rename ICD_3 hovedtilstand
rename lopenr pasientlopenr
// save 
export delimited using "\Analysis\Only_Chronic_after_CCIR.csv", replace

// histogram number of morbidties7 percent in the population 
import delimited using "\Analysis\Only_Chronic_after_CCIR.csv", clear
sort pasientlopenr hovedtilstand
duplicates drop pasientlopenr hovedtilstand, force
by pasientlopenr : gen n_diseases = _N
preserve
keep pasientlopenr n_diseases
duplicates drop
drop if n_diseases == 1
tabstat n_diseases , stats(n mean median min max)
hist n_diseases, percent

// Move to python to calculate RR and Phi
// CODE IN PYTHON 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Load data from created in Stata to calculate RR and Phi
############################################################
# filter out disease pairs that occurs in only one patient and then calulate RR and Phi
#######################################################################################
import pandas as pd
import numpy as np
# Loading data created in Stata
data= pd.read_csv(r"\Analysis\Only_Chronic_after_CCIR.csv")
data.shape
data = data.rename(columns={'pasientlopenr': 'patientid', 'hovedtilstand': 'code_system'})
# Create a count column
data["Count"] = 1

data['code_system'].nunique()
# Create a pivot table
data_pivot = data.pivot_table(index="patientid", columns="code_system", values="Count", aggfunc="max").fillna(0).astype(int)
# Calculate observed co-occurrences
observed_cooccurrences = data_pivot.T.dot(data_pivot).values  # Actual count of co-occurrences

# Convert observed_cooccurrences to a DataFrame 
disease_names = data_pivot.columns
observed_df = pd.DataFrame(observed_cooccurrences, index=disease_names, columns=disease_names)

# Filter out disease pairs with co-occurrences in only one patient
filtered_observed_df = observed_df.where(observed_df > 1)
filtered_observed_df.columns
# count unique patients and ICD codes after filtering the disease pairs that only occured once
unique_icd_codes_after_filter = filtered_observed_df.index.unique().size
filtered_diseases = filtered_observed_df.index[filtered_observed_df.any(axis=1)]
unique_patients_after_filter = (data_pivot[filtered_diseases] > 0).any(axis=1).sum()


# Calculate incidence rates for filtered disease pairs
incidence_rates = data_pivot.sum() / data_pivot.shape[0]

# Expected and observed rates for filtered pairs
expected_rates = np.outer(incidence_rates, incidence_rates)
observed_rates = (filtered_observed_df / data_pivot.shape[0]).values

# RR, Phi calculations only for filtered pairs
rr_values = np.divide(observed_rates, expected_rates, out=np.zeros_like(observed_rates), where=expected_rates!=0)
chi2_values = np.divide(np.square(observed_rates * data_pivot.shape[0] - expected_rates * data_pivot.shape[0]), expected_rates * data_pivot.shape[0], out=np.zeros_like(observed_rates), where=expected_rates!=0)
phi_values = np.sqrt(chi2_values / data_pivot.shape[0])

# Apply filters to RR and Phi to remove non-significant values
rr_filtered = np.where(rr_values > 1, rr_values, np.nan)
phi_filtered = np.where(phi_values > 0, phi_values, np.nan)

# Compile results, including patient count for each filtered pair
result = []
for i in range(len(disease_names)):
    for j in range(i+1, len(disease_names)):
        if not np.isnan(rr_filtered[i,j]) and not np.isnan(phi_filtered[i,j]):
            count = filtered_observed_df.iloc[i, j]  # The count of patients with both diseases
            if not np.isnan(count):  # Ensure we're only including pairs that weren't filtered out
                result.append((disease_names[i], disease_names[j], rr_filtered[i,j], phi_filtered[i,j], count))

result_df = pd.DataFrame(result, columns=["Disease1", "Disease2", "RR", "Phi", "PatientCount"])

# Display the resulting DataFrame with RR and Phi values for disease pairs
print(result_df)
result_df['PatientCount'].describe()
result_df['Disease1'].nunique() #641
result_df['Disease2'].nunique()#682
result_df.to_csv(r"\Analysis\reults_df_RR_Phi.csv", index= False)

# counting uniques patiets and ICD codes
unique_diseases_result_df = pd.unique(result_df[['Disease1', 'Disease2']].values.ravel()).size
# Step 1: Extract disease pairs from 'result_df'
disease_pairs = [(row['Disease1'], row['Disease2']) for index, row in result_df.iterrows()]

# Step 2: Map disease pairs to patients
# This step assumes 'data_pivot' is structured with patients as rows and diseases as columns
unique_patients = set()  # To store unique patient identifiers

for disease1, disease2 in disease_pairs:
    # Find patients with both diseases
    patients_with_pair = data_pivot[(data_pivot[disease1] > 0) & (data_pivot[disease2] > 0)].index
    # Add these patients to the set of unique patients
    unique_patients.update(patients_with_pair)

# Step 3: Count unique patients
unique_patient_count = len(unique_patients)

# to calculate the uniqiue patients and ICD codes for the exclusion flow diagram
###################################################################################
# Now, apply the additional exclusion criteria
# Exclude ICD-10 codes that occur in fewer than 50 patients
icd_counts = data_pivot.sum(axis=0)
filtered_icd_codes = icd_counts[icd_counts >= 50].index
data_pivot_filtered = data_pivot[filtered_icd_codes]

# Recalculate observed co-occurrences with the filtered data
observed_cooccurrences_filtered = data_pivot_filtered.T.dot(data_pivot_filtered).values
disease_names_filtered = data_pivot_filtered.columns
observed_df_filtered = pd.DataFrame(observed_cooccurrences_filtered, index=disease_names_filtered, columns=disease_names_filtered)

# Filter pairs of ICD-10 codes with co-occurrences in fewer than 5 patients
final_filtered_observed_df = observed_df_filtered.where(observed_df_filtered >= 5)

# Calculate incidence rates for the final filtered disease pairs
incidence_rates_filtered = data_pivot_filtered.sum() / data_pivot_filtered.shape[0]

# Expected and observed rates for final filtered pairs
expected_rates_filtered = np.outer(incidence_rates_filtered, incidence_rates_filtered)
observed_rates_filtered = (final_filtered_observed_df / data_pivot_filtered.shape[0]).values

# RR, Phi calculations for final filtered pairs
rr_values_filtered = np.divide(observed_rates_filtered, expected_rates_filtered, out=np.zeros_like(observed_rates_filtered), where=expected_rates_filtered!=0)
chi2_values_filtered = np.divide(np.square(observed_rates_filtered * data_pivot_filtered.shape[0] - expected_rates_filtered * data_pivot_filtered.shape[0]), expected_rates_filtered * data_pivot_filtered.shape[0], out=np.zeros_like(observed_rates_filtered), where=expected_rates_filtered!=0)
phi_values_filtered = np.sqrt(chi2_values_filtered / data_pivot_filtered.shape[0])

# Apply filters to RR and Phi to remove non-significant values
rr_final_filtered = np.where(rr_values_filtered > 1, rr_values_filtered, np.nan)
phi_final_filtered = np.where(phi_values_filtered > 0, phi_values_filtered, np.nan)

# Compile final results, including patient count for each filtered pair
final_result = []
for i in range(len(disease_names_filtered)):
    for j in range(i+1, len(disease_names_filtered)):
        if not np.isnan(rr_final_filtered[i,j]) and not np.isnan(phi_final_filtered[i,j]):
            count = final_filtered_observed_df.iloc[i, j]
            if not np.isnan(count):
                final_result.append((disease_names_filtered[i], disease_names_filtered[j], rr_final_filtered[i,j], phi_final_filtered[i,j], count))

final_result_df = pd.DataFrame(final_result, columns=["Disease1", "Disease2", "RR", "Phi", "PatientCount"])

# Display the final resulting DataFrame with RR and Phi values for disease pairs
print(final_result_df)
final_result_df['PatientCount'].describe()
final_result_df['Disease1'].nunique()
final_result_df['Disease2'].nunique()
final_result_df.to_csv(r"\Analysis\final_results_df_RR_Phi.csv", index=False)

# Count unique patients and ICD codes after final filtering
unique_diseases_final_result_df = pd.unique(final_result_df[['Disease1', 'Disease2']].values.ravel()).size

# Extract final disease pairs from 'final_result_df'
final_disease_pairs = [(row['Disease1'], row['Disease2']) for index, row in final_result_df.iterrows()]

# Map final disease pairs to patients
final_unique_patients = set()

for disease1, disease2 in final_disease_pairs:
    patients_with_pair = data_pivot_filtered[(data_pivot_filtered[disease1] > 0) & (data_pivot_filtered[disease2] > 0)].index
    final_unique_patients.update(patients_with_pair)

# Count unique patients after final filtering
final_unique_patient_count = len(final_unique_patients)

# Display final unique patient and disease counts
print(f"Unique ICD codes: {unique_diseases_final_result_df}")
print(f"Unique patients: {final_unique_patient_count}")

// separate files for men, women were created in the same folder

// END OF CODE IN PYTHON, BACK TO STATA
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// after keeping RR > 1 and phi >0  and saving the results in Python
import delimited "\Analysis\reults_df_RR_Phi_Filtered_50_Patients_5_Edgeweight.csv", clear 

// generate the networks 
preserve
drop if disease1 == disease2
duplicates drop
nwfromedge disease1 disease2 patientcount, undirected keeporiginal name(network_all)
nwexport, type(pajek) replace


// Group codes after ICD hirarchy to color them afterwards
import delimited using "\Analysis\Only_Chronic_after_CCIR.csv", clear
gen ICD1=substr( hovedtilstand ,1,1)
gen ICD2=substr( hovedtilstand ,1,2)
gen icdgrp=0
replace icdgrp=1 if ICD1=="A" | ICD1=="B"
replace icdgrp=2 if ICD2=="C0" | ICD2=="C1"|ICD2=="C2"| ICD2=="C3"| ICD2=="C4"| ICD2=="C5"| ICD2=="C6"| ICD2=="C7"| ICD2=="C8"| ICD2=="C9"| ICD2=="D0"|ICD2=="D1"| ICD2=="D2"| ICD2=="D3"| ICD2=="D4"
replace icdgrp=3 if ICD2=="D5"| ICD2=="D6"| ICD2=="D7"| ICD2=="D8"
replace icdgrp=4 if ICD2=="E0"| ICD2=="E1"| ICD2=="E2"| ICD2=="E3"| ICD2=="E4"| ICD2=="E5"| ICD2=="E6"| ICD2=="E7"| ICD2=="E8"| ICD2=="E9"
replace icdgrp=5 if ICD2=="F0"| ICD2=="F1"| ICD2=="F2"| ICD2=="F3"| ICD2=="F4"| ICD2=="F5"| ICD2=="F6"| ICD2=="F7"| ICD2=="F8"| ICD2=="F9"
replace icdgrp=6 if ICD2=="G0"| ICD2=="G1"| ICD2=="G2"| ICD2=="G3"| ICD2=="G4"| ICD2=="G5"| ICD2=="G6"| ICD2=="G7"| ICD2=="G8"| ICD2=="G9"
replace icdgrp=7 if ICD2=="H0"| ICD2=="H1"| ICD2=="H2"| ICD2=="H3"| ICD2=="H4"| ICD2=="H5"
replace icdgrp=8 if ICD2=="H6"| ICD2=="H7"| ICD2=="H8"| ICD2=="H9"
replace icdgrp=9 if ICD2=="I0"| ICD2=="I1"| ICD2=="I2"| ICD2=="I3"| ICD2=="I4"| ICD2=="I5"| ICD2=="I6"| ICD2=="I7"| ICD2=="I8"| ICD2=="I9"
replace icdgrp=10 if ICD2=="J0"| ICD2=="J1"| ICD2=="J2"| ICD2=="J3"| ICD2=="J4"| ICD2=="J5"| ICD2=="J6"| ICD2=="J7"| ICD2=="J8"| ICD2=="J9"
replace icdgrp=11 if ICD2=="K0"| ICD2=="K1"| ICD2=="K2"| ICD2=="K3"| ICD2=="K4"| ICD2=="K5"| ICD2=="K6"| ICD2=="K7"| ICD2=="K8"| ICD2=="K9"
replace icdgrp=12 if ICD2=="L0"| ICD2=="L1"| ICD2=="L2"| ICD2=="L3"| ICD2=="L4"| ICD2=="L5"| ICD2=="L6"| ICD2=="L7"| ICD2=="L8"| ICD2=="L9"
replace icdgrp=13 if ICD2=="M0"| ICD2=="M1"| ICD2=="M2"| ICD2=="M3"| ICD2=="M4"| ICD2=="M5"| ICD2=="M6"| ICD2=="M7"| ICD2=="M8"| ICD2=="M9"
replace icdgrp=14 if ICD2=="N0"| ICD2=="N1"| ICD2=="N2"| ICD2=="N3"| ICD2=="N4"| ICD2=="N5"| ICD2=="N6"| ICD2=="N7"| ICD2=="N8"| ICD2=="N9"
replace icdgrp=15 if ICD2=="O0"| ICD2=="O1"| ICD2=="O2"| ICD2=="O3"| ICD2=="O4"| ICD2=="O5"| ICD2=="O6"| ICD2=="O7"| ICD2=="O8"| ICD2=="O9"
replace icdgrp=16 if ICD2=="P0"| ICD2=="P1"| ICD2=="P2"| ICD2=="P3"| ICD2=="P4"| ICD2=="P5"| ICD2=="P6"| ICD2=="P7"| ICD2=="P8"| ICD2=="P9"
replace icdgrp=17 if ICD2=="Q0"| ICD2=="Q1"| ICD2=="Q2"| ICD2=="Q3"| ICD2=="Q4"| ICD2=="Q5"| ICD2=="Q6"| ICD2=="Q7"| ICD2=="Q8"| ICD2=="Q9"
replace icdgrp=18 if ICD2=="R"| ICD2=="R1"| ICD2=="R2"| ICD2=="R3"| ICD2=="R4"| ICD2=="R5"| ICD2=="R6"| ICD2=="R7"| ICD2=="R8"| ICD2=="R9"
replace icdgrp=19 if ICD2=="S0"| ICD2=="S1"| ICD2=="S2"| ICD2=="S3"| ICD2=="S4"| ICD2=="S5"| ICD2=="S6"| ICD2=="S7"| ICD2=="S8"| ICD2=="S9"|ICD2=="T0"| ICD2=="T1"| ICD2=="T2"| ICD2=="T3"| ICD2=="T4"| ICD2=="T5"| ICD2=="T6"| ICD2=="T7"| ICD2=="T8"| ICD2=="T9"

// colour the ICD groups 
gen color="."
replace color ="#4cc2c9" if icdgrp==1 // I
replace color ="#9d4f46" if icdgrp==2 // II
replace color ="#e59fc6" if icdgrp==3 // III
replace color ="#df4497" if icdgrp==4 // IV
replace color ="#4bb649" if icdgrp==5 // V
replace color ="#3a5ca9" if icdgrp==6 // VI
replace color ="#fbc084" if icdgrp==7 // VII
replace color ="#a8a5d1" if icdgrp==8 // VIII
replace color ="#8c4198" if icdgrp==9 // IX
replace color ="#625e5e" if icdgrp==10 // X
replace color ="#0b8494" if icdgrp==11 // XI
replace color ="#729556" if icdgrp==12 // XII
replace color ="#942168" if icdgrp==13 // XIII
replace color ="#e9d215" if icdgrp==14 // XIV
replace color ="#196131" if icdgrp==15 // XV
replace color ="#bf9000" if icdgrp==16 // XVI
replace color ="#f4831f" if icdgrp==17 // XVII
replace color ="#27286f" if icdgrp==18 // XVIII
replace color ="#ec2224" if icdgrp==19 // XIX

export delimited using "\Analysis\Only_Chronic_after_CCIR.csv", replace

// attributes for import to Gephi
// add patients number of each ICD code
import delimited "\Analysis\Only_Chronic_after_CCIR.csv", clear 
preserve
duplicates drop hovedtilstand pasientlopenr, force
bysort hovedtilstand :egen no_patients_each_code =count(pasientlopenr)
sort no_patients_each_code
keep hovedtilstand icdgrp color no_patients_each_code
duplicates drop
export delimited using "\Analysis\ICD_Code_Chapter_Color_No_Users.csv", replace
save "\Analysis\ICD_Code_Chapter_Color_No_Users.dta", replace
// same was done to men, women files
// move to Gephi for visulaizing  
// File reults_df_RR_Phi.csv, can we make the network by patient counts, RR or Phi 
// File Chronic_ICD_Codes_Description.csv, can we get the ICD description
// File ICD_Code_Chapter_Color_No_Users.csv, include ICD_chapter, color code, no. of patients for each ICD

// add attributes to gephi
// merge color, no_users and ICD chapter
import delimited "\Analysis\From_Gephi\Network_modules_from_Gephi.csv", clear 
merge 1:1 label using "\Analysis\ICD_Code_Chapter_Color_No_Users.dta" 
keep if _merge==3
drop _merge
// merge ICD description 
merge 1:1 label using "\Analysis\ICD_Codes_Describtion.dta"
keep if _merge==3
drop _merge
export delimited using "\Analysis\To_Gephi\Network_Modules_With_Full_Attributes.csv", replace

// determine the ICD chapter distribtion in each module
drop _merge
preserve
keep if modularity_class == 0
tab icdgrp
restore // and so on 


// make a network filtered for no. op patients (drop if < 50 ptients), and edges (drop if < 5 edges)
import delimited "\Analysis\reults_df_RR_Phi.csv", clear
drop if patientcount < 5
rename disease1 label
merge m:1 label using "\Analysis\ICD_Code_Chapter_Color_No_Users.dta"
keep if _merge == 3
drop _merge
rename  label disease1
rename disease2 label
import delimited "\reults_df_RR_Phi.csv", clear
drop if patientcount < 5
rename disease1 label
merge m:1 label using "\Analysis\ICD_Code_Chapter_Color_No_Users.dta"
keep if _merge == 3
drop _merge
rename  label disease1
rename disease2 label
rename no_patients_each_code no_patients_each_code_2
merge m:1 label using "\Analysis\ICD_Code_Chapter_Color_No_Users.dta"
keep if _merge == 3
drop _merge
sort no_patients_each_code
drop if no_patients_each_code < 50
drop if no_patients_each_code_2 < 50
sort no_patients_each_code_2 
drop icdgrp color no_patients_each_code_2 no_patients_each_code
rename label disease2
preserve
drop if disease1 == disease2
duplicates drop
nwfromedge disease1 disease2 patientcount, undirected keeporiginal name(network_all_filtered_50_patient_5_edgeweight)
nwexport, type(pajek) replace
// export from Gephi as .csv 
// Heatmap ICD-10 chapter coocurances in the network edges

import delimited "\Analysis\reults_df_RR_Phi_Filtered_50_Patients_5_Edgeweight.csv", clear 
preserve
gen ICD1=substr( disease1 ,1,1)
gen ICD2=substr( disease1 ,1,2)
gen icdgrp=0
replace icdgrp=1 if ICD1=="A" | ICD1=="B"
replace icdgrp=2 if ICD2=="C0" | ICD2=="C1"|ICD2=="C2"| ICD2=="C3"| ICD2=="C4"| ICD2=="C5"| ICD2=="C6"| ICD2=="C7"| ICD2=="C8"| ICD2=="C9"| ICD2=="D0"|ICD2=="D1"| ICD2=="D2"| ICD2=="D3"| ICD2=="D4"
replace icdgrp=3 if ICD2=="D5"| ICD2=="D6"| ICD2=="D7"| ICD2=="D8"
replace icdgrp=4 if ICD2=="E0"| ICD2=="E1"| ICD2=="E2"| ICD2=="E3"| ICD2=="E4"| ICD2=="E5"| ICD2=="E6"| ICD2=="E7"| ICD2=="E8"| ICD2=="E9"
replace icdgrp=5 if ICD2=="F0"| ICD2=="F1"| ICD2=="F2"| ICD2=="F3"| ICD2=="F4"| ICD2=="F5"| ICD2=="F6"| ICD2=="F7"| ICD2=="F8"| ICD2=="F9"
replace icdgrp=6 if ICD2=="G0"| ICD2=="G1"| ICD2=="G2"| ICD2=="G3"| ICD2=="G4"| ICD2=="G5"| ICD2=="G6"| ICD2=="G7"| ICD2=="G8"| ICD2=="G9"
replace icdgrp=7 if ICD2=="H0"| ICD2=="H1"| ICD2=="H2"| ICD2=="H3"| ICD2=="H4"| ICD2=="H5"
replace icdgrp=8 if ICD2=="H6"| ICD2=="H7"| ICD2=="H8"| ICD2=="H9"
replace icdgrp=9 if ICD2=="I0"| ICD2=="I1"| ICD2=="I2"| ICD2=="I3"| ICD2=="I4"| ICD2=="I5"| ICD2=="I6"| ICD2=="I7"| ICD2=="I8"| ICD2=="I9"
replace icdgrp=10 if ICD2=="J0"| ICD2=="J1"| ICD2=="J2"| ICD2=="J3"| ICD2=="J4"| ICD2=="J5"| ICD2=="J6"| ICD2=="J7"| ICD2=="J8"| ICD2=="J9"
replace icdgrp=11 if ICD2=="K0"| ICD2=="K1"| ICD2=="K2"| ICD2=="K3"| ICD2=="K4"| ICD2=="K5"| ICD2=="K6"| ICD2=="K7"| ICD2=="K8"| ICD2=="K9"
replace icdgrp=12 if ICD2=="L0"| ICD2=="L1"| ICD2=="L2"| ICD2=="L3"| ICD2=="L4"| ICD2=="L5"| ICD2=="L6"| ICD2=="L7"| ICD2=="L8"| ICD2=="L9"
replace icdgrp=13 if ICD2=="M0"| ICD2=="M1"| ICD2=="M2"| ICD2=="M3"| ICD2=="M4"| ICD2=="M5"| ICD2=="M6"| ICD2=="M7"| ICD2=="M8"| ICD2=="M9"
replace icdgrp=14 if ICD2=="N0"| ICD2=="N1"| ICD2=="N2"| ICD2=="N3"| ICD2=="N4"| ICD2=="N5"| ICD2=="N6"| ICD2=="N7"| ICD2=="N8"| ICD2=="N9"
replace icdgrp=15 if ICD2=="O0"| ICD2=="O1"| ICD2=="O2"| ICD2=="O3"| ICD2=="O4"| ICD2=="O5"| ICD2=="O6"| ICD2=="O7"| ICD2=="O8"| ICD2=="O9"
replace icdgrp=16 if ICD2=="P0"| ICD2=="P1"| ICD2=="P2"| ICD2=="P3"| ICD2=="P4"| ICD2=="P5"| ICD2=="P6"| ICD2=="P7"| ICD2=="P8"| ICD2=="P9"
replace icdgrp=17 if ICD2=="Q0"| ICD2=="Q1"| ICD2=="Q2"| ICD2=="Q3"| ICD2=="Q4"| ICD2=="Q5"| ICD2=="Q6"| ICD2=="Q7"| ICD2=="Q8"| ICD2=="Q9"
replace icdgrp=18 if ICD2=="R"| ICD2=="R1"| ICD2=="R2"| ICD2=="R3"| ICD2=="R4"| ICD2=="R5"| ICD2=="R6"| ICD2=="R7"| ICD2=="R8"| ICD2=="R9"
replace icdgrp=19 if ICD2=="S0"| ICD2=="S1"| ICD2=="S2"| ICD2=="S3"| ICD2=="S4"| ICD2=="S5"| ICD2=="S6"| ICD2=="S7"| ICD2=="S8"| ICD2=="S9"|ICD2=="T0"| ICD2=="T1"| ICD2=="T2"| ICD2=="T3"| ICD2=="T4"| ICD2=="T5"| ICD2=="T6"| ICD2=="T7"| ICD2=="T8"| ICD2=="T9"
drop ICD1 ICD2
rename icdgrp icdgrp1
gen ICD1=substr( disease2 ,1,1)
gen ICD2=substr( disease2 ,1,2)
gen icdgrp=0
replace icdgrp=1 if ICD1=="A" | ICD1=="B"
replace icdgrp=2 if ICD2=="C0" | ICD2=="C1"|ICD2=="C2"| ICD2=="C3"| ICD2=="C4"| ICD2=="C5"| ICD2=="C6"| ICD2=="C7"| ICD2=="C8"| ICD2=="C9"| ICD2=="D0"|ICD2=="D1"| ICD2=="D2"| ICD2=="D3"| ICD2=="D4"
replace icdgrp=3 if ICD2=="D5"| ICD2=="D6"| ICD2=="D7"| ICD2=="D8"
replace icdgrp=4 if ICD2=="E0"| ICD2=="E1"| ICD2=="E2"| ICD2=="E3"| ICD2=="E4"| ICD2=="E5"| ICD2=="E6"| ICD2=="E7"| ICD2=="E8"| ICD2=="E9"
replace icdgrp=5 if ICD2=="F0"| ICD2=="F1"| ICD2=="F2"| ICD2=="F3"| ICD2=="F4"| ICD2=="F5"| ICD2=="F6"| ICD2=="F7"| ICD2=="F8"| ICD2=="F9"
replace icdgrp=6 if ICD2=="G0"| ICD2=="G1"| ICD2=="G2"| ICD2=="G3"| ICD2=="G4"| ICD2=="G5"| ICD2=="G6"| ICD2=="G7"| ICD2=="G8"| ICD2=="G9"
replace icdgrp=7 if ICD2=="H0"| ICD2=="H1"| ICD2=="H2"| ICD2=="H3"| ICD2=="H4"| ICD2=="H5"
replace icdgrp=8 if ICD2=="H6"| ICD2=="H7"| ICD2=="H8"| ICD2=="H9"
replace icdgrp=9 if ICD2=="I0"| ICD2=="I1"| ICD2=="I2"| ICD2=="I3"| ICD2=="I4"| ICD2=="I5"| ICD2=="I6"| ICD2=="I7"| ICD2=="I8"| ICD2=="I9"
replace icdgrp=10 if ICD2=="J0"| ICD2=="J1"| ICD2=="J2"| ICD2=="J3"| ICD2=="J4"| ICD2=="J5"| ICD2=="J6"| ICD2=="J7"| ICD2=="J8"| ICD2=="J9"
replace icdgrp=11 if ICD2=="K0"| ICD2=="K1"| ICD2=="K2"| ICD2=="K3"| ICD2=="K4"| ICD2=="K5"| ICD2=="K6"| ICD2=="K7"| ICD2=="K8"| ICD2=="K9"
replace icdgrp=12 if ICD2=="L0"| ICD2=="L1"| ICD2=="L2"| ICD2=="L3"| ICD2=="L4"| ICD2=="L5"| ICD2=="L6"| ICD2=="L7"| ICD2=="L8"| ICD2=="L9"
replace icdgrp=13 if ICD2=="M0"| ICD2=="M1"| ICD2=="M2"| ICD2=="M3"| ICD2=="M4"| ICD2=="M5"| ICD2=="M6"| ICD2=="M7"| ICD2=="M8"| ICD2=="M9"
replace icdgrp=14 if ICD2=="N0"| ICD2=="N1"| ICD2=="N2"| ICD2=="N3"| ICD2=="N4"| ICD2=="N5"| ICD2=="N6"| ICD2=="N7"| ICD2=="N8"| ICD2=="N9"
replace icdgrp=15 if ICD2=="O0"| ICD2=="O1"| ICD2=="O2"| ICD2=="O3"| ICD2=="O4"| ICD2=="O5"| ICD2=="O6"| ICD2=="O7"| ICD2=="O8"| ICD2=="O9"
replace icdgrp=16 if ICD2=="P0"| ICD2=="P1"| ICD2=="P2"| ICD2=="P3"| ICD2=="P4"| ICD2=="P5"| ICD2=="P6"| ICD2=="P7"| ICD2=="P8"| ICD2=="P9"
replace icdgrp=17 if ICD2=="Q0"| ICD2=="Q1"| ICD2=="Q2"| ICD2=="Q3"| ICD2=="Q4"| ICD2=="Q5"| ICD2=="Q6"| ICD2=="Q7"| ICD2=="Q8"| ICD2=="Q9"
replace icdgrp=18 if ICD2=="R"| ICD2=="R1"| ICD2=="R2"| ICD2=="R3"| ICD2=="R4"| ICD2=="R5"| ICD2=="R6"| ICD2=="R7"| ICD2=="R8"| ICD2=="R9"
replace icdgrp=19 if ICD2=="S0"| ICD2=="S1"| ICD2=="S2"| ICD2=="S3"| ICD2=="S4"| ICD2=="S5"| ICD2=="S6"| ICD2=="S7"| ICD2=="S8"| ICD2=="S9"|ICD2=="T0"| ICD2=="T1"| ICD2=="T2"| ICD2=="T3"| ICD2=="T4"| ICD2=="T5"| ICD2=="T6"| ICD2=="T7"| ICD2=="T8"| ICD2=="T9"
rename icdgrp icdgrp2
drop ICD1 ICD2
tab icdgrp1 icdgrp2 // paste the table to excel --> move to python to make teh heat map 

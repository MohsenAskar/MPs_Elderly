#############################################################################
# Load data from created in Stata to calculate RR and Phi
############################################################
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
observed_cooccurrences = data_pivot.T.dot(data_pivot).values  

# Convert observed_cooccurrences to a DataFrame 
disease_names = data_pivot.columns
observed_df = pd.DataFrame(observed_cooccurrences, index=disease_names, columns=disease_names)

# Filter out disease pairs with co-occurrences in only one patient to avoid biad of RR and Phi towards rare conditions
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

# Compile results including patient count for each filtered pair
result = []
for i in range(len(disease_names)):
    for j in range(i+1, len(disease_names)):
        if not np.isnan(rr_filtered[i,j]) and not np.isnan(phi_filtered[i,j]):
            count = filtered_observed_df.iloc[i, j]  
            if not np.isnan(count):  
                result.append((disease_names[i], disease_names[j], rr_filtered[i,j], phi_filtered[i,j], count))

result_df = pd.DataFrame(result, columns=["Disease1", "Disease2", "RR", "Phi", "PatientCount"])

# Display the resulting DataFrame with RR and Phi values for disease pairs
print(result_df)
result_df['PatientCount'].describe()
result_df['Disease1'].nunique() #641
result_df['Disease2'].nunique()#682
result_df.to_csv(r"\Analysis\reults_df_RR_Phi.csv", index= False)

# counting uniques patiets and ICD codes (for the exculsion figure)
unique_diseases_result_df = pd.unique(result_df[['Disease1', 'Disease2']].values.ravel()).size
# 1. Extract disease pairs from 'result_df'
disease_pairs = [(row['Disease1'], row['Disease2']) for index, row in result_df.iterrows()]

# 2. Map disease pairs to patients
unique_patients = set()  # To store unique patient identifiers

for disease1, disease2 in disease_pairs:
    # Find patients with both diseases
    patients_with_pair = data_pivot[(data_pivot[disease1] > 0) & (data_pivot[disease2] > 0)].index
    # Add these patients to the set of unique patients
    unique_patients.update(patients_with_pair)

# 3. Count unique patients
unique_patient_count = len(unique_patients)

# to calculate the uniqiue patients and ICD codes for the exclusion flow diagram
###################################################################################
# Apply the additional exclusion criteria
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

# Confirm: count unique patients and ICD codes after final filtering
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

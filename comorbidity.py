# import pandas as pd
# from prefixspan import PrefixSpan

# # Load the Stata dataset into a pandas DataFrame
# df = pd.read_stata(r'C:\Users\mas082\OneDrive - UiT Office 365\Desktop\VS_Code\Comorbidity_Patterns\ICD_sequences.dta')

# # Extract the sequence variables into a list of sequences
# sequences = []
# for i in range(2, 246):
#     colname = 'hovedtilstand{}'.format(i)
#     sequence = df[colname].tolist()
#     sequence = [x for x in sequence if pd.notnull(x)]
#     if sequence:
#         sequences.append(sequence)

# # Set the minimum support threshold for frequent pattern mining
# min_support = 3

# # Instantiate the PrefixSpan algorithm with the minimum support threshold
# ps = PrefixSpan(sequences)
# ps.minlen = 3  # set the minimum pattern length to 2

# # Mine frequent patterns using PrefixSpan
# freq_patterns = ps.frequent(min_support)

# # Print the frequent patterns
# for pattern, support in freq_patterns:
#     print(pattern, support)

# #############################################
# # calculating the probailities of sequences 
# #############################################

# # Create a dictionary to store the counts of each diagnosis transition
# transition_counts = {}

# # Loop over each sequence in the DataFrame
# for index, row in df.iterrows():
#     sequence = row.filter(regex='hovedtilstand').dropna().tolist()
#     # Loop over each diagnosis in the sequence
#     for i in range(len(sequence)-1):
#         from_diagnosis = sequence[i]
#         to_diagnosis = sequence[i+1]
#         # Increment the count for this diagnosis transition in the dictionary
#         if from_diagnosis in transition_counts:
#             if to_diagnosis in transition_counts[from_diagnosis]:
#                 transition_counts[from_diagnosis][to_diagnosis] += 1
#             else:
#                 transition_counts[from_diagnosis][to_diagnosis] = 1
#         else:
#             transition_counts[from_diagnosis] = {to_diagnosis: 1}

# # Create a DataFrame to store the transition probabilities
# transition_probabilities = pd.DataFrame(columns=df.filter(regex='hovedtilstand').columns)

# # Loop over each diagnosis and calculate the probability of each transition
# for from_diagnosis in transition_counts.keys():
#     total_transitions = sum(transition_counts[from_diagnosis].values())
#     for to_diagnosis in transition_counts[from_diagnosis].keys():
#         transition_probabilities.loc[from_diagnosis, to_diagnosis] = transition_counts[from_diagnosis][to_diagnosis] / total_transitions

# # Print the transition probabilities DataFrame
# print(transition_probabilities)

# ######################################################################################################################################################
# # In case we want to calucalte the RR and Φ-correlation between diseases in the lontidunal format of the dataset. 
# # RR > 1 and Φ-correlation > 0 will indicate that
# # the disease co-occurancies are not by chance
# #################################################################################################################################################################

# import pandas as pd
# import numpy as np
# from scipy.stats import chi2_contingency

# # read in your longitudinal hospital dataset
# data = pd.read_csv(r"D:\Tauseef\for_report\only_patient_ICD.csv")
# data.columns
# data = pd.read_csv(r"C:\Users\mas082\OneDrive - UiT Office 365\Desktop\test_after_CCIR.csv")

# data = data.rename(columns={'pasientlopenr': 'patientid', 'hovedtilstand': 'code_system'})
# data.shape
# unique_values = data['code_system'].nunique() #1289

# # remove < 1% prevelant diseases 
# #################################
# # Identify unique occurrences of ICD codes per patient
# unique_icd_per_patient = data.drop_duplicates(subset=['patientid', 'code_system'])

# # Calculate the total number of unique patients
# total_patients = unique_icd_per_patient['patientid'].nunique()

# # Calculate the prevalence of each ICD code
# icd_prevalence = unique_icd_per_patient['code_system'].value_counts() / total_patients

# icd_prevalence_df = icd_prevalence.reset_index()
# icd_prevalence_df.columns = ['code_system', 'prevalence']

# # save the calucalted prevelance as .csv
# icd_prevalence_df.to_csv(r'C:\Users\mas082\OneDrive - UiT Office 365\Desktop\check_prevelances.csv', index=False)
# # Filter out ICD codes with a prevalence of less than 1%
# icd_prevalence_filtered = icd_prevalence[icd_prevalence >= 0.01]
# # Display the filtered ICD codes and their prevalence
# icd_prevalence_filtered
# # Correct calculation for the percentage of codes removed
# percent_codes_removed = (1 - len(icd_prevalence_filtered) / len(icd_prevalence)) * 100
# percent_codes_removed

# # Filter the original dataset to include only rows with ICD codes that are >= 1% prevalent
# data_filtered = data[~data['code_system'].isin(icd_prevalence_filtered.index)]

# # Save the filtered dataset to a CSV file
# data_filtered.to_csv(r'C:\Users\mas082\OneDrive - UiT Office 365\Desktop\data_after_prevelance.csv', index=False)
# # In Stata, we filtered out the ICD codes that are not chronic using CCIR

# calculate RR and phi
##########################
# data= pd.read_csv(r"C:\Users\mas082\OneDrive - UiT Office 365\Desktop\Only_Chronic_after_CCIR.csv")
# data.shape
# data = data.rename(columns={'pasientlopenr': 'patientid', 'hovedtilstand': 'code_system'})

# # Assuming your data is already loaded and pre-processed
# data["Count"] = 1
# data_pivot = data.pivot_table(index="patientid", columns="code_system", values="Count", aggfunc="max").fillna(0).astype(int)

# # Calculate incidence rates
# incidence_rates = data_pivot.sum() / data_pivot.shape[0]

# # Expected and observed rates
# expected_rates = np.outer(incidence_rates, incidence_rates)
# observed_rates = (data_pivot.T.dot(data_pivot) / data_pivot.shape[0]).values
# observed_cooccurrences = data_pivot.T.dot(data_pivot).values  # Actual count of co-occurrences

# disease_names = data_pivot.columns
# observed_df = pd.DataFrame(observed_cooccurrences, index=disease_names, columns=disease_names)
# filtered_observed_df = observed_df.where(observed_df > 1)

# # RR, Phi calculations
# rr_values = observed_rates / expected_rates
# chi2_values = (observed_rates * data_pivot.shape[0] - expected_rates * data_pivot.shape[0])**2 / (expected_rates * data_pivot.shape[0])
# phi_values = np.sqrt(chi2_values / data_pivot.shape[0])
# rr_filtered = np.where(rr_values > 1, rr_values, np.nan)
# phi_filtered = np.where(phi_values > 0, phi_values, np.nan)

# # Compile results, including patient count for each pair
# result = []
# disease_names = data_pivot.columns


# for i in range(len(disease_names)):
#     for j in range(i+1, len(disease_names)):
#         if not np.isnan(rr_filtered[i,j]) and not np.isnan(phi_filtered[i,j]):
#             count = observed_cooccurrences[i, j]  # The count of patients with both diseases
#             result.append((disease_names[i], disease_names[j], rr_filtered[i,j], phi_filtered[i,j], count))

# result_df = pd.DataFrame(result, columns=["Disease1", "Disease2", "RR", "Phi", "PatientCount"])
# print(result_df)
# result_df['Disease1'].nunique()
# result_df['Disease2'].nunique()

# result_df.to_csv(r"C:\Users\mas082\OneDrive - UiT Office 365\Desktop\reults_df_RR_Phi.csv", index= False)


#############################################################################
# Start of code of the articel
#############################################################################
# Load data from created in Stata to calculate RR and Phi
############################################################
# filter out disease pairs that occurs in only one patient and then calulate RR and Phi
#######################################################################################
import pandas as pd
import numpy as np
# Loading data created in Stata
data= pd.read_csv(r"C:\Users\mas082\OneDrive - UiT Office 365\Desktop\Ongoing_Projects\Modularity_Comorbidy_Detection\Analysis\Only_Chronic_after_CCIR.csv")
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
result_df.to_csv(r"C:\Users\mas082\OneDrive - UiT Office 365\Desktop\Ongoing_Projects\Modularity_Comorbidy_Detection\Analysis\reults_df_RR_Phi.csv", index= False)

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
final_result_df.to_csv(r"C:\Users\mas082\OneDrive - UiT Office 365\Desktop\Ongoing_Projects\Modularity_Comorbidy_Detection\Analysis\final_results_df_RR_Phi.csv", index=False)

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
#############################################################################
# End of code of the articel, what comes next is experiments
#############################################################################
# making a network from it using Disease 1, 2 as the node pairs and Patietcount as the edge thickness
import networkx as nx
import matplotlib.pyplot as plt

# Initialize a new graph
G = nx.Graph()

# Add edges and nodes to the graph from the DataFrame
for index, row in result_df.iterrows():
    G.add_edge(row['Disease1'], row['Disease2'], weight=row['PatientCount'])

total_nodes = G.number_of_nodes()
total_edges = G.number_of_edges()

print(f"Total number of nodes: {total_nodes}")
print(f"Total number of edges: {total_edges}")

# export to Gephi as pajak
nx.write_pajek(G, r"C:\Users\mas082\OneDrive - UiT Office 365\Desktop\network_multimorbidity.net")



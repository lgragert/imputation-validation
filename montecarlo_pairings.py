
import pandas as pd
import json
import requests
import time

# Create a Monte Carlo Sampling of the pairs rather than going through each of them in the API

# Clean the JSON string to have it as 'eplet_eplet_eplet'
def clean_eplet_str(eplet_list_df):
    eplet_list_df = eplet_list_df.str.replace(', ', '_')
    eplet_list_df = eplet_list_df.str.replace('\'', '')
    eplet_list_df = eplet_list_df.str.replace('[', '', regex=True).replace(']', '', regex=True)
    return eplet_list_df


# Select random pairings from the truth table
truth_filename = 'DRDQ_pairs_truth.csv'
truth_pairs = pd.read_csv(truth_filename, header=0)

# Start with 100 random pairs as that is what the API can handle: 100, 200, 500, 1,000... until API breaks
n_pairs = 100
n_truthpairs = truth_pairs.sample(n=n_pairs, random_state=2272024).reset_index(drop=True)  # Initialize state, to keep it consistent

n_truth_file = 'DRDQ_pairs_truth_' + str(n_pairs) + ".csv"
n_truthpairs.to_csv(n_truth_file, header=True, index=False)

# Create a dictionary of pairings to find them in the imputation file
rand_pairs = {}
for row in range(len(n_truthpairs)):
    don_id = n_truthpairs.loc[row, 'DON_ID']
    rec_id = n_truthpairs.loc[row, 'REC_ID']

    pairing = don_id + "+" + rec_id

    if pairing not in rand_pairs:
        rand_pairs[pairing] = 1
    else:
        rand_pairs[pairing] = rand_pairs[pairing] + 1  # Lets us know if it chose more than one of the same pairing

eplet_truth = pd.DataFrame()
for line in range(len(n_truthpairs)):
    donor_ID = n_truthpairs.loc[line, 'DON_ID']
    recip_ID = n_truthpairs.loc[line, 'REC_ID']
    donor_geno_DRDQ = n_truthpairs.loc[line, 'DON_GLString']
    recip_geno_DRDQ = n_truthpairs.loc[line, 'REC_GLString']
    id_pair = donor_ID + "+" + recip_ID

    if len(eplet_truth) % 100 == 0:  # can only handle 100 patients at a time, so let it take a break
        time.sleep(60.0)

    immunizer_alleles = donor_geno_DRDQ
    patient_alleles = recip_geno_DRDQ

    # Use the HLA eplet registry web service, compute eplet MM for a transplant pair
    request_path = 'https://api.epregistry.com.br/eplet_mismatches?from=jh35hk423khbxin3e9849r8fuicne9ne&immunizer_alleles=' + immunizer_alleles + '&patient_alleles=' + patient_alleles
    r = requests.get(request_path, headers={'Accept': 'application/json'})
    out_put = r.json()

    # Break up the JSON output and clean it
    e_df = pd.DataFrame()
    for key, value in out_put.items():
        eplet_dict = {}
        if key == 'ABC' or key == 'DP' or key == 'MICA' or key == 'version':  # Only looking at DRDQ, version cannot do the loop, at this time version=2024-01-22
            continue
        for labels, extended in value.items():
            named = key + "_" + labels
            eplet_dict[named] = str(extended)
            eplet_df = pd.DataFrame(eplet_dict, index=[id_pair], columns=[named])
            if eplet_df.columns == key + "_" + 'details':
                eplet_df[key + "_" + 'details'] = clean_eplet_str(eplet_df[key + "_" + 'details'])
            e_df = pd.concat([e_df, eplet_df], axis=1)
    eplet_truth = pd.concat([eplet_truth, e_df])
    time.sleep(1.0)

# Clean the format a little bit to only the columns you want
eplet_truth_table = eplet_truth[['ALL_quantity', 'ALL_details', 'DRB_quantity', 'DRB_details', 'DQ_quantity', 'DQ_details']].reset_index(names=['ID'])
print(eplet_truth_table)
eplet_truth_table.to_csv('eplet_truth_table' + str(n_pairs) + '.csv', index=False, header=True)


# Run through all possibilities for those pairings in the imputation file
impute_file = 'DRDQ_pairs_imputation.csv'
impute_pairs = pd.read_csv(impute_file, header=0)

n_impute_pairs = pd.DataFrame()
for pair in rand_pairs:
    don, rec = pair.split('+')
    line = impute_pairs[(impute_pairs['DON_ID'] == don) & (impute_pairs['REC_ID'] == rec)]
    # In case the donor and recipient are switched
    if line.empty is True:
        line = impute_pairs[(impute_pairs['DON_ID'] == rec) & (impute_pairs['REC_ID'] == don)]

    line = line.sort_values(by=['PairProb'], ascending=False)
    line['Rank'] = range(1, len(line) + 1)

    prob_sum = line['PairProb'].sum()
    print(pair + ": " + str(prob_sum))  # Check each pair's probability == 1

    n_impute_pairs = pd.concat([n_impute_pairs, line])

n_impute_pairs = n_impute_pairs.reset_index(drop=True)
n_impute_filename = 'DRDQ_pairs_imputation_' + str(n_pairs) + ".csv"
n_impute_pairs.to_csv(n_impute_filename, header=True, index=False)

# Put the imputation through the eplet API
eplet_dataframe = pd.DataFrame()
for line in range(len(n_impute_pairs)):
    donor_ID = n_impute_pairs.loc[line, 'DON_ID']
    recip_ID = n_impute_pairs.loc[line, 'REC_ID']
    donor_geno_DRDQ = n_impute_pairs.loc[line, 'DON_GLString']
    recip_geno_DRDQ = n_impute_pairs.loc[line, 'REC_GLString']
    pair_prob = n_impute_pairs.loc[line, 'PairProb']
    id_pair = donor_ID + "+" + recip_ID

    # Can only handle a 100 pairs at a time, so make sure to take a break
    if len(eplet_dataframe) % 100 == 0:
        time.sleep(60.0)

    # Need to know which genotypes are used for each donor-recipient pairing corresponding to the IDs
    donor_line = pd.DataFrame({'DON_ID': donor_ID, 'DON_GLString': donor_geno_DRDQ, 'REC_ID': recip_ID, 'REC_GLString': recip_geno_DRDQ, "PairProb": pair_prob}, index=[id_pair])

    immunizer_alleles = donor_geno_DRDQ
    patient_alleles = recip_geno_DRDQ

    # Use the HLA eplet registry web service, compute eplet MM for a transplant pair
    request_path = 'https://api.epregistry.com.br/eplet_mismatches?from=jh35hk423khbxin3e9849r8fuicne9ne&immunizer_alleles=' + immunizer_alleles + '&patient_alleles=' + patient_alleles
    r = requests.get(request_path, headers={'Accept': 'application/json'})
    out_put = r.json()

    e_df = pd.DataFrame()
    for key, value in out_put.items():
        eplet_dict = {}
        if key == 'ABC' or key == 'DP' or key == 'MICA' or key == 'version':  # Only looking at DRDQ, version cannot do the loop, at this time version=2024-01-22
            continue
        for labels, extended in value.items():
            named = key + "_" + labels
            eplet_dict[named] = str(extended)
            eplet_df = pd.DataFrame(eplet_dict, index=[id_pair], columns=[named])
            if eplet_df.columns == key + "_" + 'details':
                eplet_df[key + "_" + 'details'] = clean_eplet_str(eplet_df[key + "_" + 'details'])
            e_df = pd.concat([e_df, eplet_df], axis=1)
    donor_line = pd.concat([donor_line, e_df], axis=1)  # concat the donor-recipient genotypes to clean eplet line
    eplet_dataframe = pd.concat([eplet_dataframe, donor_line])
    print(eplet_dataframe)
    time.sleep(1.0)

eplet_impute = eplet_dataframe[['DON_ID', 'REC_ID', 'PairProb', 'DON_GLString', 'REC_GLString', 'ALL_quantity', 'ALL_details', 'DRB_quantity', 'DRB_details', 'DQ_quantity', 'DQ_details']].reset_index(drop=True)  # Only need these columns for now
print(eplet_impute)
eplet_impute.to_csv('eplet_lowres_impute' + str(n_pairs) + '.csv', header=True, index=False)

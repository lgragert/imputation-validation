
import pandas as pd
from collections import defaultdict
import json
import requests
import time

# Compute eplet MM for HLA DR-DQ (DRB345, DRB1, DQA1, DQB1) in a simulated dataset of 216 'donor-recipients'


# Separate GLString for truth tables for DR-DQ genotypes
def sep_glstring(file):
    file[['A', 'B', 'C', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']] = file['GLString'].str.split('^', expand=True)
    file = file.drop(columns=['GLString', 'A', 'B', 'C', 'DPA1', 'DPB1'])

    loci = ['DRB345', 'DRB1', 'DQA1', 'DQB1']
    for locus in loci:
        file[[locus + '_1', locus + '_2']] = file[locus].str.split('+', expand=True)
        file = file.drop(columns=[locus])

    return file


# Clean the JSON string to have it as 'eplet_eplet_eplet'
def clean_eplet_str(eplet_list_df):
    eplet_list_df = eplet_list_df.str.replace(', ', '_')
    eplet_list_df = eplet_list_df.str.replace('\'', '')
    eplet_list_df = eplet_list_df.str.replace('[', '', regex=True).replace(']', '', regex=True)
    return eplet_list_df


# Clean truth table to only have DR-DQ genotypes and the simulated donor-recipient pairs
truth_file = 'genotype_truth_table.csv'
truth_table = pd.read_csv(truth_file, header=0)
truth_table = sep_glstring(truth_table)
# truth_table.to_csv('DRDQ_genotype_truth_table.csv', header=False, index=False)  # not in GLString format

truth_table['GLString'] = (truth_table['DRB345_1'] + "+" + truth_table['DRB345_2'] +
                           '^' + truth_table['DRB1_1'] + "+" + truth_table['DRB1_2'] +
                           '^' + truth_table['DQA1_1'] + "+" + truth_table['DQA1_2'] +
                           '^' + truth_table['DQB1_1'] + "+" + truth_table['DQB1_2'])

# Create truth dictionary with {ID: genotype: 1}
TT_DRDQ = defaultdict(dict)
for line in range(len(truth_table)):
    subject_id = truth_table.loc[line, 'ID']
    GENO_DRDQ = truth_table.loc[line, 'GLString']

    if GENO_DRDQ not in TT_DRDQ:
        TT_DRDQ[subject_id][GENO_DRDQ] = 1
    else:
        continue


# Simulate highres "donor-recipient" pairs from the truth dataset by looping through all combinations
truth_pairs = {}
eplet_truth = pd.DataFrame()
for donor_ID in TT_DRDQ:
    for recip_ID in TT_DRDQ:
        if donor_ID == recip_ID:  # skip pairs with same ID
            continue

        # skip pairs where ID was already computed
        (donor_ID, recip_ID) = sorted([donor_ID, recip_ID])
        id_pair = "+".join([donor_ID, recip_ID])
        if id_pair in truth_pairs:
            continue
        else:
            truth_pairs[id_pair] = 1

        # loop through all possible genotypes for donor and recip pairing
        for donor_geno_DRDQ in TT_DRDQ[donor_ID]:
            for recip_geno_DRDQ in TT_DRDQ[recip_ID]:
                if len(eplet_truth) % 100 == 0:  # can only handle 100 patients at a time, so let it take a break
                    time.sleep(60.0)

                # Format donor and recip genotype for input to web service, which is a comma separated list
                immunizer_alleles = donor_geno_DRDQ
                immunizer_alleles = immunizer_alleles.replace('^', ',').replace('+', ',')
                patient_alleles = recip_geno_DRDQ
                patient_alleles = patient_alleles.replace('^', ',').replace('+', ',')

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
                time.sleep(2.0)

# Clean the format a little bit to only the columns you want
eplet_truth_table = eplet_truth[['ALL_quantity', 'ALL_details', 'DRB_quantity', 'DRB_details', 'DQ_quantity', 'DQ_details']].reset_index(names=['ID'])
print(eplet_truth_table)
eplet_truth_table.to_csv('eplet_truth_table.csv', index=False, header=True)  # work with eplets as a CSV and not worry about the JSON formatting


# Simulate imputation "donor-recipient" pairs from imputation dataset by looping through all combinations
impute_file = 'lowres_DRDQ_impute.csv'
impute_lowres = pd.read_csv(impute_file, header=0)

# Create a dictionary of the dataset such that GF_DRDQ[subject_ID][GENO_DRDQ][freq]
GF_DRDQ = defaultdict(dict)
for line in range(len(impute_lowres)):
    subject_id = impute_lowres.loc[line, 'ID']
    GENO_DRDQ = impute_lowres.loc[line, 'DRDQ_GLString']
    freq = impute_lowres.loc[line, 'DQDR_freq']

    if GENO_DRDQ not in GF_DRDQ[subject_id]:
        GF_DRDQ[subject_id][GENO_DRDQ] = freq
    else:
        GF_DRDQ[subject_id][GENO_DRDQ] = GF_DRDQ[subject_id][GENO_DRDQ] + freq


# Use imputation output, compute eplet MM for each possible donor and recipient and store it into a DataFrame
pair_computed = {}
eplet_dataframe = pd.DataFrame()
for donor_ID in GF_DRDQ:
    for recip_ID in GF_DRDQ:
        if donor_ID == recip_ID:  # skip pairs with same ID
            continue

        # skip pairs where ID was already computed
        (donor_ID, recip_ID) = sorted([donor_ID, recip_ID])
        id_pair = "+".join([donor_ID, recip_ID])
        if id_pair in pair_computed:
            continue
        else:
            pair_computed[id_pair] = 1

        # loop through all possible genotypes for donor and recip pairing
        for donor_geno_DRDQ in GF_DRDQ[donor_ID]:
            for recip_geno_DRDQ in GF_DRDQ[recip_ID]:
                # Can only handle a 100 pairs at a time, so make sure to take a break
                if len(eplet_dataframe) % 100 == 0:
                    time.sleep(60.0)

                # Need to know which genotypes are used for each donor-recipient pairing corresponding to the IDs
                donor_line = pd.DataFrame({"Donor_GLString": donor_geno_DRDQ, 'Recip_GLString': recip_geno_DRDQ}, index=[id_pair])
                # donor_line = pd.DataFrame({'DON_ID': donor_ID, 'DON_GLString': donor_geno_DRDQ, 'REC_ID': recip_ID, 'Recip_GLString': recip_geno_DRDQ}, index=[id_pair])

                # Format donor and recip genotype for input to web service, which is a comma separated list
                immunizer_alleles = donor_geno_DRDQ
                immunizer_alleles = immunizer_alleles.replace('^', ',').replace('+', ',')
                patient_alleles = recip_geno_DRDQ
                patient_alleles = patient_alleles.replace('^', ',').replace('+', ',')

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
                        e_df = pd.concat([e_df, donor_line], axis=1)  # concat the donor-recipient genotypes to clean eplet line
                eplet_dataframe = pd.concat([eplet_dataframe, e_df])
                print(eplet_dataframe)
                time.sleep(2.0)

eplet_impute = eplet_dataframe[['ALL_quantity', 'ALL_details', 'DRB_quantity', 'DRB_details', 'DQ_quantity', 'DQ_details', 'Donor_GLString', 'Recip_GLString']].reset_index(names=['ID'])  # Only need these columns for now
print(eplet_impute)
eplet_impute.to_csv('eplet_lowres_impute.csv', header=True, index=False)  # work with eplets as a CSV and not worry about the JSON formatting

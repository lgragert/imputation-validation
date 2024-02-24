
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
impute_file = 'lowres_DRDQ_impute.csv'
impute_lowres = pd.read_csv(impute_file, header=0)
truth_table = truth_table[truth_table.ID.isin(impute_lowres.ID)].reset_index(drop=True)
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
        #(donor_ID, recip_ID) = sorted([donor_ID, recip_ID])
        id_pair = "+".join([donor_ID, recip_ID])
        if "+".join([recip_ID,donor_ID]) not in truth_pairs:
            truth_pairs[id_pair] = 1

            # loop through all possible genotypes for donor and recip pairing
            for donor_geno_DRDQ in TT_DRDQ[donor_ID]:
                for recip_geno_DRDQ in TT_DRDQ[recip_ID]:

                    immunizer_alleles = donor_geno_DRDQ
                    immunizer_alleles = immunizer_alleles.replace('^', ',').replace('+', ',')
                    patient_alleles = recip_geno_DRDQ
                    patient_alleles = patient_alleles.replace('^', ',').replace('+', ',')

                    donrec_line = pd.DataFrame({'DON_ID': donor_ID, 'DON_GLString': immunizer_alleles, 'REC_ID': recip_ID,
                                                'REC_GLString': patient_alleles}, index=[id_pair])
                    eplet_truth = pd.concat([eplet_truth, donrec_line])
            print(eplet_truth)

# Clean the format a little bit to only the columns you want
eplet_truth_table = eplet_truth.reset_index(drop=True)
print(eplet_truth_table)
eplet_truth_table.to_csv('DRDQ_pairs_truth.csv', index=False, header=True)  # work with eplets as a CSV and not worry about the JSON formatting

impute_file = 'lowres_DRDQ_impute.csv'
impute_lowres = pd.read_csv(impute_file, header=0)

# Create a dictionary of the dataset such that GF_DRDQ[subject_ID][GENO_DRDQ][freq]
GF_DRDQ = defaultdict(dict)
for line in range(len(impute_lowres)):
    subject_id = impute_lowres.loc[line, 'ID']
    GENO_DRDQ = impute_lowres.loc[line, 'DRDQ_GLString']
    freq = impute_lowres.loc[line, 'DRDQ_freq']

    if GENO_DRDQ not in GF_DRDQ[subject_id]:
        GF_DRDQ[subject_id][GENO_DRDQ] = freq
    else:
        GF_DRDQ[subject_id][GENO_DRDQ] = GF_DRDQ[subject_id][GENO_DRDQ] + freq


# Use imputation output, compute eplet MM for each possible donor and recipient and store it into a DataFrame
pair_computed = {}
eplet_prob = pd.DataFrame()
for donor_ID in GF_DRDQ:
    for recip_ID in GF_DRDQ:
        if donor_ID == recip_ID:  # skip pairs with same ID
            continue

        # skip pairs where ID was already computed
        # (donor_ID, recip_ID) = sorted([donor_ID, recip_ID])
        id_pair = "+".join([donor_ID, recip_ID])
        if "+".join([recip_ID, donor_ID]) not in pair_computed:
            pair_computed[id_pair] = 1

            # loop through all possible genotypes for donor and recip pairing
            for donor_geno_DRDQ in GF_DRDQ[donor_ID]:
                for recip_geno_DRDQ in GF_DRDQ[recip_ID]:
                    immunizer_alleles = donor_geno_DRDQ
                    immunizer_alleles = immunizer_alleles.replace('^', ',').replace('+', ',')
                    patient_alleles = recip_geno_DRDQ
                    patient_alleles = patient_alleles.replace('^', ',').replace('+', ',')

                    donrec_line = pd.DataFrame({'DON_ID': donor_ID, 'DON_GLString': immunizer_alleles, 'REC_ID': recip_ID, 'REC_GLString': patient_alleles}, index=[id_pair])
                    eplet_prob = pd.concat([eplet_prob, donrec_line])

eplet_impute = eplet_prob.reset_index(drop=True)  # Only need these columns for now
print(eplet_impute)
eplet_impute.to_csv('DRDQ_pairs_imputation.csv', header=True, index=False)

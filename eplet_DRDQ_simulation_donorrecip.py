
import pandas as pd
from collections import defaultdict
import json
import requests
import gzip
from datetime import datetime
import timeit
import itertools
import time

# Compute eplet MM for HLA DR-DQ (DRB345, DRB1, DQA1, DQB1) in a simulated dataset of ~200 'donor-recipients'

# Separate GLString for truth tables for DR-DQ genotypes
def sep_glstring(file):
    file[['A', 'B', 'C', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']] = file['GLString'].str.split('^', expand=True)
    file = file.drop(columns=['GLString', 'A', 'B', 'C', 'DPA1', 'DPB1'])

    loci = ['DRB345', 'DRB1', 'DQA1', 'DQB1']
    for locus in loci:
        file[[locus + '_1', locus + '_2']] = file[locus].str.split('+', expand=True)
        file = file.drop(columns=[locus])

    return file


# Clean truth table to only have DR-DQ genotypes and create a CSV so the format is correct for the eplet web registry
truth_file = 'genotype_truth_table.csv'
truth_table = pd.read_csv(truth_file, header=0)
truth_table = sep_glstring(truth_table)
truth_table.to_csv('DRDQ_genotype_truth_table.csv', header=False, index=False)

# Simulate highres "donor-recipient" pairs from the imputation dataset by looping through all combinations
# Compute probabilities for each possible donor-recipient pair
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


# Use imputation output, compute eplet MM for each possible donor and recipient
# Go through dictionary to find the best pairings
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
                # Format donor and recip genotype for input to web service, which is a comma separated list
                immunizer_alleles = donor_geno_DRDQ
                immunizer_alleles = immunizer_alleles.replace('^', ',')
                immunizer_alleles = immunizer_alleles.replace('+', ',')
                patient_alleles = recip_geno_DRDQ
                patient_alleles = patient_alleles.replace('^', ',')
                patient_alleles = patient_alleles.replace('+', ',')

                # Use the HLA eplet registry web service, compute eplet MM for a transplant pair
                request_path = 'https://api.epregistry.com.br/eplet_mismatches?from=jh35hk423khbxin3e9849r8fuicne9ne&immunizer_alleles=' + immunizer_alleles + '&patient_alleles=' + patient_alleles
                r = requests.get(request_path, headers={'Accept': 'application/json'})
                out_put = r.json()

                e_df = pd.DataFrame()
                for key, value in out_put.items():
                    eplet_dict = {}
                    if key == 'ABC' or key == 'DP' or key == 'MICA' or key == 'version':  # Only looking at DRDQ, version cannot do the loop, at this time version=01-2024
                        continue
                    for labels, extended in value.items():
                        named = key + "_" + labels
                        eplet_dict[named] = str(extended)
                        eplet_df = pd.DataFrame(eplet_dict, index=[id_pair], columns=[named])
                        e_df = pd.concat([e_df, eplet_df], axis=1)
                        print(e_df)
                eplet_dataframe = pd.concat([eplet_dataframe, e_df])

                # Combining the calculator output and probabilities for each possible donor and recipient, compute probability distribution for each unique eplet MM
                # TODO - Use immunizer_geno_prob and patient_geno_prob to compute eplet MM probabilities
                immunizer_geno_prob = GF_DRDQ[donor_ID][donor_geno_DRDQ]
                patient_geno_prob = GF_DRDQ[recip_ID][recip_geno_DRDQ]

eplet_dataframe.to_csv('EpRegistry_lowres_impute.csv', header=True, index=True)

# Use the true data, compute eplet MM for the actual donor and recip, using the output to create truth tables


# Combining the calculator output and probabilities for each possible donor and recipient, compute probability distribution for each unique eplet MM
# Unique set of eplet MM for: DR, DQ, DR+DQ
# Eplet MM Counts for: DR, DQ, DR+DQ
# Single molecule eplet MM risk categories: high, med, low


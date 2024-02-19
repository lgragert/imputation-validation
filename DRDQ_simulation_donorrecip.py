
import pandas as pd
from collections import defaultdict
import gzip

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


# Use the HLA eplet registry web service, compute eplet MM for a transplant pair
# Use imputation output, compute eplet MM for each possible donor and recipient
# Use the true data, compute eplet MM for the actual donor and recip, using the output to create truth tables


# Combining the calculator output and probabilities for each possible donor and recipient, compute probability distribution for each unique eplet MM
# Unique set of eplet MM for: DR, DQ, DR+DQ
# Eplet MM Counts for: DR, DQ, DR+DQ
# Single molecule eplet MM risk categories: high, med, low



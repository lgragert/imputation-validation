
import pandas as pd
from collections import defaultdict
import math
from itertools import combinations
import sys
# Create pairings for the dataset of 216 subjects.


# Separate GLString for truth tables for DR-DQ genotypes
def sep_glstring(file):
    file[['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']] = file['GLString'].str.split('^',
                                                                                                         expand=True)
    file = file.drop(columns=['GLString', 'A', 'C', 'B', 'DPA1', 'DPB1'])

    loci = ['DRB345', 'DRB1', 'DQA1', 'DQB1']
    for locus in loci:
        file[[locus + '_1', locus + '_2']] = file[locus].str.split('+', expand=True)
        file = file.drop(columns=[locus])

    return file


# Makes a dictionary turn into a DataFrame
def geno_pairs(GF_dict, don_geno, rec_geno, donID, recID, pair_ids, whichimpute):
    # Reformats GLString to API format
    don_API = don_geno
    don_API = don_API.replace('^', ',').replace('+', ',')
    rec_API = rec_geno
    rec_API = rec_API.replace('^', ',').replace('+', ',')

    # Create the pair frequency
    don_freq = GF_dict[donID][don_geno]
    rec_freq = GF_dict[recID][rec_geno]
    pair_freq = don_freq * rec_freq

    line_donrec = pd.DataFrame({'DON_ID': donID, 'DON_' + whichimpute: don_API, 'REC_ID': recID, 'REC_' + whichimpute: rec_API,
                                'PairProb_' + whichimpute: pair_freq}, index=[pair_ids])
    return line_donrec


# Clean truth table to only have DR-DQ genotypes and the simulated donor-recipient pairs
which_impute = sys.argv[3]   # will change this parameter from DRDQ, DR, or DQ
truth_file = sys.argv[1]
truth_table = pd.read_csv(truth_file, header=0, dtype={"ID": str})
impute_file = sys.argv[2]
impute_lowres = pd.read_csv(impute_file, header=0, dtype={"ID": str})
truth_table = truth_table[truth_table.ID.isin(impute_lowres.ID)].reset_index(drop=True)
truth_table = sep_glstring(truth_table)
# truth_table.to_csv('DRDQ_genotype_truth_table.csv', header=False, index=False)  # not in GLString format

truth_table['DRDQ_GLString'] = (truth_table['DRB345_1'] + "+" + truth_table['DRB345_2'] +
                                '^' + truth_table['DRB1_1'] + "+" + truth_table['DRB1_2'] +
                                '^' + truth_table['DQA1_1'] + "+" + truth_table['DQA1_2'] +
                                '^' + truth_table['DQB1_1'] + "+" + truth_table['DQB1_2'])

truth_table['DR_GLString'] = (truth_table['DRB345_1'] + "+" + truth_table['DRB345_2'] +
                              '^' + truth_table['DRB1_1'] + "+" + truth_table['DRB1_2'])

truth_table['DQ_GLString'] = (truth_table['DQA1_1'] + "+" + truth_table['DQA1_2'] +
                              '^' + truth_table['DQB1_1'] + "+" + truth_table['DQB1_2'])

# Create truth dictionary with {genotype: [IDs]}
TT_DRDQ = defaultdict(list)
for line in range(len(truth_table)):
    subject_id = truth_table.loc[line, 'ID']
    GENO_DRDQ = truth_table.loc[line, which_impute + '_GLString']

    if subject_id not in TT_DRDQ:
        TT_DRDQ[GENO_DRDQ].append(subject_id)
    else:
        continue

print("Number of Unique Genotypes in Truth Table: ", len(TT_DRDQ))
print("Number of Pairings for Truth Table Genotypes: ", math.comb(len(TT_DRDQ), 2))

# Create truth dictionary with {ID: genotype: 1}
TT_DRDQ = {}
for line in range(len(truth_table)):
    subject_id = truth_table.loc[line, 'ID']
    GENO_DRDQ = truth_table.loc[line, which_impute + '_GLString']

    if GENO_DRDQ not in TT_DRDQ:
        TT_DRDQ[subject_id] = GENO_DRDQ
    else:
        continue

print("Number of IDs in Truth Table: ", len(TT_DRDQ))
truth_id = truth_table['ID'].to_list()
print("Number of Pairings for Truth Table IDs: ", math.comb(len(TT_DRDQ), 2))

# Simulate highres "donor-recipient" pairs from the truth dataset by looping through all combinations
eplet_truth = pd.DataFrame()
for don, rec in combinations(truth_id, 2):
    pair_id = "+".join([don, rec])
    don_geno = TT_DRDQ.get(don)
    rec_geno = TT_DRDQ.get(rec)

    # Clean up the strings to API format
    immunizer_alleles = don_geno
    immunizer_alleles = immunizer_alleles.replace('^', ',').replace('+', ',')
    patient_alleles = rec_geno
    patient_alleles = patient_alleles.replace('^', ',').replace('+', ',')

    donrec_line = pd.DataFrame({'DON_ID': don, 'DON_GLString': immunizer_alleles, 'REC_ID': rec, 'REC_GLString': patient_alleles}, index=[pair_id])
    eplet_truth = pd.concat([eplet_truth, donrec_line])


# Clean the format a little bit to only the columns you want
eplet_truth_table = eplet_truth.sort_index().reset_index(drop=True)
print(eplet_truth_table)
eplet_truth_table.to_csv(which_impute + '_pairs_truth.csv', header=True, index=False)  # work with eplets as a CSV and not worry about the JSON formatting

# Create file for multiple imputations
# Create a dictionary of the dataset such that GF_DRDQ[GENO_DRDQ][subject_id]
GF_DRDQ = defaultdict(list)
for line in range(len(impute_lowres)):
    subject_id = impute_lowres.loc[line, 'ID']
    GENO_DRDQ = impute_lowres.loc[line, which_impute + '_GLString']
    DRDQ_freq = impute_lowres.loc[line, which_impute + '_freq']

    if subject_id not in GF_DRDQ[GENO_DRDQ]:
        GF_DRDQ[GENO_DRDQ].append(subject_id)

print("Number of Unique Genotypes in Multiple Imputation: ", len(GF_DRDQ))
print("Number of Pairings for Multiple Imputation Genotypes: ", math.comb(len(GF_DRDQ), 2))

# Create a dictionary of the dataset such that GF_DRDQ[subject_ID][GENO_DRDQ][DRDQ_freq]
GF_DRDQ = defaultdict(dict)
for line in range(len(impute_lowres)):
    subject_id = impute_lowres.loc[line, 'ID']
    GENO_DRDQ = impute_lowres.loc[line, which_impute + '_GLString']
    DRDQ_freq = impute_lowres.loc[line, which_impute + '_freq']

    if GENO_DRDQ not in GF_DRDQ[subject_id]:
        GF_DRDQ[subject_id][GENO_DRDQ] = DRDQ_freq

print("Number of lines in Multiple Imputation: ", len(impute_lowres))
print("Number of Pairings for Multiple Imputations lines: ", math.comb(len(impute_lowres), 2))

# Use imputation output, compute eplet MM for each possible donor and recipient and store it into a DataFrame
eplet_DRDQ = pd.DataFrame()
for don, rec in combinations(truth_id, 2):
    pair_id = "+".join([don, rec])

    for donor_geno_DRDQ in GF_DRDQ[don]:
        for recip_geno_DRDQ in GF_DRDQ[rec]:
            donrec = geno_pairs(GF_DRDQ, donor_geno_DRDQ, recip_geno_DRDQ, don, rec, pair_id, which_impute)
            eplet_DRDQ = pd.concat([eplet_DRDQ, donrec])

eplet_impute = eplet_DRDQ.sort_index().reset_index(drop=True)  # Only need these columns for now
print(eplet_impute)
eplet_impute.to_csv(which_impute + '_pairs_imputation.csv', header=True, index=False)

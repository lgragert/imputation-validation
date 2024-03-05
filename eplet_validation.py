
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

# File that has looks at the highest and compares it to the truth

# Rename all the columns to what they need to be and create one dataFrame
def merge_dfs(DRDQ, DR, DQ, true_impute):
    DRDQ = DRDQ.rename(columns={'ALL_quantity': true_impute + '_ALL_quantity'})
    DR = DR.rename(columns={'DRB_quantity': true_impute + '_DR_quantity'})
    DQ = DQ.rename(columns={'DQ_quantity': true_impute + '_DQ_quantity'})

    merge_on = ['ID']

    DRDQ_DR = pd.merge(DRDQ, DR, on=merge_on, how='outer')
    DRDQ_DQ = pd.merge(DRDQ_DR, DQ, on=merge_on, how='outer')

    if true_impute == 'True':
        DRDQ_DQ[['DON_ID', 'REC_ID']] = DRDQ_DQ['ID'].str.split('+', expand=True)
        DRDQ_DQ = DRDQ_DQ.drop(columns=['ALL_details', 'DRB_details', 'DQ_details'])

    return DRDQ_DQ


# Create dictionaries for each pairing
def impute_dict(pairing_rows, DRDQ_dict, pair_ids, which_impute, which_eplet, geno_count):
    for row in range(len(pairing_rows)):
        DRDQ_quant = pairing_rows.loc[row, which_eplet + '_quantity']
        DRDQ_pairprob = pairing_rows.loc[row, 'PairProb_' + which_impute]
        if which_impute == 'DRDQ' and geno_count == 'GENO':
            DRDQ_DON = pairing_rows.loc[row, 'DON_DRDQ']
            DRDQ_REC = pairing_rows.loc[row, 'REC_DRDQ']
            DRDQ_pair = DRDQ_DON + '+' + DRDQ_REC

            if DRDQ_pair not in DRDQ_dict[pair_ids]:
                DRDQ_dict[pair_ids][DRDQ_pair] = DRDQ_pairprob
            else:
                DRDQ_dict[pair_ids][DRDQ_pair] = DRDQ_dict[pair_ids][DRDQ_pair] + DRDQ_pairprob
        else:
            if DRDQ_quant not in DRDQ_dict[pair_ids]:
                DRDQ_dict[pair_ids][DRDQ_quant] = DRDQ_pairprob
            else:
                DRDQ_dict[pair_ids][DRDQ_quant] = DRDQ_dict[pair_ids][DRDQ_quant] + DRDQ_pairprob

    return DRDQ_dict


# Top probabilities from the impute eplet files
def top_impute_df(top_impute, eplet_dict, which_impute, which_eplet, geno_count):
    top_eplets = pd.DataFrame()
    for ids, values in eplet_dict.items():
        top_eplet = max(values, key=values.get)  # only want the most probable eplets/counts
        top_freq = values[top_eplet]             # only want the highest frequency corresponding with the eplets/counts

        if geno_count == 'GENO':
            line = pd.DataFrame({'Highest_HighRes_DRDQ_GenotypePair': top_eplet, 'Highest_EpletMM_Pair_Prob': top_freq},index=[ids])
        else:
            line = pd.DataFrame({which_eplet + '_quantity': top_eplet, 'True_' + which_impute + '_Impute_Pair_Prob': top_freq}, index=[ids])

        top_eplets = pd.concat([top_eplets, line])

    top_impute = pd.concat([top_impute, top_eplets], axis=1)
    top_impute = top_impute.reset_index(names=['ID'])

    return top_impute


n_pairs = 100
# All the eplet truth table CSVs
eplet_truth_filename = 'DRDQ_eplet_truth_table' + str(n_pairs) + '.csv'
DRDQ_truth = pd.read_csv(eplet_truth_filename, header=0)
DRDQ_truth = DRDQ_truth.sort_values(by=['ID']).reset_index(drop=True)
eplet_truth_filename = 'DR_eplet_truth_table' + str(n_pairs) + '.csv'
DR_truth = pd.read_csv(eplet_truth_filename, header=0)
DR_truth = DR_truth.sort_values(by=['ID']).reset_index(drop=True)
eplet_truth_filename = 'DQ_eplet_truth_table' + str(n_pairs) + '.csv'
DQ_truth = pd.read_csv(eplet_truth_filename, header=0)
DQ_truth = DQ_truth.sort_values(by=['ID']).reset_index(drop=True)

# All the eplet impute table CSVs
eplet_impute_filename = 'DRDQ_eplet_lowres_impute' + str(n_pairs) + '.csv'
DRDQ_impute = pd.read_csv(eplet_impute_filename, header=0)
eplet_impute_filename = 'DR_eplet_lowres_impute' + str(n_pairs) + '.csv'
DR_impute = pd.read_csv(eplet_impute_filename, header=0)
eplet_impute_filename = 'DQ_eplet_lowres_impute' + str(n_pairs) + '.csv'
DQ_impute = pd.read_csv(eplet_impute_filename, header=0)

# Rename columns for truth table so you can merge it later
DRDQ_truth = merge_dfs(DRDQ_truth, DR_truth, DQ_truth, 'True')

rand_pairs = {}
for row in range(len(DRDQ_truth)):
    pairing = DRDQ_truth.loc[row, 'ID']

    if pairing not in rand_pairs:
        rand_pairs[pairing] = 1
    else:
        rand_pairs[pairing] = rand_pairs[pairing] + 1

# Get top pairings from the impute files
DRDQ_count = defaultdict(dict)
DRDQ_GENO = defaultdict(dict)
DR_count = defaultdict(dict)
DQ_count = defaultdict(dict)
for pair in rand_pairs:
    don, rec = pair.split("+")

    pairing_DRDQ = DRDQ_impute[(DRDQ_impute['DON_ID'] == don) & (DRDQ_impute['REC_ID'] == rec)]
    pairing_DR = DR_impute[(DR_impute['DON_ID'] == don) & (DR_impute['REC_ID'] == rec)]
    pairing_DQ = DQ_impute[(DQ_impute['DON_ID'] == don) & (DQ_impute['REC_ID'] == rec)]

    if pairing_DRDQ.empty is True:
        pairing_DRDQ = DRDQ_impute[(DRDQ_impute['DON_ID'] == rec) & (DRDQ_impute['REC_ID'] == don)]
    if pairing_DR.empty is True:
        pairing_DR = DR_impute[(DR_impute['DON_ID'] == rec) & (DR_impute['REC_ID'] == don)]
    if pairing_DQ.empty is True:
        pairing_DQ = DQ_impute[(DQ_impute['DON_ID'] == rec) & (DQ_impute['REC_ID'] == don)]

    pairing_DRDQ = pairing_DRDQ.reset_index(drop=True)
    pairing_DR = pairing_DR.reset_index(drop=True)
    pairing_DQ = pairing_DQ.reset_index(drop=True)

    DRDQ_count = impute_dict(pairing_DRDQ, DRDQ_count, pair, 'DRDQ', 'ALL', 'count')
    DRDQ_GENO = impute_dict(pairing_DRDQ, DRDQ_GENO, pair, 'DRDQ', 'ALL', 'GENO')
    DR_count = impute_dict(pairing_DR, DR_count, pair, 'DR', 'DRB', 'count')
    DQ_count = impute_dict(pairing_DQ, DQ_count, pair, 'DQ', 'DQ', 'count')


# Get top counts and top eplet MMs for each level: DRDQ, DR, DQ
top_DRDQ_eplets = pd.DataFrame()
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DRDQ_count, 'DRDQ', 'ALL', 'count')
top_DRDQ_geno = pd.DataFrame()
top_DRDQ_geno = top_impute_df(top_DRDQ_geno, DRDQ_GENO, 'DRDQ', 'ALL', 'GENO')
top_DR_eplets = pd.DataFrame()
top_DR_eplets = top_impute_df(top_DR_eplets, DR_count, 'DR', 'DRB', 'count')
top_DQ_eplets = pd.DataFrame()
top_DQ_eplets = top_impute_df(top_DQ_eplets, DQ_count, 'DQ', 'DQ', 'count')

top_impute = merge_dfs(top_DRDQ_eplets, top_DR_eplets, top_DQ_eplets, 'Highest')
top_impute = pd.merge(top_impute, top_DRDQ_geno, on=['ID'], how='outer')

# Merge truth with top impute dfs
truth_impute = pd.merge(DRDQ_truth, top_impute, on=['ID'], how='outer')

# Absolute Difference Columns
truth_impute['True_Highest_ALL_abs_diff'] = abs(truth_impute['Highest_ALL_quantity'] - truth_impute['True_ALL_quantity'])
truth_impute['True_Highest_DR_abs_diff'] = abs(truth_impute['Highest_DR_quantity'] - truth_impute['True_DR_quantity'])
truth_impute['True_Highest_DQ_abs_diff'] = abs(truth_impute['Highest_DQ_quantity'] - truth_impute['True_DQ_quantity'])


truth_impute = truth_impute[['DON_ID', 'REC_ID', 'Highest_EpletMM_Pair_Prob', 'Highest_ALL_quantity', 'True_ALL_quantity', 'True_DRDQ_Impute_Pair_Prob',
                             'Highest_DR_quantity', 'True_DR_quantity', 'True_DR_Impute_Pair_Prob', 'Highest_DQ_quantity',
                             'True_DQ_quantity', 'True_DQ_Impute_Pair_Prob', 'True_Highest_ALL_abs_diff',
                             'True_Highest_DR_abs_diff', 'True_Highest_DQ_abs_diff', 'Highest_HighRes_DRDQ_GenotypePair']]

# Get truth genotype pairs from DRDQ_pairs_truth_100.csv
truth_filename = 'DRDQ_pairs_truth_100.csv'
truth_geno = pd.read_csv(truth_filename, header=0)
truth_geno['True_HighRes_DRDQ_GenotypePair'] = truth_geno['DON_GLString'] + '+' + truth_geno['REC_GLString']
truth_geno = truth_geno[['DON_ID', 'REC_ID', 'True_HighRes_DRDQ_GenotypePair']]

truth_impute = pd.merge(truth_impute, truth_geno, on=['DON_ID', 'REC_ID'], how='outer')
truth_impute = truth_impute.sort_values(by=['Highest_EpletMM_Pair_Prob'], ascending=False)
truth_impute.to_csv('truth_impute_eplet_validation.csv', header=True, index=False)

# Create a histogram of the absolute difference columns
abs_DRDQ = truth_impute['True_Highest_ALL_abs_diff']
abs_DRDQ = abs_DRDQ.to_numpy()
ax1 = plt.hist(abs_DRDQ, bins=30, edgecolor='black', color='blue')
ax1 = plt.title('Histogram Plot of Absolute Difference in Eplet Counts for DRDQ')
ax1 = plt.xlabel('Absolute Difference in Counts')
ax1 = plt.ylabel('Frequency')
ax1 = plt.savefig('Histo_DRDQ_absdiff.png')

abs_DR = truth_impute['True_Highest_DR_abs_diff']
abs_DR = abs_DR.to_numpy()
ax2 = plt.hist(abs_DR, bins=30, edgecolor='black', color='blue')
ax2 = plt.title('Histogram Plot of Absolute Difference in Eplet Counts for DR')
ax2 = plt.xlabel('Absolute Difference in Counts')
ax2 = plt.ylabel('Frequency')
ax2 = plt.savefig('Histo_DR_absdiff.png')

abs_DQ = truth_impute['True_Highest_DQ_abs_diff']
abs_DQ = abs_DQ.to_numpy()
ax3 = plt.hist(abs_DQ, bins=30, edgecolor='black', color='blue')
ax3 = plt.title('Histogram Plot of Absolute Difference in Eplet Counts for DQ')
ax3 = plt.xlabel('Absolute Difference in Counts')
ax3 = plt.ylabel('Frequency')
ax3 = plt.savefig('Histo_DQ_absdiff.png')

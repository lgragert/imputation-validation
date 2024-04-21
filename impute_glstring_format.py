
import pandas as pd
from collections import defaultdict
import gzip
import sys
import warnings
warnings.filterwarnings("ignore")
import pyard
ard = pyard.init("3520")

# Pick top probability from the low resolution imputation files
# Take the imputation files and get the highest ranked haplotype pair probability for each patient
filename_list = []
for input_file in sys.argv[1:]:
    filename_list.append(input_file)

# Dict to store the multilocus with their frequency per patient
multiloc_freq = defaultdict(dict)
ClassI_freq = defaultdict(dict)
DRDQ_freq = defaultdict(dict)
DR_freq = defaultdict(dict)
DQ_freq = defaultdict(dict)
sevenloc_freq = defaultdict(dict)
# Dicts to store allele level probs for genotypes
GF_A = defaultdict(dict)
GF_B = defaultdict(dict)
GF_C = defaultdict(dict)
GF_DRB345 = defaultdict(dict)
GF_DRB1 = defaultdict(dict)
GF_DQA1 = defaultdict(dict)
GF_DQB1 = defaultdict(dict)
GF_DPA1 = defaultdict(dict)
GF_DPB1 = defaultdict(dict)
for filename in filename_list:
    file_name = filename
    impute_outfile = gzip.open(file_name, "rt")

    # Overall frequency already calculated for haplo pairs in 9loc impute data, but we want to preprocess it
    happair_id_total = {}  # cumulative genotype frequency total
    for line in impute_outfile:
        (subject_id, rank, hap1, hap2, freq) = line.split(',')

        happair_freq = 0
        happair_freq = float(freq)

        if subject_id not in happair_id_total:
            happair_id_total[subject_id] = 0
        happair_id_total[subject_id] += float(happair_freq)

    impute_outfile.close()

    ndonor = len(happair_id_total)

    # reopen impute file for 2nd pass
    impute_outfile = gzip.open(file_name, "rt")

    # compute probabilties for each haplotype pair
    happair_probs = {}  # HLA probability distribution
    for line in impute_outfile:
        (subject_id, rank, hap1, hap2, freq) = line.split(',')

        (A_1, C_1, B_1, DRB345_1, DRB1_1, DQA1_1, DQB1_1, DPA1_1, DPB1_1) = hap1.split('~')
        (A_2, C_2, B_2, DRB345_2, DRB1_2, DQA1_2, DQB1_2, DPA1_2, DPB1_2) = hap2.split('~')

        happair_freq = float(freq)

        # sort strings for the glstring format
        (A_1, A_2) = sorted([A_1, A_2])
        (B_1, B_2) = sorted([B_1, B_2])
        (C_1, C_2) = sorted([C_1, C_2])
        (DRB345_1, DRB345_2) = sorted([DRB345_1, DRB345_2])
        (DRB1_1, DRB1_2) = sorted([DRB1_1, DRB1_2])
        (DQA1_1, DQA1_2) = sorted([DQA1_1, DQA1_2])
        (DQB1_1, DQB1_2) = sorted([DQB1_1, DQB1_2])
        (DPA1_1, DPA1_2) = sorted([DPA1_1, DPA1_2])
        (DPB1_1, DPB1_2) = sorted([DPB1_1, DPB1_2])

        GENO_A = '+'.join([A_1, A_2])
        GENO_C = '+'.join([C_1, C_2])
        GENO_B = '+'.join([B_1, B_2])
        GENO_DRB345 = '+'.join([DRB345_1, DRB345_2])
        GENO_DRB1 = '+'.join([DRB1_1, DRB1_2])
        GENO_DQA1 = '+'.join([DQA1_1, DQA1_2])
        GENO_DQB1 = '+'.join([DQB1_1, DQB1_2])
        GENO_DPA1 = '+'.join([DPA1_1, DPA1_2])
        GENO_DPB1 = '+'.join([DPB1_1, DPB1_2])

        happair = '^'.join([GENO_A, GENO_C, GENO_B, GENO_DRB345, GENO_DRB1, GENO_DQA1, GENO_DQB1, GENO_DPA1, GENO_DPB1])

        prob = happair_freq / happair_id_total[subject_id]  # For genotype probability

        # For multilocus analysis for 9-loci, use the haplotype frequency
        if happair not in multiloc_freq[subject_id]:
            multiloc_freq[subject_id][happair] = happair_freq
        else:
            multiloc_freq[subject_id][happair] = multiloc_freq[subject_id][happair] + happair_freq

        # For singlelocus analysis, use genotype frequency
        if GENO_A not in GF_A[subject_id]:
            GF_A[subject_id][GENO_A] = prob
        else:
            GF_A[subject_id][GENO_A] = GF_A[subject_id][GENO_A] + prob

        if GENO_B not in GF_B[subject_id]:
            GF_B[subject_id][GENO_B] = prob
        else:
            GF_B[subject_id][GENO_B] = GF_B[subject_id][GENO_B] + prob

        if GENO_C not in GF_C[subject_id]:
            GF_C[subject_id][GENO_C] = prob
        else:
            GF_C[subject_id][GENO_C] = GF_C[subject_id][GENO_C] + prob

        if GENO_DRB345 not in GF_DRB345[subject_id]:
            GF_DRB345[subject_id][GENO_DRB345] = prob
        else:
            GF_DRB345[subject_id][GENO_DRB345] = GF_DRB345[subject_id][GENO_DRB345] + prob

        if GENO_DRB1 not in GF_DRB1[subject_id]:
            GF_DRB1[subject_id][GENO_DRB1] = prob
        else:
            GF_DRB1[subject_id][GENO_DRB1] = GF_DRB1[subject_id][GENO_DRB1] + prob

        if GENO_DQA1 not in GF_DQA1[subject_id]:
            GF_DQA1[subject_id][GENO_DQA1] = prob
        else:
            GF_DQA1[subject_id][GENO_DQA1] = GF_DQA1[subject_id][GENO_DQA1] + prob

        if GENO_DQB1 not in GF_DQB1[subject_id]:
            GF_DQB1[subject_id][GENO_DQB1] = prob
        else:
            GF_DQB1[subject_id][GENO_DQB1] = GF_DQB1[subject_id][GENO_DQB1] + prob

        if GENO_DPA1 not in GF_DPA1[subject_id]:
            GF_DPA1[subject_id][GENO_DPA1] = prob
        else:
            GF_DPA1[subject_id][GENO_DPA1] = GF_DPA1[subject_id][GENO_DPA1] + prob

        if GENO_DPB1 not in GF_DPB1[subject_id]:
            GF_DPB1[subject_id][GENO_DPB1] = prob
        else:
            GF_DPB1[subject_id][GENO_DPB1] = GF_DPB1[subject_id][GENO_DPB1] + prob

        # Create 7-loci probability
        sevenloc = '^'.join([GENO_A, GENO_C, GENO_B, GENO_DRB345, GENO_DRB1, GENO_DQA1, GENO_DQB1])
        if sevenloc not in sevenloc_freq[subject_id]:
            sevenloc_freq[subject_id][sevenloc] = prob
        else:
            sevenloc_freq[subject_id][sevenloc] = sevenloc_freq[subject_id][sevenloc] + prob

        # Create Class I probability
        classI = '^'.join([GENO_A, GENO_C, GENO_B])
        if classI not in ClassI_freq[subject_id]:
            ClassI_freq[subject_id][classI] = prob
        else:
            ClassI_freq[subject_id][classI] = ClassI_freq[subject_id][classI] + prob

        # Create DQ-DR probability
        dq_dr = '^'.join([GENO_DRB345, GENO_DRB1, GENO_DQA1, GENO_DQB1])
        if dq_dr not in DRDQ_freq[subject_id]:
            DRDQ_freq[subject_id][dq_dr] = prob
        else:
            DRDQ_freq[subject_id][dq_dr] = DRDQ_freq[subject_id][dq_dr] + prob

        # Create DR probability
        dr_only = '^'.join([GENO_DRB345, GENO_DRB1])
        if dr_only not in DR_freq[subject_id]:
            DR_freq[subject_id][dr_only] = prob
        else:
            DR_freq[subject_id][dr_only] = DR_freq[subject_id][dr_only] + prob

        # Create DQ probability
        dq_only = '^'.join([GENO_DQA1, GENO_DQB1])
        if dq_only not in DQ_freq[subject_id]:
            DQ_freq[subject_id][dq_only] = prob
        else:
            DQ_freq[subject_id][dq_only] = DQ_freq[subject_id][dq_only] + prob


# Reformat the Genotype Frequency dictionaries into DataFrames
def top_impute_df(top_sl_impute, geno_dict, locus, which_impute):
    top_singleloc = pd.DataFrame()
    for id, values in geno_dict.items():
        top_genotype = max(values, key=values.get)  # only want the haplotype highest frequency for now
        top_freq = values[top_genotype]             # only want the highest frequency corresponding with the haplotype
        if which_impute == 'singleloc':
            line = pd.DataFrame({'GENO_' + locus: top_genotype, 'GENO_' + locus + '_Prob': top_freq}, index=[id])
        elif which_impute == '9loc':
            line = pd.DataFrame({'9loc_GLString': top_genotype, 'HapPair_Prob': top_freq}, index=[id])
        elif which_impute == '7loc':
            line = pd.DataFrame({'7loc_GLString': top_genotype, '7loc_Prob': top_freq}, index=[id])
        elif which_impute == 'class I':
            line = pd.DataFrame({'ClassI_GLString': top_genotype, 'ClassI_Prob': top_freq}, index=[id])
        elif which_impute == 'DR-DQ':
            line = pd.DataFrame({'DRDQ_GLString': top_genotype, 'DRDQ_Prob': top_freq}, index=[id])
        elif which_impute == 'DR':
            line = pd.DataFrame({'DR_GLString': top_genotype, 'DR_Prob': top_freq}, index=[id])
        elif which_impute == 'DQ':
            line = pd.DataFrame({'DQ_GLString': top_genotype, 'DQ_Prob': top_freq}, index=[id])

        top_singleloc = pd.concat([top_singleloc, line])

    top_sl_impute = pd.concat([top_sl_impute, top_singleloc], axis=1)

    return top_sl_impute


# Keeps all imputations rather than the top probability
def all_impute(impute_all, geno_dict, which_impute):
    line_df = pd.DataFrame()
    for key, values in geno_dict.items():
        ids = key
        for geno, freq in values.items():
            if which_impute == 'DRDQ':
                line = pd.DataFrame({'ID': ids, 'DRDQ_GLString': geno, 'DRDQ_freq': freq}, index=[ids])
            elif which_impute == 'DR':
                line = pd.DataFrame({'ID': ids, 'DR_GLString': geno, 'DR_freq': freq}, index=[ids])
            else:
                line = pd.DataFrame({'ID': ids, 'DQ_GLString': geno, 'DQ_freq': freq}, index=[ids])
            line_df = pd.concat([line_df, line])

    impute_all = pd.concat([impute_all, line_df], axis=1)

    return impute_all


# For Class II eplet MM, keep all imputations and not just the top-imputation for DR-DQ
# Keep genotype probabilities of DRB345, DRB1, DQA1, DQB1
DRDQ_impute = pd.DataFrame()
DRDQ_impute = all_impute(DRDQ_impute, DRDQ_freq, 'DRDQ')
DR_impute = pd.DataFrame()
DR_impute = all_impute(DR_impute, DR_freq, 'DR')
DQ_impute = pd.DataFrame()
DQ_impute = all_impute(DQ_impute, DQ_freq, 'DQ')

print('Head of HLA DR-DQ imputations: \n', DRDQ_impute.head())
DRDQ_impute.to_csv('lowres_DRDQ_impute.csv', header=True, index=False)
DR_impute.to_csv('lowres_DR_impute.csv', header=True, index=False)
DQ_impute.to_csv('lowres_DQ_impute.csv', header=True, index=False)


# Get top impute with the new 9-loci level multilocus frequency
top_9loc_impute = pd.DataFrame()
top_9loc_impute = top_impute_df(top_9loc_impute, multiloc_freq, '', '9loc')
top_7loc_impute = pd.DataFrame()
top_7loc_impute = top_impute_df(top_7loc_impute, sevenloc_freq, '', '7loc')
top_classi_impute = pd.DataFrame()
top_classi_impute = top_impute_df(top_classi_impute, ClassI_freq, '', 'class I')
top_drdq_impute = pd.DataFrame()
top_drdq_impute = top_impute_df(top_drdq_impute, DRDQ_freq, '', 'DR-DQ')
top_dr_impute = pd.DataFrame()
top_dr_impute = top_impute_df(top_dr_impute, DR_freq, '', 'DR')
top_dq_impute = pd.DataFrame()
top_dq_impute = top_impute_df(top_dq_impute, DQ_freq, '', 'DQ')

# Get top impute with the new genotype frequencies
# The most probable genotype at a locus might be different than what was in the multilocus genotype
top_singleloc_impute = pd.DataFrame()
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_A, 'A', 'singleloc')
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_B, 'B', 'singleloc')
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_C, 'C', 'singleloc')
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_DRB1, 'DRB1', 'singleloc')
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_DRB345, 'DRB345', 'singleloc')
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_DQA1, 'DQA1', 'singleloc')
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_DQB1, 'DQB1', 'singleloc')
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_DPA1, 'DPA1', 'singleloc')
top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_DPB1, 'DPB1', 'singleloc')

top_singleloc_impute['SLUG_GLString'] = (top_singleloc_impute['GENO_A'] + '^' + top_singleloc_impute['GENO_C'] + '^' +
                                         top_singleloc_impute['GENO_B'] + '^' + top_singleloc_impute['GENO_DRB345'] + '^' +
                                         top_singleloc_impute['GENO_DRB1'] + '^' + top_singleloc_impute['GENO_DQA1'] + '^' +
                                         top_singleloc_impute['GENO_DQB1'] + '^' + top_singleloc_impute['GENO_DPA1'] + '^' +
                                         top_singleloc_impute['GENO_DPB1'])
top_singleloc_impute = top_singleloc_impute[['SLUG_GLString', 'GENO_A_Prob', 'GENO_B_Prob', 'GENO_C_Prob', 'GENO_DRB345_Prob',
                                             'GENO_DRB1_Prob', 'GENO_DQA1_Prob', 'GENO_DQB1_Prob', 'GENO_DPA1_Prob', 'GENO_DPB1_Prob']]

top_impute = pd.concat([top_9loc_impute, top_singleloc_impute, top_7loc_impute, top_classi_impute, top_drdq_impute, top_dr_impute, top_dq_impute], axis=1)
top_impute = top_impute.reset_index(names=['ID'])
print('Head of top imputations file: \n', top_impute.head())

# See if GLStrings are similar for MUG and SLUG analysis
boolean = pd.DataFrame()
boolean['Boolean'] = top_impute['9loc_GLString'] == top_impute['SLUG_GLString']
print('Amount of GLStrings that are the same for 9loc-MUG and SLUG: ', str(boolean['Boolean'].sum()) + "/" + str(len(top_impute)))

top_impute.to_csv('lowres_topprob_impute.csv', header=True, index=False)

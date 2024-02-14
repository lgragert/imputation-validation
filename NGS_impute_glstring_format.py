
import pandas as pd
from collections import defaultdict
import gzip

# Clean true genoytope high resolution data to standard format
# Pick top probability from the low resolution imputation files

# Add the letter in front of typing for each loci
def add_locus(ngs, letter):
    typing1 = ngs['NGS_' + letter + 'i']
    typing2 = ngs['NGS_' + letter + 'ii']

    if letter == 'DRB345':
        typing1 = typing1.fillna('X*NNNN')
        typing2 = typing2.fillna('X*NNNN')
        ngs['DON_' + letter] = 'DRB' + typing1 + '+' + 'DRB' + typing2
    else:
        ngs['DON_' + letter] = letter + "*" + typing1 + "+" + letter + '*' + typing2

    return ngs


filename = 'n217_withNGS.csv'
ngs_file = pd.read_csv(filename, header=0)

# Grab NGS typing and donor_race columns only
NGS = ngs_file[['unos', 'NGS_Ai', 'NGS_Aii', 'NGS_Bi', 'NGS_Bii', 'NGS_Ci', 'NGS_Cii', 'NGS_DRB1i', 'NGS_DRB1ii',
                'NGS_DRB345i', 'NGS_DRB345ii', 'NGS_DQA1i', 'NGS_DQA1ii', 'NGS_DQB1i', 'NGS_DQB1ii',
                'NGS_DPA1i', 'NGS_DPA1ii', 'NGS_DPB1i', 'NGS_DPB1ii', 'donorrace']]


# Add the appropriate letters together
loci = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
for locus in loci:
    NGS = add_locus(NGS, locus)

# Put the loci together in the GLString Format
NGS['GLString'] = NGS['DON_A'] + '^' + NGS['DON_B'] + '^' + NGS['DON_C'] + '^' + NGS['DON_DRB1'] + '^' + NGS['DON_DRB345'] + '^' + NGS['DON_DQA1'] + '^' + NGS['DON_DQB1'] + '^' + NGS['DON_DPA1'] + '^' + NGS['DON_DPB1']

# # Split the race column by race and ethnicity to get it into correct format
# NGS[['DON_RACE', 'DON_ETH']] = NGS['donorrace'].str.split(', ', expand=True)
# NGS['DON_RACE'] = NGS['DON_RACE'].replace({'White': 'CAU', 'Hispanic/Latino': 'HIS', 'Black': 'AFA', 'Asian': 'API', 'Amer Ind/Alaska Native': 'NAM'})
high_res = NGS[['unos', 'GLString']]
high_res = high_res.rename(columns={'unos': 'ID'})

high_res.to_csv('genotype_truth_table.csv', header=True, index=False)

# Check to see how many patients are in the dataset
print('Duplicates in NGS: ', len(high_res[high_res.duplicated(subset='ID')]))
print('Number of patients in NGS: ', len(high_res))


# Take the imputation files and get the highest ranked haplotype pair probability for each patient
pops = ['AFA', 'API', 'CAU', 'HIS', 'NAM']
# Dict to store the multilocus with their frequency per patient
multiloc_freq = defaultdict(dict)
for pop in pops:
    file_name = 'impute.nyulowres.' + pop + '.csv.gz'
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

        # happair = hap1 + '+' + hap2
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

        happair = '^'.join([GENO_A, GENO_B, GENO_C, GENO_DRB1, GENO_DRB345, GENO_DQA1, GENO_DQB1, GENO_DPA1, GENO_DPB1])

        if happair not in multiloc_freq[subject_id]:
            multiloc_freq[subject_id][happair] = happair_freq
        else:
            multiloc_freq[subject_id][happair] = multiloc_freq[subject_id][happair] + happair_freq

# TODO - Add more columns with most probable genotype per locus and the probability
# The most probable genotype at a locus might be different than what was in the multilocus genotype
# Expect many of the false negatives to become true positives

# Get top impute with the new multilocus frequency
top_impute = pd.DataFrame()
for key, value in multiloc_freq.items():
    id = key
    haplotype = max(value, key=value.get)  # only want the haplotype highest frequency for now
    freq = value[haplotype]                # only want the highest frequency corresponding with the haplotype
    line = pd.DataFrame({'GLString': haplotype, 'HapPair_Prob': freq}, index=[id])

    top_impute = pd.concat([top_impute, line])

top_impute = top_impute.reset_index(names=['ID'])

top_impute.to_csv('lowres_topprob_impute.csv', header=True, index=False)

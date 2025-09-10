import pandas as pd
from collections import defaultdict
import gzip
import sys
import pyard
import os
ard = pyard.init("3520")

filename_list = []
for input_file in sys.argv[1:]:
    filename_list.append(input_file)

# Define all possible loci in order
ALL_LOCI = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

# Dicts to store genotype frequencies for each locus
GF_dicts = {locus: defaultdict(dict) for locus in ALL_LOCI}
multiloc_freq = defaultdict(dict)
other_freqs = defaultdict(lambda: defaultdict(dict))  # for 7loc, classI, DRDQ, DR, DQ

# Track all loci found in all input files, in the order they appear in the first haplotype string
all_loci_found = []

def open_impute_file(filepath):
    if filepath.endswith('.gz'):
        return gzip.open(filepath, "rt")
    else:
        return open(filepath, "r")

for filename in filename_list:
    file_name = filename
    impute_outfile = open_impute_file(file_name)

    # Detect loci from first data line (order preserved)
    for line in impute_outfile:
        if line.strip() == "" or line.startswith("ID"):
            continue
        (_, _, hap1, _, _) = line.split(',', 4)
        loci_in_file = [h.split('*')[0] for h in hap1.strip().split('~')]
        # Add new loci in order of appearance
        for locus in loci_in_file:
            if locus not in all_loci_found:
                all_loci_found.append(locus)
        break
    impute_outfile.close()

    # 1st pass: sum total freq per subject
    happair_id_total = {}
    impute_outfile = open_impute_file(file_name)
    for line in impute_outfile:
        if line.strip() == "" or line.startswith("ID"):
            continue
        (subject_id, _, _, _, freq) = line.split(',', 4)
        happair_id_total.setdefault(subject_id, 0)
        happair_id_total[subject_id] += float(freq)
    impute_outfile.close()

    # 2nd pass: process haplotype pairs
    impute_outfile = open_impute_file(file_name)
    for line in impute_outfile:
        if line.strip() == "" or line.startswith("ID"):
            continue
        (subject_id, rank, hap1, hap2, freq) = line.strip().split(',', 4)
        hap1_alleles = hap1.split('~')
        hap2_alleles = hap2.split('~')
        happair_freq = float(freq)
        prob = happair_freq / happair_id_total[subject_id]

        # Sort alleles for each locus
        geno_locus = {}
        for idx, locus in enumerate(loci_in_file):
            alleles = sorted([hap1_alleles[idx], hap2_alleles[idx]])
            geno_locus[locus] = '+'.join(alleles)

        # Build multilocus GLString for present loci
        multiloc_glstring = '^'.join([geno_locus[locus] for locus in loci_in_file])
        multiloc_freq[subject_id][multiloc_glstring] = multiloc_freq[subject_id].get(multiloc_glstring, 0) + happair_freq

        # Single locus genotype probabilities
        for locus in loci_in_file:
            GF_dicts[locus][subject_id][geno_locus[locus]] = GF_dicts[locus][subject_id].get(geno_locus[locus], 0) + prob

        # Build other combinations only if loci are present
        # 7-loci
        sevenloc = [l for l in ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1'] if l in loci_in_file]
        if sevenloc:
            sevenloc_gl = '^'.join([geno_locus[l] for l in sevenloc])
            other_freqs['sevenloc'][subject_id][sevenloc_gl] = other_freqs['sevenloc'][subject_id].get(sevenloc_gl, 0) + prob
        # Class I
        classI = [l for l in ['A', 'C', 'B'] if l in loci_in_file]
        if classI:
            classI_gl = '^'.join([geno_locus[l] for l in classI])
            other_freqs['classI'][subject_id][classI_gl] = other_freqs['classI'][subject_id].get(classI_gl, 0) + prob
        # DRDQ
        drdq = [l for l in ['DRB345', 'DRB1', 'DQA1', 'DQB1'] if l in loci_in_file]
        if drdq:
            drdq_gl = '^'.join([geno_locus[l] for l in drdq])
            other_freqs['DRDQ'][subject_id][drdq_gl] = other_freqs['DRDQ'][subject_id].get(drdq_gl, 0) + prob
        # DR
        dr = [l for l in ['DRB345', 'DRB1'] if l in loci_in_file]
        if dr:
            dr_gl = '^'.join([geno_locus[l] for l in dr])
            other_freqs['DR'][subject_id][dr_gl] = other_freqs['DR'][subject_id].get(dr_gl, 0) + prob
        # DQ
        dq = [l for l in ['DQA1', 'DQB1'] if l in loci_in_file]
        if dq:
            dq_gl = '^'.join([geno_locus[l] for l in dq])
            other_freqs['DQ'][subject_id][dq_gl] = other_freqs['DQ'][subject_id].get(dq_gl, 0) + prob
    impute_outfile.close()

# Helper to get top impute for any dictionary
def top_impute_df(top_df, geno_dict, locus, which_impute):
    top_singleloc = pd.DataFrame()
    for id, values in geno_dict.items():
        top_genotype = max(values, key=values.get)
        top_freq = values[top_genotype]
        if which_impute == 'singleloc':
            line = pd.DataFrame({'GENO_' + locus: top_genotype, 'GENO_' + locus + '_Prob': top_freq}, index=[id])
        elif which_impute == 'multiloc':
            line = pd.DataFrame({'GLString': top_genotype, 'HapPair_Prob': top_freq}, index=[id])
        else:
            line = pd.DataFrame({which_impute + '_GLString': top_genotype, which_impute + '_Prob': top_freq}, index=[id])
        top_singleloc = pd.concat([top_singleloc, line])
    top_df = pd.concat([top_df, top_singleloc], axis=1)
    return top_df

# Helper to output all imputations for eplet analysis
def all_impute(impute_all, geno_dict, which_impute):
    line_df = pd.DataFrame()
    for key, values in geno_dict.items():
        for geno, freq in values.items():
            line = pd.DataFrame({'ID': key, which_impute + '_GLString': geno, which_impute + '_freq': freq}, index=[key])
            line_df = pd.concat([line_df, line])
    impute_all = pd.concat([impute_all, line_df], axis=1)
    return impute_all

# Output all DRDQ, DR, DQ imputations for eplet analysis
for which in ['DRDQ', 'DR', 'DQ']:
    if other_freqs[which]:
        df = pd.DataFrame()
        df = all_impute(df, other_freqs[which], which)
        df.to_csv(f'lowres_{which}_impute.csv', header=True, index=False)

# Build top impute DataFrame for multilocus
top_multiloc_impute = pd.DataFrame()
top_multiloc_impute = top_impute_df(top_multiloc_impute, multiloc_freq, '', 'multiloc')

# Build top impute DataFrame for single loci present
top_singleloc_impute = pd.DataFrame()
for locus in all_loci_found:
    if any(GF_dicts[locus].values()):
        top_singleloc_impute = top_impute_df(top_singleloc_impute, GF_dicts[locus], locus, 'singleloc')

# Build GLString for single locus genotype (SLUG) for present loci only (preserve order)
if not top_singleloc_impute.empty:
    slug_cols = [f'GENO_{locus}' for locus in all_loci_found if f'GENO_{locus}' in top_singleloc_impute.columns]
    top_singleloc_impute['SLUG_GLString'] = top_singleloc_impute[slug_cols].agg('^'.join, axis=1)
    prob_cols = [f'GENO_{locus}_Prob' for locus in all_loci_found if f'GENO_{locus}_Prob' in top_singleloc_impute.columns]
    top_singleloc_impute = top_singleloc_impute[['SLUG_GLString'] + prob_cols]

# Build top impute DataFrames for other combinations present
other_top_imputes = []
for which in ['sevenloc', 'classI', 'DRDQ', 'DR', 'DQ']:
    if other_freqs[which]:
        df = pd.DataFrame()
        df = top_impute_df(df, other_freqs[which], '', which)
        other_top_imputes.append(df)

# Concatenate all top impute DataFrames
dfs_to_concat = [df for df in [top_multiloc_impute, top_singleloc_impute] + other_top_imputes if not df.empty]
if dfs_to_concat:
    top_impute = pd.concat(dfs_to_concat, axis=1)
    top_impute = top_impute.reset_index(names=['ID'])
    print('Head of top imputations file: \n', top_impute.head())
    top_impute.to_csv('lowres_topprob_impute.csv', header=True, index=False)
else:
    print("No imputations found for any loci in input files.")

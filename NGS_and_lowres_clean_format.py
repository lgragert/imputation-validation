
import pandas as pd
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


# Split the race column by race and ethnicity to get it into correct format
NGS[['DON_RACE', 'DON_ETH']] = NGS['donorrace'].str.split(', ', expand=True)
NGS['DON_RACE'] = NGS['DON_RACE'].replace({'White': 'CAU', 'Hispanic/Latino': 'HIS', 'Black': 'AFA', 'Asian': 'API', 'Amer Ind/Alaska Native': 'NAM'})

high_res = NGS[['unos', 'GLString', 'DON_RACE']]
high_res = high_res.rename(columns={'unos': 'ID'})

high_res.to_csv('true_genotype_glstring.csv', header=True, index=False)

# Check to see how many patients are in the dataset
print('Duplicates in NGS: ', len(high_res[high_res.duplicated(subset='ID')]))
print('Number of patients in NGS: ', len(high_res))


# Take the imputation files and get the highest ranked probability for each patient
top_impute = pd.DataFrame()
pops = ['AFA', 'API', 'CAU', 'HIS', 'NAM']
for pop in pops:
    file_name = 'impute.nyulowres.' + pop + '.csv.gz'
    impute_file = pd.read_csv(file_name, header=None, names=['ID', 'Rank', 'Hap1', 'Hap2', 'HapPair_Prob'])

    top_prob = impute_file[impute_file['Rank'] == 1]
    top_impute = pd.concat([top_impute, top_prob])


top_impute = top_impute.reset_index(drop=True)
print('Number of patients in LowRes: ', len(top_impute))
print('Duplicates in top imputation: ', len(top_impute[top_impute.duplicated(subset=['ID'])]))

# Check to see if there are any missing IDs from highres dataset
missing_ids = high_res[~high_res.ID.isin(top_impute.ID)]
print('Missing in imputation files: ', missing_ids['ID'].to_list())

top_impute.to_csv('lowres_top_prob_imputation.csv', header=True, index=False)

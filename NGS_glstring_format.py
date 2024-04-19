
import pandas as pd
from collections import defaultdict
import gzip
import sys
import warnings
warnings.filterwarnings("ignore")
import pyard
ard = pyard.init("3520")

# Clean true genoytope from NGS high resolution data into standard format
# Assumes Next Generation Sequencing. Cleans data up to go from '01:01' to 'A*01:01' data, and then to GLString format

# Add the letter in front of typing for each loci, roll back to lgx typing, and sort it
def add_locus(ngs, letter):

    for line in range(len(ngs)):
        typing1 = ngs.loc[line, 'NGS_' + letter + 'i']
        typing1 = str(typing1)
        typing2 = ngs.loc[line, 'NGS_' + letter + 'ii']
        typing2 = str(typing2)

        if letter == 'DRB345':
            if typing1 == 'nan':
                typing1 = 'DRBX*NNNN'
            elif typing1 == '5*01:08N':   # When using redux it creates string of 'DRB5*01:02/DRB5*01:08' -> DRB5*01:08 does not work in HLAGenie, so chose DRB5*01:02
                typing1 = 'DRB5*01:02'
            else:
                typing1 = 'DRB' + typing1
                typing1 = ard.redux(typing1, 'lgx')
            if typing2 == 'nan':
                typing2 = 'DRBX*NNNN'
            elif typing2 == '5*01:08N':
                typing2 = 'DRB5*01:02'
            else:
                typing2 = 'DRB' + typing2
                typing2 = ard.redux(typing2, 'lgx')
            (typing1, typing2) = sorted([typing1, typing2])
            ngs.loc[line, 'DON_' + letter] = typing1 + '+' + typing2
        else:
            typing1 = letter + '*' + typing1
            typing2 = letter + '*' + typing2
            (typing1, typing2) = sorted([typing1, typing2])

            if typing1.endswith('NEW'):
                typing1 = typing1
            else:
                typing1 = ard.redux(typing1, 'lgx')

            if typing2.endswith('NEW'):
                typing2 = typing2
            else:
                typing2 = ard.redux(typing2, 'lgx')

            ngs.loc[line, 'DON_' + letter] = typing1 + "+" + typing2

    return ngs


filename = sys.argv[1]  # High resolution NGS data file
NGS = pd.read_csv(filename, header=0)

# Add the appropriate letters to the right locus. This takes data that is '01:01' and turns it into 'A*01:01'
loci = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
for locus in loci:
    NGS = add_locus(NGS, locus)

# Put the loci together in the GLString Format
NGS['GLString'] = (NGS['DON_A'] + '^' + NGS['DON_C'] + '^' + NGS['DON_B'] + '^' + NGS['DON_DRB345'] + '^' + NGS['DON_DRB1'] +
                   '^' + NGS['DON_DQA1'] + '^' + NGS['DON_DQB1'] + '^' + NGS['DON_DPA1'] + '^' + NGS['DON_DPB1'])

# Only keep the GLString and the ID column
high_res = NGS[['unos', 'GLString']]
high_res = high_res.rename(columns={'unos': 'ID'})

# Check to see how many patients are in the dataset
print('Duplicates in NGS: ', len(high_res[high_res.duplicated(subset='ID')]))
high_res = high_res.drop_duplicates(subset=['ID']).reset_index(drop=True)
print('Number of patients in NGS: ', len(high_res))

high_res.to_csv('genotype_truth_table.csv', header=True, index=False)  # Will use this file as the truth data


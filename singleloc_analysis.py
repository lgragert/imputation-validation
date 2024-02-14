
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, RocCurveDisplay
import matplotlib.pyplot as plt
import pyard


ard = pyard.init("3520")

# Compare the imputation to the truth table for single locus unphased genotype analysis
# If both alleles are the same, then that is a positive prediction=1, and if one/both are incorrect then that is a negative prediction=-1


# Separate the GLString out into locus pairings
def sep_glstring(file, high_res):
    file[['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']] = file['GLString'].str.split('^', expand=True)
    file = file.drop(columns=['GLString'])

    loci = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    for locus in loci:
        file[[locus + '_1', locus + '_2']] = file[locus].str.split('+', expand=True)
        file = file.drop(columns=[locus])

    # If your genotype truth table is not actually high resolution, then drop DPA1 and DPB1
    if high_res is False:
        file = file.drop(columns=['DRB345_1', 'DRB345_2', 'DQA1_1', 'DQA1_2', 'DPA1_1', 'DPA1_2', 'DPB1_1', 'DPB1_2'])
    return file


# Count the number of incorrect predictions and then make it in terms that we use to make negative(-1) and positive(1) predictions
def neg_prediction(truth_typ1,truth_typ2,impute_typ1,impute_typ2):
    if (truth_typ1 == "MISSING"):
        return "NA"
    if (impute_typ1 == "MISSING"):
        return "NA"

    # print ("True Genotype: " + truth_typ1 + "+" + truth_typ2)
    # print ("Top Imputed Genotype: " + impute_typ1 + "+" + impute_typ2) 

    donor_homoz = 0
    if (truth_typ1 == truth_typ2):
        donor_homoz = 1

    neg_count = 0
    if ((truth_typ1 != impute_typ1) & (truth_typ1 != impute_typ2)):
        neg_count += 1

    if ((truth_typ2 != impute_typ1) & (truth_typ2 != impute_typ2)):
        neg_count += 1

    if ((neg_count == 2) & (donor_homoz == 1)):
        neg_count = 1

    # Want to make it to where 1=correct prediction and -1=incorrect prediction
    if neg_count >= 1:
        neg_count = -1
    else:
        neg_count = 1

    return neg_count

truth_filename = 'genotype_truth_table.csv'
impute_filename = 'lowres_topprob_impute.csv'
truth_table = pd.read_csv(truth_filename, header=0)
impute = pd.read_csv(impute_filename, header=0)

truth_table = truth_table.drop_duplicates(subset=['ID'])  # TODO - bring this over to clean data script
truth_table = truth_table[truth_table.ID.isin(impute.ID)].reset_index(drop=True)  # Makes sure they are the same length and looking at the same patients
truth_table = truth_table.sort_values(by=['ID']).reset_index(drop=True)  # Sorts the patients, so each patient is in the same row as the imputation rows
impute = impute.sort_values(by=['ID']).reset_index(drop=True)

# print (truth_table)

high_res = False
truth_table = sep_glstring(truth_table, high_res)
impute = sep_glstring(impute, high_res)

if high_res == False:
    loci = ['A', 'B', 'C', 'DRB1', 'DQB1']
else:
    loci = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

for line in range(len(truth_table)):
    # print ("ID: " + truth_table.loc[line, 'ID'])
    for locus in loci:
        truth1 = truth_table.loc[line, locus + '_1']
        truth2 = truth_table.loc[line, locus + '_2']
        (truth1, truth2) = sorted([truth1, truth2])

        # py-ARD lgx rollup on truth
        truth1 = ard.redux(truth1,'lgx')
        truth2 = ard.redux(truth2,'lgx')

        impute1 = impute.loc[line, locus + '_1']
        impute2 = impute.loc[line, locus + '_2']
        (impute1, impute2) = sorted([impute1, impute2])

        impute.loc[line, locus + '_True'] = neg_prediction(truth1, truth2, impute1, impute2)

        # Create column for the probability taken from the probability
        probability = impute.loc[line, 'HapPair_Prob']
        if probability >= 0.5:
            impute.loc[line, 'Prob'] = 1
        else:
            impute.loc[line, 'Prob'] = -1

probability = impute['HapPair_Prob'].to_numpy()
probs = impute['Prob'].to_numpy()
for locus in loci:
    print('Confusion Matrix for Locus: ', locus)
    true_pred = impute[locus + '_True'].to_numpy()
    print(confusion_matrix(true_pred, probs))
    tn, fp, fn, tp = confusion_matrix(true_pred, probs).ravel()
    print('TP: ' + str(tp) + ', TN: ' + str(tn) + ', FP: ' + str(fp) + ', FN: ' + str(fn))
    print(classification_report(true_pred, probs))
    print('Brier Score Loss: ', brier_score_loss(true_pred, probability > 0.5))
    print('ROC-AUC Score: ', roc_auc_score(true_pred, probs))


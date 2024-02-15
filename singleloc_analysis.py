
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, RocCurveDisplay
from sklearn.calibration import CalibratedClassifierCV, CalibrationDisplay, calibration_curve
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pyard
ard = pyard.init("3520")

# Compare the imputation to the truth table for single locus unphased genotype analysis
# If both alleles are the same, then that is a positive prediction=1, and if one/both are incorrect then that is a negative prediction=-1


# Separate the GLString out into locus pairings
def sep_glstring(file, high_res, glstring):
    if glstring == 'GLString':
        file[['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']] = file['GLString'].str.split('^', expand=True)
        file = file.drop(columns=['GLString'])
    else:
        file[['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']] = file['SLUG_GLString'].str.split('^', expand=True)
        file = file.drop(columns=['SLUG_GLString', 'MUG_GLString', 'HapPair_Prob'])

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
    if truth_typ1 == "MISSING":
        return "NA"
    if impute_typ1 == "MISSING":
        return "NA"

    # print ("True Genotype: " + truth_typ1 + "+" + truth_typ2)
    # print ("Top Imputed Genotype: " + impute_typ1 + "+" + impute_typ2)

    # Cannot assume they are in exact order, so create a count of negative predictions
    donor_homoz = 0
    if truth_typ1 == truth_typ2:
        donor_homoz = 1

    neg_count = 0
    if (truth_typ1 != impute_typ1) & (truth_typ1 != impute_typ2):
        neg_count += 1

    if (truth_typ2 != impute_typ1) & (truth_typ2 != impute_typ2):
        neg_count += 1

    if (neg_count == 2) & (donor_homoz == 1):
        neg_count = 1

    # Want to make it to where 1=correct prediction and -1=incorrect prediction
    if neg_count >= 1:
        neg_count = 0
    else:
        neg_count = 1

    return neg_count


truth_filename = 'genotype_truth_table.csv'
impute_filename = 'lowres_topprob_impute.csv'
truth_table = pd.read_csv(truth_filename, header=0)
impute = pd.read_csv(impute_filename, header=0)

truth_table = truth_table[truth_table.ID.isin(impute.ID)].reset_index(drop=True)  # Makes sure they are the same length and looking at the same patients
truth_table = truth_table.sort_values(by=['ID']).reset_index(drop=True)  # Sorts the patients, so each patient is in the same row as the imputation rows
impute = impute.sort_values(by=['ID']).reset_index(drop=True)

# print (truth_table)

high_res = False
truth_table = sep_glstring(truth_table, high_res, 'GLString')
impute = sep_glstring(impute, high_res, 'SLUG_GLString')

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
        probability = impute.loc[line, 'GENO_' + locus + '_Prob']
        threshold = 0.95
        if probability > threshold:
            impute.loc[line, locus + '_Prob'] = 1
        else:
            impute.loc[line, locus + '_Prob'] = 0


confusion_mat = pd.DataFrame()
brier_loss = {}
roc_auc = {}
for locus in loci:
    # Need everything to be in numpy array for sklearn metrics
    probs = impute[locus + '_Prob'].to_numpy()                  # Probability of correctness in terms of -1/1
    probability = impute['GENO_' + locus + '_Prob'].to_numpy()  # Probability of correctness in [0,1] terms
    true_pred = impute[locus + '_True'].to_numpy()              # Actual prediction in terms of -1/1

    # Confusion Matrix per Locus
    print('Confusion Matrix for Locus: ', locus)
    print(confusion_matrix(true_pred, probs))
    tn, fp, fn, tp = confusion_matrix(true_pred, probs).ravel()
    print('TP: ' + str(tp) + ', TN: ' + str(tn) + ', FP: ' + str(fp) + ', FN: ' + str(fn))

    confusion_mat_line = pd.DataFrame({'TN': tn, 'FP': fp, 'FN': fn, 'TP': tp}, index=[locus])
    confusion_mat = pd.concat([confusion_mat, confusion_mat_line])

    # Classification Report
    print(classification_report(true_pred, probs))

    # Brier Score Loss
    brier_loss[locus] = brier_score_loss(true_pred, probability > threshold)
    print('Brier Score Loss: ', brier_loss[locus])

    # ROC-AUC Score
    roc_auc[locus] = roc_auc_score(true_pred, probability)
    print('ROC-AUC Score: ', roc_auc[locus])

    # Sort the Dataset based on low to high predictions and create four points
    impute_sort = impute.sort_values(by=['GENO_' + locus + '_Prob'])
    sort_probability = impute_sort['GENO_' + locus + '_Prob'].to_numpy()
    sort_true_pred = impute_sort[locus + '_True'].sort_values().to_numpy()
    # Create four points by taking the average in each bin
    prob1,prob2,prob3,prob4 = np.array_split(sort_probability, 4)
    true1,true2,true3,true4 = np.array_split(sort_true_pred, 4)
    probability_avg = [prob1.mean(), prob2.mean(), prob3.mean(), prob4.mean()]
    min_prob = np.min(probability_avg)
    true_avg = [true1.mean(), true2.mean(), true3.mean(), true4.mean()]
    min_true = np.min(true_avg)

    plt.figure(figsize=(8, 8))
    plt.plot(probability_avg, true_avg, marker='o', linestyle='-', label='Calibration Curve')
    plt.plot([0, 1], linestyle='--', label='Ideal Calibration')
    plt.xlabel('Mean Predicted Probability for ' + locus)
    plt.ylabel('Fraction of Predictions Correct')
    plt.title('Probability Calibration Curve for ' + locus)
    plt.legend()
    plt.show()

    RocCurveDisplay.from_predictions(true_pred, probability, plot_chance_level=True)
    plt.title('ROC-AUC for ' + locus)
    plt.show()


bf = pd.DataFrame({'Brier_Loss_Score': brier_loss}, index=brier_loss.keys())
ra = pd.DataFrame({'ROC-AUC': roc_auc}, index=roc_auc.keys())
confusion_mat = pd.concat([confusion_mat, bf], axis=1)
confusion_mat = pd.concat([confusion_mat, ra], axis=1)
print(confusion_mat)


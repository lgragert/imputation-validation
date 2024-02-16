
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, RocCurveDisplay, roc_curve
import matplotlib.pyplot as plt
import sys
import pyard
ard = pyard.init("3520")

# Compare the imputation to the truth table for single locus unphased genotype analysis
# If both alleles are the same, then that is a positive prediction=1, and if one/both are incorrect then that is a negative prediction=0


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


# Count the number of incorrect predictions and then make it in terms that we use to make negative(0) and positive(1) predictions
def neg_prediction(truth_typ1,truth_typ2,impute_typ1,impute_typ2):
    if truth_typ1 == "MISSING":
        return "NA"
    if impute_typ1 == "MISSING":
        return "NA"
    # print ("True Genotype: " + truth_typ1 + "+" + truth_typ2)
    # print ("Top Imputed Genotype: " + impute_typ1 + "+" + impute_typ2)

    if (truth_typ1 == impute_typ1) & (truth_typ2 == impute_typ2):
        neg_count = 1
    else:
        neg_count = 0

    # # Cannot assume they are in exact order, so create a count of negative predictions
    # donor_homoz = 0
    # if truth_typ1 == truth_typ2:
    #     donor_homoz = 1
    # neg_count = 0
    # if (truth_typ1 != impute_typ1) & (truth_typ1 != impute_typ2):
    #     neg_count += 1
    # if (truth_typ2 != impute_typ1) & (truth_typ2 != impute_typ2):
    #     neg_count += 1
    # if (neg_count == 2) & (donor_homoz == 1):
    #     neg_count = 1
    # # Want to make it to where 1=correct prediction and 0=incorrect prediction
    # if neg_count >= 1:
    #     neg_count = 0
    # else:
    #     neg_count = 1

    return neg_count


truth_filename = 'genotype_truth_table.csv'  # sys.argv[1]
impute_filename = 'lowres_topprob_impute.csv'   # sys.argv[2]
truth_table = pd.read_csv(truth_filename, header=0)
impute = pd.read_csv(impute_filename, header=0)

truth_table = truth_table[truth_table.ID.isin(impute.ID)].reset_index(drop=True)  # Makes sure they are the same length and looking at the same patients
truth_table = truth_table.sort_values(by=['ID']).reset_index(drop=True)  # Sorts the patients, so each patient is in the same row as the imputation rows
impute = impute.sort_values(by=['ID']).reset_index(drop=True)

# print (truth_table)

high_res = True
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

        impute1 = impute.loc[line, locus + '_1']
        impute2 = impute.loc[line, locus + '_2']
        (impute1, impute2) = sorted([impute1, impute2])

        # py-ARD lgx rollup on truth
        if truth1 == 'DRBX*NNNN':
            truth1 = 'DRBX*NNNN'
        elif truth1 == 'DPB1*NEW':
            truth1 = 'DPB1*NEW'
        elif truth1 == 'DQA1*NEW':
            truth1 = 'DQA1*NEW'
        else:
            truth1 = ard.redux(truth1, 'lgx')

        if truth2 == 'DRBX*NNNN':
            truth2 = 'DRBX*NNNN'
        elif truth2 == 'DPB1*NEW':
            truth2 = 'DPB1*NEW'
        elif truth2 == 'DQA1*NEW':
            truth2 = 'DQA1*NEW'
        else:
            truth2 = ard.redux(truth2,'lgx')

        impute.loc[line, locus + '_True'] = neg_prediction(truth1, truth2, impute1, impute2)

        # Create column for the probability taken from the probability
        probability = impute.loc[line, 'GENO_' + locus + '_Prob']
        threshold = 0.95
        if probability >= threshold:
            impute.loc[line, locus + '_Prob'] = 1
        else:
            impute.loc[line, locus + '_Prob'] = 0


confusion_mat = pd.DataFrame()
brier_loss = {}
roc_auc = {}
for locus in loci:
    # Need everything to be in numpy array for sklearn metrics
    probs = impute[locus + '_Prob'].to_numpy()                  # Probability of correctness in terms of 0/1
    probability = impute['GENO_' + locus + '_Prob'].to_numpy()  # Probability of correctness in [0,1] terms
    true_pred = impute[locus + '_True'].to_numpy()              # Actual prediction in terms of 0/1

    # Confusion Matrix per Locus
    print('Classification Report for Locus ', locus)
    # print(confusion_matrix(true_pred, probs))
    tn, fp, fn, tp = confusion_matrix(true_pred, probs).ravel()
    # print('TP: ' + str(tp) + ', TN: ' + str(tn) + ', FP: ' + str(fp) + ', FN: ' + str(fn))

    confusion_mat_line = pd.DataFrame({'TN': tn, 'FP': fp, 'FN': fn, 'TP': tp}, index=[locus])
    confusion_mat = pd.concat([confusion_mat, confusion_mat_line])

    # Classification Report
    print(classification_report(true_pred, probs))

    # Brier Score Loss
    brier_loss[locus] = brier_score_loss(true_pred, probability > threshold)
    # print('Brier Score Loss: ', brier_loss[locus])

    # ROC-AUC Score
    roc_auc[locus] = roc_auc_score(true_pred, probability)
    # print('ROC-AUC Score: ', roc_auc[locus])

    # Sort the Dataset based on probabilities to allow the low/high probabilities to be together in the split
    impute_sort = impute.sort_values(by=['GENO_' + locus + '_Prob'])
    sort_probability = impute_sort['GENO_' + locus + '_Prob'].to_numpy()
    sort_true_pred = impute_sort[locus + '_True'].to_numpy()

    # Create four points by taking the average in each bin
    n_bins = 4
    prob_q1,prob_q2,prob_q3,prob_q4 = np.array_split(sort_probability, n_bins)
    true_q1,true_q2,true_q3,true_q4 = np.array_split(sort_true_pred, n_bins)
    probability_avg = [prob_q1.mean(), prob_q2.mean(), prob_q3.mean(), prob_q4.mean()]
    min_prob = np.min(probability_avg)
    max_prob = np.max(probability_avg)
    min_prob_in_bin = [np.min(prob_q1), np.min(prob_q2), np.min(prob_q3), np.min(prob_q4)]
    max_prob_in_bin = [np.max(prob_q1), np.max(prob_q2), np.max(prob_q3), np.max(prob_q4)]
    true_avg = [true_q1.mean(), true_q2.mean(), true_q3.mean(), true_q4.mean()]
    min_true = np.min(true_avg)

    # Create the standard error bars for each point
    snd_err1 = np.std(prob_q1, ddof=1) / np.sqrt(np.size(prob_q1))
    snd_err2 = np.std(prob_q2, ddof=1) / np.sqrt(np.size(prob_q2))
    snd_err3 = np.std(prob_q3, ddof=1) / np.sqrt(np.size(prob_q3))
    snd_err4 = np.std(prob_q4, ddof=1) / np.sqrt(np.size(prob_q4))
    snd_err = [snd_err1, snd_err2, snd_err3, snd_err4]

    # Compute the city-block distance between the quantiles and the diagonal (x=y)
    city_block_dst = (abs(probability_avg[0] - true_avg[0]) + abs(probability_avg[1] - true_avg[1]) + abs(probability_avg[2] - true_avg[2]) + abs(probability_avg[3] - true_avg[3])) / n_bins


    # Create a table exactly like the print statements above to add to the bottom of the calibration plot
    table_data = [["Quantile", "Prob Avg", "True Fraction", "Min Prob", "Max Prob", "Standard Error"],
                  ['Q1', str(round(probability_avg[0], 4)), str(round(true_avg[0], 4)), str(round(min_prob_in_bin[0], 4)), str(round(max_prob_in_bin[0], 4)), str(round(snd_err1, 4))],
                  ['Q2', str(round(probability_avg[1], 4)), str(round(true_avg[1], 4)), str(round(min_prob_in_bin[1], 4)), str(round(max_prob_in_bin[1], 4)), str(round(snd_err2, 4))],
                  ['Q3', str(round(probability_avg[2], 4)), str(round(true_avg[2], 4)), str(round(min_prob_in_bin[2], 4)), str(round(max_prob_in_bin[2], 4)), str(round(snd_err3, 4))],
                  ['Q4', str(round(probability_avg[3], 4)), str(round(true_avg[3], 4)), str(round(min_prob_in_bin[3], 4)), str(round(max_prob_in_bin[3], 4)), str(round(snd_err4, 4))]]
    print('Quantile Statistics for Locus: ', locus)
    print(table_data[0])
    print(table_data[1])
    print(table_data[2])
    print(table_data[3])
    print(table_data[4])

    # Create a bar plot where it shows the distribution of predictions, have to separate it from plot so it does not get added in
    counts, bins, _ = plt.hist(sort_probability, bins=20)
    num_IDs = len(impute)
    fract_counts = counts / num_IDs  # This allows us to get the fraction (* 100 = %) of cases for each count from the histogram

    calibrat_plot = plt.figure(figsize=(8, 8))
    calibrat_plot = plt.errorbar(probability_avg, true_avg, yerr=snd_err, marker='o', linestyle='', label='True Fraction vs Probability Average for Quartile', color='red', ecolor='black', capsize=7)
    calibrat_plot = plt.plot([0,1], linestyle='--', label='Ideal Calibration', color='blue')
    calibrat_plot = plt.xlabel('Mean Predicted Probability for Quantile')
    calibrat_plot = plt.ylabel('Fraction of Predictions Correct')
    # plt.yscale('log')
    # plt.xscale('log')
    # calibrat_plot = plt.suptitle('Calibration Plot and Prediction Probability Distribution for HLA-' + locus + " Locus\n" + 'Brier Score Loss: ' + str(round(brier_loss[locus], 4)) + '\n                        \n')
    calibrat_plot = plt.bar(bins[:-1], fract_counts, width=np.diff(bins), edgecolor='black', color='grey')
    calibrat_plot = plt.xlim(0,1.05)
    calibrat_plot = plt.ylim(0,1.05)
    table = plt.table(cellText=table_data)
    table = table.scale(1, 1.5)
    ax = plt.gca()
    ax.xaxis.set_label_position('top')  # have to add x-axis to the top because of the table at the bottom
    ax.xaxis.set_ticks_position('top')
    ax.set_title('Calibration Plot and Prediction Probability Distribution for HLA-' + locus + " Locus\n" + 'Brier Score Loss: ' + str(round(brier_loss[locus], 4)) + '\n' + 'City-Block Distance: ' + str(round(city_block_dst, 4)), pad=20)  # Space between x-axis and title
    calibrat_plot = plt.legend()
    calibrat_plot = plt.savefig("Calibration_" + locus + ".png", bbox_inches='tight')
    # plt.show()

plt.figure(0)
for locus in loci:
    # Display ROC plots
    probability = impute['GENO_' + locus + '_Prob'].to_numpy()  # Probability of correctness in [0,1] terms
    true_pred = impute[locus + '_True'].to_numpy()  # Actual prediction in terms of 0/1
    # Shows the old code for how to do one ROC per locus
    # roc_plot = RocCurveDisplay.from_predictions(true_pred, probability, plot_chance_level=True)
    # roc_plot = plt.title('ROC Curve and AUC for ' + locus)
    # roc_plot = plt.savefig("ROC_" + locus + ".png", bbox_inches='tight')
    # # plt.show()
    fpr, tpr, thresholds = roc_curve(true_pred, probability)
    # Calculate Area under the curve to display on the plot
    auc = roc_auc_score(true_pred, probability)
    # Now, plot the computed values
    roc_plot = plt.plot(fpr, tpr, label='%s (AUC = %0.2f)' % (locus, auc))
    # Custom settings for the plot
roc_plot = plt.plot([0, 1], 'r--', color='black')
roc_plot = plt.xlabel('Specificity (False Positive Rate)')
roc_plot = plt.ylabel('Sensitivity (True Positive Rate)')
roc_plot = plt.title('ROC-AUC Curves for each Locus Predictor')
roc_plot = plt.legend(loc="lower right")
roc_plot = plt.savefig("ROC_AUC_9loc.png", bbox_inches='tight')

bf = pd.DataFrame({'Brier_Loss_Score': brier_loss}, index=brier_loss.keys())
ra = pd.DataFrame({'ROC-AUC': roc_auc}, index=roc_auc.keys())
confusion_mat = pd.concat([confusion_mat, bf], axis=1)
confusion_mat = pd.concat([confusion_mat, ra], axis=1)
print(confusion_mat)


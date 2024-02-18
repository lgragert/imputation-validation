
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, RocCurveDisplay
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# Compare the top imputation to the truth table for multiple multilocus unphased genotype analysis
# At 9-loci, 7-loci, DR-DQ, and Class I levels
# If the GLString is the same, then it is a correct match


# Use this definition to make negative(0) and positive(1) predictions
def neg_prediction(truth_typ,impute_typ, whichtype):
    for row in range(len(impute_typ)):
        impute_string = impute_typ.loc[row, whichtype + '_GLString']
        truth_string = truth_typ.loc[row, 'GLString']

        if truth_string == "MISSING":
            return "NA"
        if impute_string == "MISSING":
            return "NA"

        if truth_string == impute_string:
            neg_pred = 1
        else:
            neg_pred = 0

        impute_typ.loc[row, whichtype + '_True'] = neg_pred

    return impute_typ


# Separate GLString for truth tables for calculating predictions for anything other than 9-loci
def sep_glstring(file, whichlevel):
    file[['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']] = file['GLString'].str.split('^', expand=True)
    file = file.drop(columns=['GLString'])

    if whichlevel == '7-loc':
        file['GLString'] = file['A'] + '^' + file['B'] + '^' + file['C'] + '^' + file['DRB1'] + '^' + file['DRB345'] + '^' + file['DQA1'] + '^' + file['DQB1']
    elif whichlevel == 'DR-DQ':
        file['GLString'] = file['DRB1'] + '^' + file['DRB345'] + '^' + file['DQA1'] + '^' + file['DQB1']
    elif whichlevel == 'Class I':
        file['GLString'] = file['A'] + '^' + file['B'] + '^' + file['C']

    loci = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    for locus in loci:
        file = file.drop(columns=[locus])

    return file


truth_filename = 'genotype_truth_table.csv'  # sys.argv[1]
impute_filename = 'lowres_topprob_impute.csv'   # sys.argv[2]
truth_table = pd.read_csv(truth_filename, header=0)
impute = pd.read_csv(impute_filename, header=0)

truth_table = truth_table[truth_table.ID.isin(impute.ID)].reset_index(drop=True)  # Makes sure they are the same length and looking at the same patients
truth_table = truth_table.sort_values(by=['ID']).reset_index(drop=True)  # Sorts the patients, so each patient is in the same row as the imputation rows
impute = impute.sort_values(by=['ID']).reset_index(drop=True)

# 9-loci predition use HapPair Probability
impute_9loc = impute[['ID', '9loc_GLString', 'HapPair_Prob']]
impute_9loc = neg_prediction(truth_table, impute_9loc, '9loc')
print('Positive and Negative Predictions for 9-loci level:\n', impute_9loc['9loc_True'].value_counts())

# 7-loci prediction use A+B+C+DRB345+DRB1+DQA1+DQB1 Probabilities Combined Together
truth_7loc = sep_glstring(truth_table, '7-loc')
impute_7loc = impute[['ID', '7loc_GLString', '7loc_Prob']]
impute_7loc = neg_prediction(truth_7loc, impute_7loc, '7loc')
print('Positive and Negative Predictions for 7-loci level:\n', impute_7loc['7loc_True'].value_counts())

# DR+DQ prediction use DRB345+DRB1+DQA1+DQB1 Probabilities Combined Together
truth_DRDQ = sep_glstring(truth_table, 'DR-DQ')
impute_DRDQ = impute[['ID', 'DRDQ_GLString', 'DRDQ_Prob']]
impute_DRDQ = neg_prediction(truth_DRDQ, impute_DRDQ, 'DRDQ')
print('Positive and Negative Predictions for DR-DQ loci level:\n', impute_DRDQ['DRDQ_True'].value_counts())

# Class I prediction use A+B+C Probabilities Combined Together
truth_classI = sep_glstring(truth_table, 'Class I')
impute_classI = impute[['ID', 'ClassI_GLString', 'ClassI_Prob']]
impute_classI = neg_prediction(truth_classI, impute_classI, 'ClassI')
print('Positive and Negative Predictions for Class I loci level:\n', impute_classI['ClassI_True'].value_counts())


# Create Calibration Plots for each case by creating four points
def calibration_plot(impute_typ, which_typ):
    if which_typ == '9loc':
        probability = impute_typ['HapPair_Prob'].to_numpy()  # Probability of correctness in [0,1] terms
        impute_sort = impute_typ.sort_values(by=['HapPair_Prob'])
        sort_probability = impute_sort['HapPair_Prob'].to_numpy()
    else:
        impute_sort = impute_typ.sort_values(by=[which_typ + '_Prob'])
        sort_probability = impute_sort[which_typ + '_Prob'].to_numpy()
        probability = impute_typ[which_typ + '_Prob'].to_numpy()

    true_pred = impute_typ[which_typ + '_True'].to_numpy()
    threshold = 0.5
    brier_loss = brier_score_loss(true_pred, probability > threshold)
    sort_true_pred = impute_sort[which_typ + '_True'].to_numpy()

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

    # Mean Squared Error (MSE) for the bin averages
    mse_bins = np.square((abs(probability_avg[0] - true_avg[0]) + abs(probability_avg[1] - true_avg[1]) + abs(probability_avg[2] - true_avg[2]) + abs(probability_avg[3] - true_avg[3]))) / n_bins

    # Create a table exactly like the print statements above to add to the bottom of the calibration plot
    table_data = [["Quantile", "Prob Avg", "True Fraction", "Min Prob", "Max Prob", "Standard Error"],
                    ['Q1', str(round(probability_avg[0], 4)), str(round(true_avg[0], 4)), str(round(min_prob_in_bin[0], 4)), str(round(max_prob_in_bin[0], 4)), str(round(snd_err1, 4))],
                    ['Q2', str(round(probability_avg[1], 4)), str(round(true_avg[1], 4)), str(round(min_prob_in_bin[1], 4)), str(round(max_prob_in_bin[1], 4)), str(round(snd_err2, 4))],
                    ['Q3', str(round(probability_avg[2], 4)), str(round(true_avg[2], 4)), str(round(min_prob_in_bin[2], 4)), str(round(max_prob_in_bin[2], 4)), str(round(snd_err3, 4))],
                    ['Q4', str(round(probability_avg[3], 4)), str(round(true_avg[3], 4)), str(round(min_prob_in_bin[3], 4)), str(round(max_prob_in_bin[3], 4)), str(round(snd_err4, 4))]]

    print('Quantile Statistics for ' + which_typ + ':')
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
    calibrat_plot = plt.errorbar(probability_avg, true_avg, yerr=snd_err, marker='o', linestyle='', label='True Fraction vs Probability Average for Quantile', color='red', ecolor='black', capsize=7)
    calibrat_plot = plt.plot([0, 1], linestyle='--', label='Ideal Calibration', color='blue')
    calibrat_plot = plt.xlabel('Mean Predicted Probability for Quantile')
    calibrat_plot = plt.ylabel('Fraction of Predictions Correct')
    calibrat_plot = plt.bar(bins[:-1], fract_counts, width=np.diff(bins), edgecolor='black', color='grey')
    calibrat_plot = plt.xlim(0, 1.05)
    calibrat_plot = plt.ylim(0, 1.05)
    table = plt.table(cellText=table_data)
    table = table.scale(1, 1.5)
    ax = plt.gca()
    ax.xaxis.set_label_position('top')  # have to add x-axis to the top because of the table at the bottom
    ax.xaxis.set_ticks_position('top')
    if which_typ == '9loc':
        title = 'Calibration Plot and Prediction Probability Distribution for 9-loci\n'
    elif which_typ == '7loc':
        title = 'Calibration Plot and Prediction Probability Distribution for 7-loci\n'
    elif which_typ == 'ClassI':
        title = 'Calibration Plot and Prediction Probability Distribution for Class I loci\n'
    elif which_typ == 'DRDQ':
        title = 'Calibration Plot and Prediction Probability Distribution for DR-DQ loci\n'
    ax.set_title(title +
                 'Brier: ' + str(round(brier_loss, 4)) + ', Bin Avg MSE: ' +
                 str(round(mse_bins, 4)) + ', Bin Avg City-Block Dist: ' + str(round(city_block_dst, 4)), pad=20)  # Space between x-axis and title
    calibrat_plot = plt.legend()

    return calibrat_plot

calibrat_plot9loc = calibration_plot(impute_9loc, '9loc')
calibrat_plot9loc = plt.savefig("Calibration_9loc.png", bbox_inches='tight')

calibrat_plot7loc = calibration_plot(impute_7loc, '7loc')
calibrat_plot7loc = plt.savefig("Calibration_7loc.png", bbox_inches='tight')

calibrat_plotclassI = calibration_plot(impute_classI, 'ClassI')
calibrat_plotclassI = plt.savefig("Calibration_ClassI.png", bbox_inches='tight')

calibrat_plotDRDQ = calibration_plot(impute_DRDQ, 'DRDQ')
calibrat_plotDRDQ = plt.savefig("Calibration_DRDQ.png", bbox_inches='tight')


import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, RocCurveDisplay
import sys
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
    file[['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']] = file['GLString'].str.split('^', expand=True)
    file = file.drop(columns=['GLString'])

    if whichlevel == '7-loc':
        file['GLString'] = file['A'] + '^' + file['C'] + '^' + file['B'] + '^' + file['DRB345'] + '^' + file['DRB1'] + '^' + file['DQA1'] + '^' + file['DQB1']
    elif whichlevel == 'DR-DQ':
        file['GLString'] = file['DRB345'] + '^' + file['DRB1'] + '^' + file['DQA1'] + '^' + file['DQB1']
    elif whichlevel == 'Class I':
        file['GLString'] = file['A'] + '^' + file['C'] + '^' + file['B']

    loci = ['A', 'C', 'B', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    for locus in loci:
        file = file.drop(columns=[locus])

    return file


truth_filename = sys.argv[1]
impute_filename = sys.argv[2]
num_bins = int(sys.argv[3])  # amount of bins you want to create
truth_table = pd.read_csv(truth_filename, header=0)
impute = pd.read_csv(impute_filename, header=0)

truth_table = truth_table[truth_table.ID.isin(impute.ID)].reset_index(drop=True)  # Makes sure they are the same length and looking at the same patients
truth_table = truth_table.sort_values(by=['ID']).reset_index(drop=True)  # Sorts the patients, so each patient is in the same row as the imputation rows
impute = impute.sort_values(by=['ID']).reset_index(drop=True)

same_ids = truth_table['ID']
impute = pd.merge(same_ids, impute, how='inner', on='ID')
print('Amount of IDs in the imputation file: ', len(impute))
print('Amount of IDs in the truth table file: ', len(truth_table))
print(impute.head())

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
def calibration_plot(impute_typ, which_typ, n_bins):
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

    # Create how many quantiles you want by taking the average in each bin
    split_prob = np.array_split(sort_probability, n_bins)
    split_true = np.array_split(sort_true_pred, n_bins)
    probability_avg = [quantiles.mean() for quantiles in split_prob]
    min_prob = np.min(probability_avg)
    max_prob = np.max(probability_avg)
    min_prob_in_bin = [np.min(quantiles) for quantiles in split_prob]
    max_prob_in_bin = [np.max(quantiles) for quantiles in split_prob]
    true_avg = [quantiles.mean() for quantiles in split_true]
    min_true = np.min(true_avg)

    # Create the standard error bars for each point
    snd_err = [(np.std(quantiles, ddof=1) / np.sqrt(np.size(quantiles))) for quantiles in split_prob]

    # Compute the city-block distance between the quantiles and the diagonal (x=y)
    city_block_dst = (sum((abs(probability_avg[quantile] - true_avg[quantile]) for quantile in range(0, n_bins)))) / n_bins

    # Mean Squared Error (MSE) for the bin averages
    mse_bins = np.square((sum((abs(probability_avg[quantile] - true_avg[quantile]) for quantile in range(0, n_bins))))) / n_bins

    # Create a table exactly like the print statements above to add to the bottom of the calibration plot
    table_data = [
        ['Q' + str(quantiles + 1), str(round(probability_avg[quantiles], 4)), str(round(true_avg[quantiles], 4)),
         str(round(min_prob_in_bin[quantiles], 4)), str(round(max_prob_in_bin[quantiles], 4)),
         str(round(snd_err[quantiles], 4))] for quantiles in range(0, n_bins)]
    heading = ["Quantile", "Prob Avg", "True Fraction", "Min Prob", "Max Prob", "Standard Error"]
    table_data.insert(0, heading)

    print('Quantile Statistics for ' + which_typ + ':')
    for quantiles in range(0,n_bins + 1):
        print(table_data[quantiles])

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

calibrat_plot9loc = calibration_plot(impute_9loc, '9loc', num_bins)
calibrat_plot9loc = plt.savefig("Calibration_9loc.png", bbox_inches='tight')

calibrat_plot7loc = calibration_plot(impute_7loc, '7loc', num_bins)
calibrat_plot7loc = plt.savefig("Calibration_7loc.png", bbox_inches='tight')

calibrat_plotclassI = calibration_plot(impute_classI, 'ClassI', num_bins)
calibrat_plotclassI = plt.savefig("Calibration_ClassI.png", bbox_inches='tight')

calibrat_plotDRDQ = calibration_plot(impute_DRDQ, 'DRDQ', num_bins)
calibrat_plotDRDQ = plt.savefig("Calibration_DRDQ.png", bbox_inches='tight')

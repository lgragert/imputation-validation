
import pandas as pd
from collections import defaultdict
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import sys

# Take the Monte Carlo output and compare the truth to the probabilities and create calibration curves

# Reformat the Eplet Frequency dictionaries into DataFrames
def top_impute_df(top_impute, eplet_dict, count_str, which_impute):
    top_eplets = pd.DataFrame()
    for ids, values in eplet_dict.items():
        top_eplet = max(values, key=values.get)  # only want the most probable eplets
        top_freq = values[top_eplet]             # only want the highest frequency corresponding with the eplets

        line = pd.DataFrame({which_impute + "_" + count_str: top_eplet, which_impute + "_" + count_str + 'Prob': top_freq}, index=[ids])

        top_eplets = pd.concat([top_eplets, line])

    top_impute = pd.concat([top_impute, top_eplets], axis=1)

    return top_impute


# Count the number of incorrect predictions and then make it in terms that we use to make negative(0) and positive(1) predictions
def neg_prediction(truth_typ, impute_typ):
    # print ("True Genotype: " + truth_typ1 + "+" + truth_typ2)
    # print ("Top Imputed Genotype: " + impute_typ1 + "+" + impute_typ2)
    truth_typ = str(truth_typ)
    impute_typ = str(impute_typ)  # for anything that is NaN, better this way

    if truth_typ == impute_typ:
        neg_count = 1
    else:
        neg_count = 0

    return neg_count


truth_filename = "eplet_truth_table100.csv"
eplet_truth = pd.read_csv(truth_filename, header=0)

# Get all pairs from the truth table
rand_pairs = {}
for row in range(len(eplet_truth)):
    pair_id = eplet_truth.loc[row, 'ID']

    if pair_id not in rand_pairs:
        rand_pairs[pair_id] = 1
    else:
        rand_pairs[pair_id] = rand_pairs[pair_id] + 1

impute_filename = "eplet_lowres_impute100.csv"
eplet_imptue = pd.read_csv(impute_filename, header=0)

# Start with counts and go from there {DON+REC: {Count: prob}}
# Unique eplets do similar thing such as {DON+REC: {Eplet_str: prob}}
DRDQ_count = defaultdict(dict)
DR_count = defaultdict(dict)
DQ_count = defaultdict(dict)
DRDQ_eplet = defaultdict(dict)
DR_eplet = defaultdict(dict)
DQ_eplet = defaultdict(dict)
for pair in rand_pairs:
    don, rec = pair.split("+")

    pairing_lines = eplet_imptue[(eplet_imptue['DON_ID'] == don) & (eplet_imptue['REC_ID'] == rec)]

    if pairing_lines.empty is True:
        pairing_lines = eplet_imptue[(eplet_imptue['DON_ID'] == rec) & (eplet_imptue['REC_ID'] == don)]

    pairing_lines = pairing_lines.reset_index(drop=True)

    for row in range(len(pairing_lines)):
        DRDQ_quant = pairing_lines.loc[row, 'ALL_quantity']
        DR_quant = pairing_lines.loc[row, 'DRB_quantity']
        DQ_quant = pairing_lines.loc[row, 'DQ_quantity']
        pairprob = pairing_lines.loc[row, 'PairProb']

        DRDQ_str = pairing_lines.loc[row, 'ALL_details']
        DR_str = pairing_lines.loc[row, 'DRB_details']
        DQ_str = pairing_lines.loc[row, 'DQ_details']

        # Fill the dictionaries for counts
        if DRDQ_quant not in DRDQ_count[pair]:
            DRDQ_count[pair][DRDQ_quant] = pairprob
        else:
            DRDQ_count[pair][DRDQ_quant] = DRDQ_count[pair][DRDQ_quant] + pairprob
        if DR_quant not in DR_count[pair]:
            DR_count[pair][DR_quant] = pairprob
        else:
            DR_count[pair][DR_quant] = DR_count[pair][DR_quant] + pairprob
        if DQ_quant not in DQ_count[pair]:
            DQ_count[pair][DQ_quant] = pairprob
        else:
            DQ_count[pair][DQ_quant] = DQ_count[pair][DQ_quant] + pairprob

        # Fill the dictionaries for eplet strings
        if DRDQ_str not in DRDQ_eplet[pair]:
            DRDQ_eplet[pair][DRDQ_str] = pairprob
        else:
            DRDQ_eplet[pair][DRDQ_str] = DRDQ_eplet[pair][DRDQ_str] + pairprob
        if DR_str not in DR_eplet[pair]:
            DR_eplet[pair][DR_str] = pairprob
        else:
            DR_eplet[pair][DR_str] = DR_eplet[pair][DR_str] + pairprob
        if DQ_str not in DQ_eplet[pair]:
            DQ_eplet[pair][DQ_str] = pairprob
        else:
            DQ_eplet[pair][DQ_str] = DQ_eplet[pair][DQ_str] + pairprob

# Get top counts and top eplet MMs for each level: DRDQ, DR, DQ
top_DRDQ_eplets = pd.DataFrame()
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DRDQ_count, 'quantity', 'ALL')
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DRDQ_eplet, 'details', 'ALL')
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DR_count, 'quantity', 'DRB')
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DR_eplet, 'details', 'DRB')
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DQ_count, 'quantity', 'DQ')
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DQ_eplet, 'details', 'DQ')

top_DRDQ_eplets = top_DRDQ_eplets.reset_index(names=['ID'])
top_DRDQ_eplets = top_DRDQ_eplets.sort_values(by=['ID']).reset_index(drop=True)
eplet_truth = eplet_truth.sort_values(by=['ID']).reset_index(drop=True)

class_ii_headers = ['ALL_quantity', 'ALL_details', 'DRB_quantity', 'DRB_details', 'DQ_quantity', 'DQ_details']

for line in range(len(eplet_truth)):
    for heading in class_ii_headers:
        truth_eplets = eplet_truth.loc[line, heading]
        impute_eplets = top_DRDQ_eplets.loc[line, heading]

        top_DRDQ_eplets.loc[line, heading + '_True'] = neg_prediction(truth_eplets, impute_eplets)

        prob = top_DRDQ_eplets.loc[line, heading + 'Prob']
        threshold = 0.95
        if prob >= threshold:
            top_DRDQ_eplets.loc[line, heading + '_Prediction'] = 1
        else:
            top_DRDQ_eplets.loc[line, heading + '_Prediction'] = 0

for heading in class_ii_headers:
    true_pred = top_DRDQ_eplets[heading + '_True'].to_numpy()  # Actual prediction in terms of 0/1
    probability = top_DRDQ_eplets[heading + 'Prob'].to_numpy()  # Probability of correctness in range of [0,1] terms

    impute_sort = top_DRDQ_eplets.sort_values(by=[heading + 'Prob'])
    sort_probability = impute_sort[heading + 'Prob'].to_numpy()
    sort_true_pred = impute_sort[heading+ '_True'].to_numpy()

    n_pairs = 200

    # Create four points by taking the average in each bin
    n_bins = 4
    prob_q1, prob_q2, prob_q3, prob_q4 = np.array_split(sort_probability, n_bins)
    true_q1, true_q2, true_q3, true_q4 = np.array_split(sort_true_pred, n_bins)
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
    city_block_dst = (abs(probability_avg[0] - true_avg[0]) + abs(probability_avg[1] - true_avg[1]) + abs(
        probability_avg[2] - true_avg[2]) + abs(probability_avg[3] - true_avg[3])) / n_bins

    # Mean Squared Error (MSE) for the bin averages
    mse_bins = np.square((abs(probability_avg[0] - true_avg[0]) + abs(probability_avg[1] - true_avg[1]) + abs(
        probability_avg[2] - true_avg[2]) + abs(probability_avg[3] - true_avg[3]))) / n_bins

    # Create a table exactly like the print statements above to add to the bottom of the calibration plot
    table_data = [["Quantile", "Prob Avg", "True Fraction", "Min Prob", "Max Prob", "Standard Error"],
                  ['Q1', str(round(probability_avg[0], 4)), str(round(true_avg[0], 4)),
                   str(round(min_prob_in_bin[0], 4)), str(round(max_prob_in_bin[0], 4)), str(round(snd_err1, 4))],
                  ['Q2', str(round(probability_avg[1], 4)), str(round(true_avg[1], 4)),
                   str(round(min_prob_in_bin[1], 4)), str(round(max_prob_in_bin[1], 4)), str(round(snd_err2, 4))],
                  ['Q3', str(round(probability_avg[2], 4)), str(round(true_avg[2], 4)),
                   str(round(min_prob_in_bin[2], 4)), str(round(max_prob_in_bin[2], 4)), str(round(snd_err3, 4))],
                  ['Q4', str(round(probability_avg[3], 4)), str(round(true_avg[3], 4)),
                   str(round(min_prob_in_bin[3], 4)), str(round(max_prob_in_bin[3], 4)), str(round(snd_err4, 4))]]

    print('Quantile Statistics for Locus: ', heading)
    print(table_data[0])
    print(table_data[1])
    print(table_data[2])
    print(table_data[3])
    print(table_data[4])

    # Create a bar plot where it shows the distribution of predictions, have to separate it from plot so it does not get added in
    counts, bins, _ = plt.hist(sort_probability, bins=20)
    num_IDs = len(top_DRDQ_eplets)
    fract_counts = counts / num_IDs  # This allows us to get the fraction (* 100 = %) of cases for each count from the histogram

    if heading == 'ALL_quantity':
        TITLE = 'DRDQ eplet counts'
        FILE = 'DRDQ_counts_'
    elif heading == 'ALL_details':
        TITLE = 'DRDQ unique eplets'
        FILE = 'DRDQ_eplets_'
    elif heading == 'DRB_quantity':
        TITLE = 'DR eplet counts'
        FILE = 'DR_counts_'
    elif heading == 'DRB_details':
        TITLE = 'DR unique eplets'
        FILE = 'DR_eplets_'
    elif heading == 'DQ_quantity':
        TITLE = 'DQ eplet counts'
        FILE = 'DQ_counts_'
    elif heading == 'DQ_details':
        TITLE = 'DQ unique eplets'
        FILE = 'DQ_eplets_'

    calibrat_plot = plt.figure(figsize=(8, 8))
    calibrat_plot = plt.errorbar(probability_avg, true_avg, yerr=snd_err, marker='o', linestyle='',
                                 label='True Fraction vs Probability Average for Quantile', color='red', ecolor='black',
                                 capsize=7)
    calibrat_plot = plt.plot([0, 1], linestyle='--', label='Ideal Calibration', color='blue')
    calibrat_plot = plt.xlabel('Mean Predicted Probability for Quantile')
    calibrat_plot = plt.ylabel('Fraction of Predictions Correct')
    # plt.yscale('log')
    # plt.xscale('log')
    # calibrat_plot9loc = plt.suptitle('Calibration Plot and Prediction Probability Distribution for HLA-' + locus + " Locus\n" + 'Brier Score Loss: ' + str(round(brier_loss[locus], 4)) + '\n                        \n')
    calibrat_plot = plt.bar(bins[:-1], fract_counts, width=np.diff(bins), edgecolor='black', color='grey')
    calibrat_plot = plt.xlim(0, 1.05)
    calibrat_plot = plt.ylim(0, 1.05)
    table = plt.table(cellText=table_data)
    table = table.scale(1, 1.5)
    ax = plt.gca()
    ax.xaxis.set_label_position('top')  # have to add x-axis to the top because of the table at the bottom
    ax.xaxis.set_ticks_position('top')
    ax.set_title('Calibration Plot and Prediction Probability Distribution in ' + str(n_pairs) + " pairs for " + TITLE + " \n" + 'Bin Avg MSE: ' +
                 str(round(mse_bins, 4)) + ', Bin Avg City-Block Dist: ' + str(round(city_block_dst, 4)),
                 pad=20)  # Space between x-axis and title
    calibrat_plot = plt.legend()
    calibrat_plot = plt.savefig("Calibration_" + FILE + str(n_pairs) + ".png", bbox_inches='tight')


import pandas as pd
from collections import defaultdict
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss
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

    if truth_typ == impute_typ:
        neg_count = 1
    else:
        neg_count = 0

    return neg_count


n_pairs = '100'
which_impute = 'DRDQ'
truth_filename = which_impute + "_eplet_truth_table" + n_pairs + ".csv"
eplet_truth = pd.read_csv(truth_filename, header=0)

# Get all pairs from the truth table
rand_pairs = {}
for row in range(len(eplet_truth)):
    pair_id = eplet_truth.loc[row, 'ID']

    if pair_id not in rand_pairs:
        rand_pairs[pair_id] = 1
    else:
        rand_pairs[pair_id] = rand_pairs[pair_id] + 1

impute_filename = which_impute + "_eplet_lowres_impute" + n_pairs + ".csv"
eplet_imptue = pd.read_csv(impute_filename, header=0)

# Start with counts and go from there {DON+REC: {Count: prob}}
# Unique eplets do similar thing such as {DON+REC: {Eplet_str: prob}}
DRDQ_count = defaultdict(dict)
DR_count = defaultdict(dict)
DQ_count = defaultdict(dict)
DRDQ_eplet = defaultdict(dict)
DR_eplet = defaultdict(dict)
DQ_eplet = defaultdict(dict)

if which_impute == 'DRDQ':
    which_eplet = 'ALL'
elif which_impute == 'DR':
    which_eplet = 'DRB'
else:
    which_eplet = 'DQ'
for pair in rand_pairs:
    don, rec = pair.split("+")

    pairing_lines = eplet_imptue[(eplet_imptue['DON_ID'] == don) & (eplet_imptue['REC_ID'] == rec)]

    if pairing_lines.empty is True:
        pairing_lines = eplet_imptue[(eplet_imptue['DON_ID'] == rec) & (eplet_imptue['REC_ID'] == don)]

    pairing_lines = pairing_lines.reset_index(drop=True)

    for row in range(len(pairing_lines)):
        DRDQ_quant = pairing_lines.loc[row, which_eplet + '_quantity']
        pairprob = pairing_lines.loc[row, which_impute + '_PairProb']

        DRDQ_str = pairing_lines.loc[row, which_eplet + '_details']

        # Fill the dictionaries for counts
        if DRDQ_quant not in DRDQ_count[pair]:
            DRDQ_count[pair][DRDQ_quant] = pairprob
        else:
            DRDQ_count[pair][DRDQ_quant] = DRDQ_count[pair][DRDQ_quant] + pairprob

        # Fill the dictionaries for eplet strings
        if DRDQ_str not in DRDQ_eplet[pair]:
            DRDQ_eplet[pair][DRDQ_str] = pairprob
        else:
            DRDQ_eplet[pair][DRDQ_str] = DRDQ_eplet[pair][DRDQ_str] + pairprob


# Get top counts and top eplet MMs for each level: DRDQ, DR, DQ
top_DRDQ_eplets = pd.DataFrame()
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DRDQ_count, 'quantity', which_eplet)
top_DRDQ_eplets = top_impute_df(top_DRDQ_eplets, DRDQ_eplet, 'details', which_eplet)

top_DRDQ_eplets = top_DRDQ_eplets.reset_index(names=['ID'])
top_DRDQ_eplets = top_DRDQ_eplets.sort_values(by=['ID']).reset_index(drop=True).fillna('None')
eplet_truth = eplet_truth.sort_values(by=['ID']).reset_index(drop=True).fillna('None')

class_ii_headers = [which_eplet + '_quantity', which_eplet + '_details']

for line in range(len(eplet_truth)):
    for heading in class_ii_headers:
        truth_eplets = eplet_truth.loc[line, heading]
        impute_eplets = top_DRDQ_eplets.loc[line, heading]

        if '_details' in heading:
            truth_eplet_list = truth_eplets.split("_")
            truth_eplet_list.sort()
            impute_eplet_list = impute_eplets.split("_")
            impute_eplet_list.sort()
            top_DRDQ_eplets.loc[line, heading + '_True'] = neg_prediction(truth_eplet_list, impute_eplet_list)
        else:
            truth_eplets = int(truth_eplets)
            impute_eplets = int(impute_eplets)
            top_DRDQ_eplets.loc[line, heading + '_True'] = neg_prediction(truth_eplets, impute_eplets)

        prob = top_DRDQ_eplets.loc[line, heading + 'Prob']
        threshold = 0.9
        if prob >= threshold:
            top_DRDQ_eplets.loc[line, heading + '_Pred'] = 1
        else:
            top_DRDQ_eplets.loc[line, heading + '_Pred'] = 0

for heading in class_ii_headers:
    true_pred = top_DRDQ_eplets[heading + '_True'].to_numpy()  # Actual prediction in terms of 0/1
    probability = top_DRDQ_eplets[heading + 'Prob'].to_numpy()  # Probability of correctness in range of [0,1] terms

    impute_sort = top_DRDQ_eplets.sort_values(by=[heading + 'Prob'])
    sort_probability = impute_sort[heading + 'Prob'].to_numpy()
    sort_true_pred = impute_sort[heading + '_True'].to_numpy()

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

    if heading == which_eplet + '_quantity':
        TITLE = which_impute + ' eplet counts'
        FILE = which_impute + '_counts_'
    elif heading == which_eplet + '_details':
        TITLE = which_impute + ' unique eplets'
        FILE = which_impute + '_eplets_'

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
    ax.set_title('Calibration Plot and Prediction Probability Distribution in ' + n_pairs + " pairs for " + TITLE + " \n" + 'Bin Avg MSE: ' +
                 str(round(mse_bins, 4)) + ', Bin Avg City-Block Dist: ' + str(round(city_block_dst, 4)),
                 pad=20)  # Space between x-axis and title
    calibrat_plot = plt.legend()
    calibrat_plot = plt.savefig("Calibration_" + FILE + n_pairs + ".png", bbox_inches='tight')

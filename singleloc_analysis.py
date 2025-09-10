import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, roc_curve, RocCurveDisplay
import matplotlib.pyplot as plt
import sys

# Compare the imputation to the truth table for single locus unphased genotype analysis
# If both alleles are the same, then that is a positive prediction=1, and if one/both are incorrect then that is a negative prediction=0


# Separate the GLString out into locus pairings
def sep_glstring(file, high_res, glstring, loci=None):
    # If loci not provided, try to infer from columns
    if loci is None:
        # Try to infer loci from GLString columns
        if glstring == 'GLString':
            split_cols = file['GLString'].str.split('^', expand=True)
        else:
            split_cols = file['SLUG_GLString'].str.split('^', expand=True)
        # Use number of columns to infer loci
        default_loci = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
        loci = default_loci[:split_cols.shape[1]]

    # Split GLString into loci columns
    if glstring == 'GLString':
        file[loci] = file['GLString'].str.split('^', expand=True)
        file = file.drop(columns=['GLString'])
    else:
        file[loci] = file['SLUG_GLString'].str.split('^', expand=True)
        file = file.drop(columns=['SLUG_GLString'])
        # Drop columns that may not exist
        for col in ['9loc_GLString', 'HapPair_Prob']:
            if col in file.columns:
                file = file.drop(columns=[col])

    # Split alleles for each locus
    for locus in loci:
        file[[locus + '_1', locus + '_2']] = file[locus].str.split('+', expand=True)
        file = file.drop(columns=[locus])

    # If your genotype truth table is not actually high resolution, then drop extra loci
    if high_res is False:
        keep_loci = ['A', 'C', 'B', 'DRB1', 'DQB1']
        drop_loci = [l for l in loci if l not in keep_loci]
        drop_cols = []
        for locus in drop_loci:
            drop_cols += [locus + '_1', locus + '_2']
        file = file.drop(columns=drop_cols, errors='ignore')

    return file, loci


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

    return neg_count


truth_filename = sys.argv[1]  # 'genotype_truth_table.csv'
impute_filename = sys.argv[2]  # 'lowres_topprob_impute.csv'
truth_table = pd.read_csv(truth_filename, header=0, sep=',', dtype={"ID": str, "GLString": str})
impute = pd.read_csv(impute_filename, header=0, sep=',', dtype={"ID": str})

truth_table = truth_table[truth_table.ID.isin(impute.ID)].reset_index(drop=True)  # Makes sure they are the same length and looking at the same patients
truth_table = truth_table.sort_values(by=['ID']).reset_index(drop=True)  # Sorts the patients, so each patient is in the same row as the imputation rows
impute = impute.sort_values(by=['ID']).reset_index(drop=True)

# Drops duplicates in truth file, if there are any
truth_table = truth_table.drop_duplicates().reset_index(drop=True)
print(len(truth_table))

high_res = True

# Dynamically detect loci from impute file columns
slug_cols = [col for col in impute.columns if col.startswith('GENO_') and col.endswith('_Prob')]
loci = [col.replace('GENO_', '').replace('_Prob', '') for col in slug_cols]

# Parse GLStrings for only present loci
truth_table, loci_truth = sep_glstring(truth_table, high_res, 'GLString', loci)
impute, loci_impute = sep_glstring(impute, high_res, 'SLUG_GLString', loci)

# Use intersection of loci present in both files
loci = [l for l in loci if l in loci_truth and l in loci_impute]

# Filter loci to only those with non-missing allele columns in both truth and impute
filtered_loci = []
for locus in loci:
    truth_cols = [locus + '_1', locus + '_2']
    impute_cols = [locus + '_1', locus + '_2', 'GENO_' + locus + '_Prob']
    # Check columns exist and are not all missing
    if all(col in truth_table.columns for col in truth_cols) and all(col in impute.columns for col in impute_cols):
        # Check for at least one non-missing value
        if not (truth_table[truth_cols].isnull().all().any() or impute[impute_cols].isnull().all().any()):
            filtered_loci.append(locus)
loci = filtered_loci

for line in range(len(truth_table)):
    for locus in loci:
        truth1 = truth_table.loc[line, locus + '_1']
        truth2 = truth_table.loc[line, locus + '_2']
        (truth1, truth2) = sorted([truth1, truth2])

        impute1 = impute.loc[line, locus + '_1']
        impute2 = impute.loc[line, locus + '_2']
        (impute1, impute2) = sorted([impute1, impute2])

        impute.loc[line, locus + '_True'] = neg_prediction(truth1, truth2, impute1, impute2)

        probability = impute.loc[line, 'GENO_' + locus + '_Prob']
        threshold = 0.90
        if probability >= threshold:
            impute.loc[line, locus + '_Prob'] = 1
        else:
            impute.loc[line, locus + '_Prob'] = 0


confusion_mat = pd.DataFrame()
brier_loss = {}
roc_auc = {}
for locus in loci:
    # Need everything to be in numpy array for sklearn metrics
    probs = impute[locus + '_Prob'].to_numpy()
    probability = impute['GENO_' + locus + '_Prob'].to_numpy()
    true_pred = impute[locus + '_True'].to_numpy()

    print('Classification Report for Locus ', locus)
    cm = confusion_matrix(true_pred, probs, labels=[0, 1])
    if cm.shape == (2, 2):
        tn, fp, fn, tp = cm.ravel()
    else:
        tn = fp = fn = tp = 0
        if cm.shape == (1, 1):
            if np.unique(true_pred)[0] == 1:
                tp = cm[0, 0]
            else:
                tn = cm[0, 0]
        elif cm.shape == (1, 2):
            tn, fp = cm[0, 0], cm[0, 1]
        elif cm.shape == (2, 1):
            tn, fn = cm[0, 0], cm[1, 0]
    print('TP: ' + str(tp) + ', TN: ' + str(tn) + ', FP: ' + str(fp) + ', FN: ' + str(fn))

    confusion_mat_line = pd.DataFrame({'TN': tn, 'FP': fp, 'FN': fn, 'TP': tp}, index=[locus])
    # confusion_mat = pd.concat([confusion_mat, confusion_mat_line])

    # Classification Report
    print(classification_report(true_pred, probs, labels=[0, 1]))

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

    # Create how many quantiles you want by taking the average in each bin
    n_bins = int(sys.argv[3])
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

    print('Quantile Statistics for Locus: ', locus)
    for quantiles in range(0, n_bins + 1):
        print(table_data[quantiles])

    # Create a bar plot where it shows the distribution of predictions, have to separate it from plot so it does not get added in
    counts, bins, _ = plt.hist(sort_probability, bins=20)
    num_IDs = len(impute)
    fract_counts = counts / num_IDs  # This allows us to get the fraction (* 100 = %) of cases for each count from the histogram

    calibrat_plot = plt.figure(figsize=(8, 8))
    calibrat_plot = plt.errorbar(probability_avg, true_avg, yerr=snd_err, marker='o', linestyle='', label='True Fraction vs Probability Average for Quantile', color='red', ecolor='black', capsize=7)
    calibrat_plot = plt.plot([0,1], linestyle='--', label='Ideal Calibration', color='blue')
    calibrat_plot = plt.xlabel('Mean Predicted Probability for Quantile')
    calibrat_plot = plt.ylabel('Fraction of Predictions Correct')
    # plt.yscale('log')
    # plt.xscale('log')
    # calibrat_plot9loc = plt.suptitle('Calibration Plot and Prediction Probability Distribution for HLA-' + locus + " Locus\n" + 'Brier Score Loss: ' + str(round(brier_loss[locus], 4)) + '\n                        \n')
    calibrat_plot = plt.bar(bins[:-1], fract_counts, width=np.diff(bins), edgecolor='black', color='grey')
    calibrat_plot = plt.xlim(0,1.05)
    calibrat_plot = plt.ylim(0,1.05)
    table = plt.table(cellText=table_data)
    table = table.scale(1, 1.5)
    ax = plt.gca()
    ax.xaxis.set_label_position('top')  # have to add x-axis to the top because of the table at the bottom
    ax.xaxis.set_ticks_position('top')
    ax.set_title('Calibration Plot and Prediction Probability Distribution for HLA-' + locus + " Locus\n" +
                 'Brier: ' + str(round(brier_loss[locus], 4)) + ', Bin Avg MSE: ' +
                 str(round(mse_bins, 4)) + ', Bin Avg City-Block Dist: ' + str(round(city_block_dst, 4)), pad=20)  # Space between x-axis and title
    calibrat_plot = plt.legend()
    calibrat_plot = plt.savefig("Calibration_" + locus + ".png", bbox_inches='tight')
    # plt.show()

plt.figure(0)
for locus in loci:
    # Display ROC plots for each loci
    probability = impute['GENO_' + locus + '_Prob'].to_numpy()  # Probability of correctness in [0,1] terms
    true_pred = impute[locus + '_True'].to_numpy()  # Actual prediction in terms of 0/1
    oneloc_ROC = RocCurveDisplay.from_predictions(true_pred, probability, plot_chance_level=True)
    oneloc_ROC = plt.title('ROC Curve and AUC for ' + locus)
    oneloc_ROC = plt.savefig("ROC_" + locus + ".png", bbox_inches='tight')
    # # plt.show()

for locus in loci:
    # Display ROC plots for each loci on one plot
    probability = impute['GENO_' + locus + '_Prob'].to_numpy()  # Probability of correctness in [0,1] terms
    true_pred = impute[locus + '_True'].to_numpy()
    fpr, tpr, thresholds = roc_curve(true_pred, probability)
    # Calculate Area under the curve to display on the plot
    auc = roc_auc_score(true_pred, probability)
    # Now, plot the computed values
    roc_plot = plt.plot(fpr, tpr, label='%s (AUC = %0.2f)' % (locus, auc))
    # Custom settings for the plot
roc_plot = plt.plot([0, 1], 'r--', color='black')
roc_plot = plt.xlabel('Specificity (False Positive Rate)')
roc_plot = plt.ylabel('Sensitivity (True Positive Rate)')
roc_plot = plt.title('ROC-AUC Curves for each Locus in SLUG Analysis')
roc_plot = plt.legend(loc="lower right")
roc_plot = plt.savefig("ROC_AUC_9loc.png", bbox_inches='tight')

bf = pd.DataFrame({'Brier_Loss_Score': brier_loss}, index=brier_loss.keys())
ra = pd.DataFrame({'ROC-AUC': roc_auc}, index=roc_auc.keys())
confusion_mat = pd.concat([confusion_mat, bf], axis=1)
confusion_mat = pd.concat([confusion_mat, ra], axis=1)
print(confusion_mat)

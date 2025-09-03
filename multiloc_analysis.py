import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, RocCurveDisplay
import sys
import matplotlib.pyplot as plt

# Compare the top imputation to the truth table for multiple multilocus unphased genotype analysis
# At 9-loci, 7-loci, DR-DQ, and Class I levels
# If the GLString is the same, then it is a correct match


# Use this definition to make negative(0) and positive(1) predictions
def neg_prediction(truth_typ, impute_typ, gl_col, true_col):
    for row in range(len(impute_typ)):
        impute_string = impute_typ.loc[row, gl_col]
        truth_string = truth_typ.loc[row, gl_col]
        if truth_string == "MISSING" or impute_string == "MISSING":
            impute_typ.loc[row, true_col] = "NA"
        elif truth_string == impute_string:
            impute_typ.loc[row, true_col] = 1
        else:
            impute_typ.loc[row, true_col] = 0
    return impute_typ


# Separate GLString for truth tables for calculating predictions for anything other than 9-loci
def sep_glstring(file, loci, source_col='GLString'):
    # Split the source GLString column into loci columns
    file[loci] = file[source_col].str.split('^', expand=True)
    # Rebuild a GLString for just the selected loci
    file['GLString_temp'] = file[loci].agg('^'.join, axis=1)
    # Drop the loci columns
    for locus in loci:
        file = file.drop(columns=[locus])
    return file


truth_filename = sys.argv[1]
impute_filename = sys.argv[2]
num_bins = int(sys.argv[3])
truth_table = pd.read_csv(truth_filename, header=0, dtype={"ID": str})
impute = pd.read_csv(impute_filename, header=0, dtype={"ID": str})

truth_table = truth_table[truth_table.ID.isin(impute.ID)].reset_index(drop=True)
truth_table = truth_table.sort_values(by=['ID']).reset_index(drop=True)
impute = impute.sort_values(by=['ID']).reset_index(drop=True)
truth_table = truth_table.drop_duplicates().reset_index(drop=True)
print(len(truth_table))

# Dynamically detect loci from imputation file GLString columns
gl_col = None
prob_col = []

if 'SLUG_GLString' in impute.columns:
    slug_example = impute.loc[0, 'SLUG_GLString']
    loci = [allele.split('*')[0] for allele in slug_example.split('^')]
    gl_col = 'SLUG_GLString'
    prob_col = [f'GENO_{locus}_Prob' for locus in loci]
elif '9loc_GLString' in impute.columns:
    slug_example = impute.loc[0, '9loc_GLString']
    loci = [allele.split('*')[0] for allele in slug_example.split('^')]
    gl_col = '9loc_GLString'
    prob_col = ['HapPair_Prob']
elif 'GLString' in impute.columns:
    slug_example = impute.loc[0, 'GLString']
    loci = [allele.split('*')[0] for allele in slug_example.split('^')]
    gl_col = 'GLString'
    prob_col = [col for col in impute.columns if col.endswith('_Prob') and any(locus in col for locus in loci)]
else:
    raise Exception("No GLString column found in imputation file.")

# Build truth GLString for detected loci (always from 'GLString' column)
truth_table = sep_glstring(truth_table, loci, source_col='GLString')
# Rename the constructed column to match impute's GLString column for comparison
truth_table[gl_col] = truth_table['GLString_temp']
truth_table = truth_table.drop(columns=['GLString_temp'])

# Prepare impute DataFrame for analysis
impute = impute[['ID', gl_col] + prob_col]

# Run prediction and analysis for detected loci set
impute = neg_prediction(truth_table, impute, gl_col, gl_col.replace('GLString', 'True'))
print(f'Positive and Negative Predictions for loci {",".join(loci)}:\n', impute[gl_col.replace('GLString', 'True')].value_counts())


# Create Calibration Plots for each case by creating four points
def calibration_plot(impute_typ, gl_col, prob_col, n_bins):
    probability = impute_typ[prob_col[0]].to_numpy()
    impute_sort = impute_typ.sort_values(by=[prob_col[0]])
    sort_probability = impute_sort[prob_col[0]].to_numpy()
    true_pred = impute_typ[gl_col.replace('GLString', 'True')].to_numpy()
    threshold = 0.5
    brier_loss = brier_score_loss(true_pred, probability > threshold)
    sort_true_pred = impute_sort[gl_col.replace('GLString', 'True')].to_numpy()

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

    print('Quantile Statistics for loci ' + '^'.join(loci) + ':')
    for quantiles in range(0, n_bins + 1):
        print(table_data[quantiles])

    # Create a bar plot where it shows the distribution of predictions, have to separate it from plot so it does not get added in
    counts, bins, _ = plt.hist(sort_probability, bins=20)
    num_IDs = len(impute_typ)
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
    title = f'Calibration Plot and Prediction Probability Distribution for loci {"^".join(loci)}\n'
    ax.set_title(title +
                 'Brier: ' + str(round(brier_loss, 4)) + ', Bin Avg MSE: ' +
                 str(round(mse_bins, 4)) + ', Bin Avg City-Block Dist: ' + str(round(city_block_dst, 4)), pad=20)  # Space between x-axis and title
    calibrat_plot = plt.legend()

    return calibrat_plot

calibrat_plot = calibration_plot(impute, gl_col, prob_col, num_bins)
plt.savefig(f"Calibration_{'_'.join(loci)}.png", bbox_inches='tight')

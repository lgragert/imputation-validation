"""Plotting functions for calibration and ROC curves."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, classification_report, brier_score_loss, roc_auc_score, roc_curve

class CalibrationPlotter:
    """Create calibration and ROC plots for HLA imputation validation."""

    def __init__(self, n_bins: int = 4):
        self.n_bins = n_bins

    def calibration_plot(self, analysis_results, locus, threshold=0.9, title=None, save_path=None):
        # Extract relevant columns
        preds = analysis_results['y_pred'].to_numpy()  # Binary prediction (0/1)
        probability = analysis_results['y_prob'].to_numpy() # Predicted probability
        true_pred = analysis_results['y_true'].to_numpy() # True labels (0/1)
        id_col = analysis_results['ID']

        # Confusion Matrix per Locus
        print('Classification Report for Locus', locus)
        cm = confusion_matrix(true_pred, preds)
        if cm.size == 4:
            tn, fp, fn, tp = cm.ravel()
        else:
            tn, fp, fn, tp = 0, 0, 0, 0  # fallback for single-class
        print(f'TP: {tp}, TN: {tn}, FP: {fp}, FN: {fn}')
        confusion_mat_line = pd.DataFrame({'TN': tn, 'FP': fp, 'FN': fn, 'TP': tp}, index=[locus])
        print(classification_report(true_pred, preds))

        # Brier Score Loss
        brier = brier_score_loss(true_pred, probability > threshold)
        print('Brier Score Loss:', brier)

        # ROC-AUC Score
        roc_auc = roc_auc_score(true_pred, probability)
        print('ROC-AUC Score:', roc_auc)

        # Sort the Dataset based on probabilities
        impute_sort = analysis_results.sort_values(by=['y_prob'])
        sort_probability = impute_sort['y_prob'].to_numpy()
        sort_true_pred = impute_sort['y_true'].to_numpy()

        # Quantile statistics
        n_bins = self.n_bins
        split_prob = np.array_split(sort_probability, n_bins)
        split_true = np.array_split(sort_true_pred, n_bins)
        probability_avg = [quantiles.mean() for quantiles in split_prob]
        min_prob_in_bin = [np.min(quantiles) for quantiles in split_prob]
        max_prob_in_bin = [np.max(quantiles) for quantiles in split_prob]
        true_avg = [quantiles.mean() for quantiles in split_true]
        snd_err = [(np.std(quantiles, ddof=1) / np.sqrt(np.size(quantiles))) for quantiles in split_prob]
        city_block_dst = sum(abs(np.array(probability_avg) - np.array(true_avg))) / n_bins
        mse_bins = np.square(sum(abs(np.array(probability_avg) - np.array(true_avg)))) / n_bins

        # Table data
        table_data = [
            ['Q' + str(q + 1), str(round(probability_avg[q], 4)), str(round(true_avg[q], 4)),
             str(round(min_prob_in_bin[q], 4)), str(round(max_prob_in_bin[q], 4)),
             str(round(snd_err[q], 4))]
            for q in range(n_bins)
        ]
        heading = ["Quantile", "Prob Avg", "True Fraction", "Min Prob", "Max Prob", "Standard Error"]
        table_data.insert(0, heading)

        print('Quantile Statistics for Locus:', locus)
        for quantiles in range(0, n_bins + 1):
            print(table_data[quantiles])

        # Histogram
        counts, bins = np.histogram(sort_probability, bins=20)
        num_IDs = len(analysis_results)
        fract_counts = counts / num_IDs

        # Plot
        plt.figure(figsize=(8, 8))
        plt.errorbar(probability_avg, true_avg, yerr=snd_err, marker='o', linestyle='', label='True Fraction vs Probability Average for Quantile', color='red', ecolor='black', capsize=7)
        plt.plot([0, 1], [0, 1], linestyle='--', label='Ideal Calibration', color='blue')
        plt.bar(bins[:-1], fract_counts, width=np.diff(bins), edgecolor='black', color='grey', alpha=0.5)
        plt.xlim(0, 1.05)
        plt.ylim(0, 1.05)
        plt.xlabel('Mean Predicted Probability for Quantile')
        plt.ylabel('Fraction of Predictions Correct')
        table = plt.table(cellText=table_data, loc='bottom', cellLoc='center')
        table.scale(1, 1.5)
        ax = plt.gca()
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_ticks_position('top')
        ax.set_title(
            f'{title} \n'
            f'Brier: {round(brier, 4)}, Bin Avg MSE: {round(mse_bins, 4)}, Bin Avg City-Block Dist: {round(city_block_dst, 4)}',
            pad=20
        )
        plt.legend()
        if save_path:
            plt.savefig(save_path, bbox_inches='tight', dpi=300)
        return plt.gcf(), confusion_mat_line, brier, roc_auc

    def roc_plot(self, y_true, y_prob, title, save_path=None):
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        auc = roc_auc_score(y_true, y_prob)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)
        if save_path:
            fig.savefig(save_path, bbox_inches='tight', dpi=300)
        return fig

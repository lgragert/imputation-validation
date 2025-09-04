"""Plotting functions for calibration and ROC curves."""

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score, brier_score_loss


class CalibrationPlotter:
    """Create calibration and ROC plots for HLA imputation validation."""

    def __init__(self, n_bins: int = 10):
        self.n_bins = n_bins

    def calibration_plot(self, y_true, y_prob, title, save_path=None):
        # Sort by probability
        sorted_idx = np.argsort(y_prob)
        sorted_prob = y_prob[sorted_idx]
        sorted_true = y_true[sorted_idx]

        split_prob = np.array_split(sorted_prob, self.n_bins)
        split_true = np.array_split(sorted_true, self.n_bins)
        probability_avg = [quantiles.mean() for quantiles in split_prob]
        true_avg = [quantiles.mean() for quantiles in split_true]
        snd_err = [(np.std(quantiles, ddof=1) / np.sqrt(np.size(quantiles))) for quantiles in split_prob]

        brier_loss = brier_score_loss(y_true, y_prob)
        city_block_dst = np.mean([abs(p - t) for p, t in zip(probability_avg, true_avg)])

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.errorbar(probability_avg, true_avg, yerr=snd_err, marker='o', linestyle='', color='red', capsize=7)
        ax.plot([0, 1], [0, 1], linestyle='--', color='blue')
        counts, bins = np.histogram(y_prob, bins=20, range=(0, 1))
        bin_width = bins[1] - bins[0]
        ax.bar(bins[:-1], counts / len(y_prob), width=bin_width, alpha=0.3, color='gray', edgecolor='black')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Mean Predicted Probability')
        ax.set_ylabel('Fraction of Positives')
        ax.set_title(f'{title}\nBrier Score: {brier_loss:.4f}, City Block Distance: {city_block_dst:.4f}')
        ax.legend(['Calibration', 'Perfect Calibration'])
        if save_path:
            fig.savefig(save_path, bbox_inches='tight', dpi=300)
        return fig

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
"""Plotting functions for calibration and ROC curves."""

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score


class CalibrationPlotter:
    """Create calibration plots for HLA imputation validation."""

    def __init__(self, n_bins: int = 10):
        """Initialize with number of bins for calibration."""
        self.n_bins = n_bins

    def create_calibration_plot(self, y_true: np.ndarray, y_prob: np.ndarray,
                              title: str, save_path: str = None) -> plt.Figure:
        """Create a calibration plot."""
        # Sort by probability
        sorted_indices = np.argsort(y_prob)
        sorted_prob = y_prob[sorted_indices]
        sorted_true = y_true[sorted_indices]

        # Create bins
        bin_boundaries = np.linspace(0, len(sorted_prob) - 1, self.n_bins + 1)
        bin_lowers = bin_boundaries[:-1].astype(int)
        bin_uppers = bin_boundaries[1:].astype(int)

        bin_centers = []
        bin_true_fractions = []
        bin_errors = []

        for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
            bin_probs = sorted_prob[bin_lower:bin_upper]
            bin_trues = sorted_true[bin_lower:bin_upper]

            if len(bin_probs) > 0:
                bin_centers.append(np.mean(bin_probs))
                bin_true_fractions.append(np.mean(bin_trues))
                bin_errors.append(np.std(bin_probs) / np.sqrt(len(bin_probs)))

        # Create figure
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot calibration curve
        ax.errorbar(bin_centers, bin_true_fractions, yerr=bin_errors,
                   marker='o', linestyle='', color='red', capsize=5,
                   label='Observed')

        # Plot perfect calibration line
        ax.plot([0, 1], [0, 1], linestyle='--', color='blue',
                label='Perfect Calibration')

        # Plot probability distribution
        counts, bins = np.histogram(y_prob, bins=20, range=(0, 1))
        bin_width = bins[1] - bins[0]
        ax.bar(bins[:-1], counts / len(y_prob), width=bin_width,
               alpha=0.3, color='gray', edgecolor='black')

        # Calculate metrics
        brier_score = np.mean((y_prob - y_true) ** 2)
        city_block_dist = np.mean([abs(pred - true) for pred, true
                                  in zip(bin_centers, bin_true_fractions)])

        ax.set_xlabel('Mean Predicted Probability')
        ax.set_ylabel('Fraction of Positives')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(f'{title}\nBrier Score: {brier_score:.4f}, '
                    f'City Block Distance: {city_block_dist:.4f}')
        ax.legend()

        if save_path:
            fig.savefig(save_path, bbox_inches='tight', dpi=300)

        return fig

    def create_roc_plot(self, y_true: np.ndarray, y_prob: np.ndarray,
                       label: str = None, save_path: str = None) -> plt.Figure:
        """Create ROC curve plot."""
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        auc = roc_auc_score(y_true, y_prob)

        fig, ax = plt.subplots(figsize=(8, 6))

        if label:
            ax.plot(fpr, tpr, label=f'{label} (AUC = {auc:.3f})')
        else:
            ax.plot(fpr, tpr, label=f'AUC = {auc:.3f}')

        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title('ROC Curve')
        ax.legend()
        ax.grid(True, alpha=0.3)

        if save_path:
            fig.savefig(save_path, bbox_inches='tight', dpi=300)

        return fig
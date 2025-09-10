"""Core functionality for HLA imputation validation."""

import pandas as pd
import numpy as np
from sklearn.metrics import (
    confusion_matrix, classification_report, brier_score_loss,
    roc_auc_score, roc_curve
)


class BaseAnalysis:
    """Base class for HLA imputation analysis."""

    def __init__(self, truth_file: str, impute_file: str):
        """Initialize with truth and imputation data files."""
        self.truth_data = pd.read_csv(truth_file, dtype={"ID": str})
        self.impute_data = pd.read_csv(impute_file, dtype={"ID": str})
        self._align_data()

    def _align_data(self):
        """Align truth and imputation data by ID."""
        # Keep only common IDs
        common_ids = self.truth_data['ID'].isin(self.impute_data['ID'])
        self.truth_data = self.truth_data[common_ids].reset_index(drop=True)

        # Sort both by ID
        self.truth_data = self.truth_data.sort_values('ID').reset_index(drop=True)
        self.impute_data = self.impute_data.sort_values('ID').reset_index(drop=True)

        # Remove duplicates
        self.truth_data = self.truth_data.drop_duplicates().reset_index(drop=True)

    def calculate_metrics(self, y_true: np.ndarray, y_pred: np.ndarray,
                         y_prob: np.ndarray) -> dict:
        """Calculate standard classification metrics."""
        threshold = 0.5

        # Confusion matrix
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

        # Brier score
        brier_score = brier_score_loss(y_true, y_prob > threshold)

        # ROC-AUC
        roc_auc = roc_auc_score(y_true, y_prob)

        return {
            'confusion_matrix': {'TP': tp, 'TN': tn, 'FP': fp, 'FN': fn},
            'brier_score': brier_score,
            'roc_auc': roc_auc
        }
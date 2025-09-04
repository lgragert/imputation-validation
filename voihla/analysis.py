"""Analysis classes for HLA imputation validation."""

import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, brier_score_loss, roc_auc_score


class SingleLocusAnalysis:
    """Single locus analysis for arbitrary loci."""

    def __init__(self, truth_df: pd.DataFrame, impute_df: pd.DataFrame):
        self.truth_df = truth_df.copy()
        self.impute_df = impute_df.copy()
        self.loci = self._detect_loci()

    def _detect_loci(self):
        # Detect loci from SLUG_GLString column
        slug_example = self.impute_df.loc[0, 'SLUG_GLString']
        loci = [allele.split('*')[0] for allele in slug_example.split('^')]
        return loci

    def run(self):
        results = {}
        for locus in self.loci:
            prob_col = f'GENO_{locus}_Prob'
            if prob_col not in self.impute_df.columns:
                continue
            y_true = []
            for i in range(len(self.truth_df)):
                truth_alleles = sorted(self.truth_df.loc[i, 'GLString'].split('^')[self.loci.index(locus)].split('+'))
                impute_alleles = sorted(self.impute_df.loc[i, 'SLUG_GLString'].split('^')[self.loci.index(locus)].split('+'))
                y_true.append(1 if truth_alleles == impute_alleles else 0)
            y_true = np.array(y_true)
            y_prob = self.impute_df[prob_col].values
            y_pred = (y_prob > 0.5).astype(int)
            results[locus] = self._metrics(y_true, y_pred, y_prob)
        return results

    def get_results_df(self):
        """Return a dict of DataFrames with columns y_true, y_pred, y_prob for each locus."""
        results_df = {}
        for locus in self.loci:
            prob_col = f'GENO_{locus}_Prob'
            if prob_col not in self.impute_df.columns:
                continue
            y_true = []
            for i in range(len(self.truth_df)):
                truth_alleles = sorted(self.truth_df.loc[i, 'GLString'].split('^')[self.loci.index(locus)].split('+'))
                impute_alleles = sorted(self.impute_df.loc[i, 'SLUG_GLString'].split('^')[self.loci.index(locus)].split('+'))
                y_true.append(1 if truth_alleles == impute_alleles else 0)
            y_true = np.array(y_true)
            y_prob = self.impute_df[prob_col].values
            y_pred = (y_prob > 0.5).astype(int)
            df = pd.DataFrame({'y_true': y_true, 'y_pred': y_pred, 'y_prob': y_prob})
            results_df[locus] = df
        return results_df

    def _metrics(self, y_true, y_pred, y_prob):
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
        brier = brier_score_loss(y_true, y_prob)
        roc_auc = roc_auc_score(y_true, y_prob)
        return {'TP': tp, 'TN': tn, 'FP': fp, 'FN': fn, 'brier': brier, 'roc_auc': roc_auc}

    def pred_incorrect_high_prob(self, locus, threshold: float = 0.9):
        """Cases above threshold where prediction was incorrect for a given locus."""
        return self._threshold_pred_df(locus, threshold, correct=False, high=True)

    def pred_incorrect_low_prob(self, locus, threshold: float = 0.9):
        """Cases below threshold where prediction was incorrect for a given locus."""
        return self._threshold_pred_df(locus, threshold, correct=False, high=False)

    def pred_correct_high_prob(self, locus, threshold: float = 0.9):
        """Cases above threshold where prediction was correct for a given locus."""
        return self._threshold_pred_df(locus, threshold, correct=True, high=True)

    def pred_correct_low_prob(self, locus, threshold: float = 0.9):
        """Cases below threshold where prediction was correct for a given locus."""
        return self._threshold_pred_df(locus, threshold, correct=True, high=False)

    def _threshold_pred_df(self, locus, threshold, correct, high):
        prob_col = f'GENO_{locus}_Prob'
        if prob_col not in self.impute_df.columns:
            raise ValueError(f"{prob_col} not found in imputation DataFrame.")
        id_col = self.impute_df['ID']
        y_true = []
        imputed_loci = []
        true_loci = []
        for i in range(len(self.truth_df)):
            truth_alleles = sorted(self.truth_df.loc[i, 'GLString'].split('^')[self.loci.index(locus)].split('+'))
            impute_alleles = sorted(self.impute_df.loc[i, 'SLUG_GLString'].split('^')[self.loci.index(locus)].split('+'))
            y_true.append(1 if truth_alleles == impute_alleles else 0)
            imputed_loci.append('+'.join(impute_alleles))
            true_loci.append('+'.join(truth_alleles))
        imp_prob = self.impute_df[prob_col]
        concat_pred = pd.DataFrame({
            'ID': id_col,
            'imp_prob': imp_prob,
            'true_pred': y_true,
            'imputed_loci': imputed_loci,
            'true_loci': true_loci
        })
        if high:
            df = concat_pred[concat_pred['imp_prob'] >= threshold]
        else:
            df = concat_pred[concat_pred['imp_prob'] < threshold]
        if correct:
            return df[df['true_pred'] == 1]
        else:
            return df[df['true_pred'] == 0]


class MultiLocusAnalysis:
    """Multi-locus analysis for arbitrary sets of loci."""

    def __init__(self, truth_df: pd.DataFrame, impute_df: pd.DataFrame):
        self.truth_df = truth_df.copy()
        self.impute_df = impute_df.copy()
        self.analysis_types = self._detect_analysis_types()

    def _detect_analysis_types(self):
        # Detect available GLString columns for analysis
        types = []
        for col in self.impute_df.columns:
            if col.endswith('_GLString'):
                types.append(col.replace('_GLString', ''))
        return types

    def run(self):
        results = {}
        for analysis_type in self.analysis_types:
            gl_col = f'{analysis_type}_GLString'
            prob_col = f'{analysis_type}_Prob' if f'{analysis_type}_Prob' in self.impute_df.columns else 'HapPair_Prob'
            if gl_col not in self.impute_df.columns or prob_col not in self.impute_df.columns:
                continue
            y_true = []
            for i in range(len(self.truth_df)):
                truth_gl = self.truth_df.loc[i, 'GLString']
                impute_gl = self.impute_df.loc[i, gl_col]
                y_true.append(1 if truth_gl == impute_gl else 0)
            y_true = np.array(y_true)
            y_prob = self.impute_df[prob_col].values
            y_pred = (y_prob > 0.5).astype(int)
            results[analysis_type] = self._metrics(y_true, y_pred, y_prob)
        return results

    def get_results_df(self):
        """Return a dict of DataFrames with columns y_true, y_pred, y_prob for each analysis type."""
        results_df = {}
        for analysis_type in self.analysis_types:
            gl_col = f'{analysis_type}_GLString'
            prob_col = f'{analysis_type}_Prob' if f'{analysis_type}_Prob' in self.impute_df.columns else 'HapPair_Prob'
            if gl_col not in self.impute_df.columns or prob_col not in self.impute_df.columns:
                continue
            y_true = []
            for i in range(len(self.truth_df)):
                truth_gl = self.truth_df.loc[i, 'GLString']
                impute_gl = self.impute_df.loc[i, gl_col]
                y_true.append(1 if truth_gl == impute_gl else 0)
            y_true = np.array(y_true)
            y_prob = self.impute_df[prob_col].values
            y_pred = (y_prob > 0.5).astype(int)
            df = pd.DataFrame({'y_true': y_true, 'y_pred': y_pred, 'y_prob': y_prob})
            results_df[analysis_type] = df
        return results_df

    def _metrics(self, y_true, y_pred, y_prob):
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
        brier = brier_score_loss(y_true, y_prob)
        roc_auc = roc_auc_score(y_true, y_prob)
        return {'TP': tp, 'TN': tn, 'FP': fp, 'FN': fn, 'brier': brier, 'roc_auc': roc_auc}

    def pred_incorrect_high_prob(self, analysis_type, threshold: float = 0.9):
        """Cases above threshold where prediction was incorrect for a given analysis type."""
        return self._threshold_pred_df(analysis_type, threshold, correct=False, high=True)

    def pred_incorrect_low_prob(self, analysis_type, threshold: float = 0.9):
        """Cases below threshold where prediction was incorrect for a given analysis type."""
        return self._threshold_pred_df(analysis_type, threshold, correct=False, high=False)

    def pred_correct_high_prob(self, analysis_type, threshold: float = 0.9):
        """Cases above threshold where prediction was correct for a given analysis type."""
        return self._threshold_pred_df(analysis_type, threshold, correct=True, high=True)

    def pred_correct_low_prob(self, analysis_type, threshold: float = 0.9):
        """Cases below threshold where prediction was correct for a given analysis type."""
        return self._threshold_pred_df(analysis_type, threshold, correct=True, high=False)

    def _threshold_pred_df(self, analysis_type, threshold, correct, high):
        gl_col = f'{analysis_type}_GLString'
        prob_col = f'{analysis_type}_Prob' if f'{analysis_type}_Prob' in self.impute_df.columns else 'HapPair_Prob'
        if gl_col not in self.impute_df.columns or prob_col not in self.impute_df.columns:
            raise ValueError(f"{gl_col} or {prob_col} not found in imputation DataFrame.")
        id_col = self.impute_df['ID']
        y_true = []
        imputed_loci = []
        true_loci = []
        for i in range(len(self.truth_df)):
            truth_gl = self.truth_df.loc[i, 'GLString']
            impute_gl = self.impute_df.loc[i, gl_col]
            y_true.append(1 if truth_gl == impute_gl else 0)
            imputed_loci.append(impute_gl)
            true_loci.append(truth_gl)
        imp_prob = self.impute_df[prob_col]
        concat_pred = pd.DataFrame({
            'ID': id_col,
            'imp_prob': imp_prob,
            'true_pred': y_true,
            'imputed_loci': imputed_loci,
            'true_loci': true_loci
        })
        if high:
            df = concat_pred[concat_pred['imp_prob'] >= threshold]
        else:
            df = concat_pred[concat_pred['imp_prob'] < threshold]
        if correct:
            return df[df['true_pred'] == 1]
        else:
            return df[df['true_pred'] == 0]

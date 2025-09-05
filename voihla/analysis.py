"""Analysis classes for HLA imputation validation."""

import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, brier_score_loss, roc_auc_score
from .core import BaseAnalysis


class SingleLocusAnalysis:
    """Single locus analysis for arbitrary loci."""

    def __init__(self, truth_df: pd.DataFrame, impute_df: pd.DataFrame):
        self.truth_df = truth_df.copy()
        self.impute_df = impute_df.copy()
        # Use BaseAnalysis._align_data to align truth and impute DataFrames
        base = BaseAnalysis.__new__(BaseAnalysis)
        base.truth_data = self.truth_df
        base.impute_data = self.impute_df
        base._align_data()
        self.truth_df = base.truth_data
        self.impute_df = base.impute_data
        self.loci = self._detect_loci()

    def _detect_loci(self):
        # Detect loci from SLUG_GLString column and map DRB3/4/5/X to DRB345
        slug_example = self.impute_df.loc[0, 'SLUG_GLString']
        loci = []
        for allele in slug_example.split('^'):
            locus = allele.split('*')[0]
            if locus.startswith('DRB3') or locus.startswith('DRB4') or locus.startswith('DRB5') or locus.startswith('DRBX'):
                locus = 'DRB345'
            loci.append(locus)
        return loci

    def run(self):
        results = {}
        for locus in self.loci:
            prob_col = f'GENO_{locus}_Prob'
            if prob_col not in self.impute_df.columns:
                continue
            y_true = []
            for i in range(len(self.truth_df)):
                truth_allele = self.truth_df.loc[i, 'GLString'].split('^')[self.loci.index(locus)]
                impute_allele = self.impute_df.loc[i, 'SLUG_GLString'].split('^')[self.loci.index(locus)]
                truth_locus = truth_allele.split('*')[0]
                impute_locus = impute_allele.split('*')[0]
                # Only map DRB3/4/5/X to DRB345
                if truth_locus.startswith('DRB3') or truth_locus.startswith('DRB4') or truth_locus.startswith('DRB5') or truth_locus.startswith('DRBX'):
                    truth_allele = truth_allele.replace(truth_locus + '*', 'DRB345*')
                if impute_locus.startswith('DRB3') or impute_locus.startswith('DRB4') or impute_locus.startswith('DRB5') or impute_locus.startswith('DRBX'):
                    impute_allele = impute_allele.replace(impute_locus + '*', 'DRB345*')
                truth_alleles = sorted(truth_allele.split('+'))
                impute_alleles = sorted(impute_allele.split('+'))
                y_true.append(1 if truth_alleles == impute_alleles else 0)
            y_true = np.array(y_true)
            y_prob = self.impute_df[prob_col].values
            y_pred = (y_prob > 0.5).astype(int)
            results[locus] = self._metrics(y_true, y_pred, y_prob)
        return results

    def get_results_df(self):
        """Return a dict of DataFrames with columns ID, y_true, y_pred, y_prob for each locus."""
        results_df = {}
        for locus in self.loci:
            prob_col = f'GENO_{locus}_Prob'
            if prob_col not in self.impute_df.columns:
                continue
            y_true = []
            ids = self.impute_df['ID'].values
            for i in range(len(self.truth_df)):
                truth_allele = self.truth_df.loc[i, 'GLString'].split('^')[self.loci.index(locus)]
                impute_allele = self.impute_df.loc[i, 'SLUG_GLString'].split('^')[self.loci.index(locus)]
                truth_locus = truth_allele.split('*')[0]
                impute_locus = impute_allele.split('*')[0]
                # Only map DRB3/4/5/X to DRB345
                if truth_locus.startswith('DRB3') or truth_locus.startswith('DRB4') or truth_locus.startswith('DRB5') or truth_locus.startswith('DRBX'):
                    truth_allele = truth_allele.replace(truth_locus + '*', 'DRB345*')
                if impute_locus.startswith('DRB3') or impute_locus.startswith('DRB4') or impute_locus.startswith('DRB5') or impute_locus.startswith('DRBX'):
                    impute_allele = impute_allele.replace(impute_locus + '*', 'DRB345*')
                truth_alleles = sorted(truth_allele.split('+'))
                impute_alleles = sorted(impute_allele.split('+'))
                y_true.append(1 if truth_alleles == impute_alleles else 0)
            y_true = np.array(y_true)
            y_prob = self.impute_df[prob_col].values
            y_pred = (y_prob > 0.5).astype(int)
            df = pd.DataFrame({'ID': ids, 'y_true': y_true, 'y_pred': y_pred, 'y_prob': y_prob})
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
            truth_allele = self.truth_df.loc[i, 'GLString'].split('^')[self.loci.index(locus)]
            impute_allele = self.impute_df.loc[i, 'SLUG_GLString'].split('^')[self.loci.index(locus)]
            truth_locus = truth_allele.split('*')[0]
            impute_locus = impute_allele.split('*')[0]
            # Only map DRB3/4/5/X to DRB345
            if truth_locus.startswith('DRB3') or truth_locus.startswith('DRB4') or truth_locus.startswith('DRB5') or truth_locus.startswith('DRBX'):
                truth_allele = truth_allele.replace(truth_locus + '*', 'DRB345*')
            if impute_locus.startswith('DRB3') or impute_locus.startswith('DRB4') or impute_locus.startswith('DRB5') or impute_locus.startswith('DRBX'):
                impute_allele = impute_allele.replace(impute_locus + '*', 'DRB345*')
            truth_alleles = sorted(truth_allele.split('+'))
            impute_alleles = sorted(impute_allele.split('+'))
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
        # Use BaseAnalysis._align_data to align truth and impute DataFrames
        base = BaseAnalysis.__new__(BaseAnalysis)
        base.truth_data = self.truth_df
        base.impute_data = self.impute_df
        base._align_data()
        self.truth_df = base.truth_data
        self.impute_df = base.impute_data
        self.analysis_types = self._detect_analysis_types()

    def _detect_analysis_types(self):
        # Only use multilocus columns, ignore SLUG_GLString and GENO_*_Prob
        types = ['multiloc']
        for col in self.impute_df.columns:
            if col.endswith('_GLString') and not col.startswith('SLUG'):
                t = col.replace('_GLString', '')
                if t != 'multiloc':
                    types.append(t)
        return types

    def _normalize_glstring(self, glstring):
        """Map DRB3/4/5/X loci to DRB345 in a GLString."""
        alleles = glstring.split('^')
        norm_alleles = []
        for allele in alleles:
            locus = allele.split('*')[0]
            if locus.startswith('DRB3') or locus.startswith('DRB4') or locus.startswith('DRB5') or locus.startswith('DRBX'):
                norm_alleles.append(allele.replace(locus + '*', 'DRB345*'))
            else:
                norm_alleles.append(allele)
        return '^'.join(norm_alleles)

    def _get_loci_for_analysis_type(self, analysis_type):
        # Map analysis type to loci
        mapping = {
            'multiloc': ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1'],
            'sevenloc': ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1'],
            'classI': ['A', 'C', 'B'],
            'DRDQ': ['DRB345', 'DRB1', 'DQA1', 'DQB1'],
            'DR': ['DRB345', 'DRB1'],
            'DQ': ['DQA1', 'DQB1'],
        }
        return mapping.get(analysis_type, None)

    def _extract_loci_from_glstring(self, glstring):
        """Return a dict mapping locus to allele string, mapping DRB3/4/5/X to DRB345 and grouping all DRB3/4/5/X alleles together. Sort alleles within each locus."""
        alleles = glstring.split('^')
        locus_map = {}
        drb345_alleles = []
        for allele in alleles:
            locus = allele.split('*')[0]
            if locus.startswith('DRB3') or locus.startswith('DRB4') or locus.startswith('DRB5') or locus.startswith('DRBX'):
                drb345_alleles.append(allele.replace(locus + '*', 'DRB345*'))
            else:
                locus_map[locus] = allele
        if drb345_alleles:
            # Sort alleles within DRB345
            drb345_sorted = '+'.join(sorted(drb345_alleles))
            locus_map['DRB345'] = drb345_sorted
        # Sort alleles within each locus
        for locus in locus_map:
            alleles = locus_map[locus].split('+')
            locus_map[locus] = '+'.join(sorted(alleles))
        return locus_map

    def run(self):
        results = {}
        for analysis_type in self.analysis_types:
            if analysis_type == 'multiloc':
                gl_col = 'GLString'
                prob_col = 'HapPair_Prob'
            else:
                gl_col = f'{analysis_type}_GLString'
                prob_col = f'{analysis_type}_Prob' if f'{analysis_type}_Prob' in self.impute_df.columns else 'HapPair_Prob'
            if gl_col not in self.impute_df.columns or prob_col not in self.impute_df.columns:
                continue
            loci_for_type = self._get_loci_for_analysis_type(analysis_type)
            y_true = []
            for i in range(len(self.truth_df)):
                truth_locus_map = self._extract_loci_from_glstring(self.truth_df.loc[i, 'GLString'])
                if loci_for_type:
                    truth_gl = '^'.join([truth_locus_map[locus] for locus in loci_for_type if locus in truth_locus_map])
                else:
                    truth_gl = self._normalize_glstring(self.truth_df.loc[i, 'GLString'])
                impute_gl = self._normalize_glstring(self.impute_df.loc[i, gl_col])
                y_true.append(1 if truth_gl == impute_gl else 0)
            y_true = np.array(y_true)
            y_prob = self.impute_df[prob_col].values
            y_pred = (y_prob > 0.5).astype(int)
            results[analysis_type] = self._metrics(y_true, y_pred, y_prob)
        return results

    def get_results_df(self):
        """Return a dict of DataFrames with columns ID, y_true, y_pred, y_prob for each analysis type."""
        results_df = {}
        for analysis_type in self.analysis_types:
            if analysis_type == 'multiloc':
                gl_col = 'GLString'
                prob_col = 'HapPair_Prob'
            else:
                gl_col = f'{analysis_type}_GLString'
                prob_col = f'{analysis_type}_Prob' if f'{analysis_type}_Prob' in self.impute_df.columns else 'HapPair_Prob'
            if gl_col not in self.impute_df.columns or prob_col not in self.impute_df.columns:
                continue
            loci_for_type = self._get_loci_for_analysis_type(analysis_type)
            ids = self.impute_df['ID'].values
            y_true = []
            for i in range(len(self.truth_df)):
                truth_locus_map = self._extract_loci_from_glstring(self.truth_df.loc[i, 'GLString'])
                if loci_for_type:
                    truth_gl = '^'.join([truth_locus_map[locus] for locus in loci_for_type if locus in truth_locus_map])
                else:
                    truth_gl = self._normalize_glstring(self.truth_df.loc[i, 'GLString'])
                impute_gl = self._normalize_glstring(self.impute_df.loc[i, gl_col])
                y_true.append(1 if truth_gl == impute_gl else 0)
            y_true = np.array(y_true)
            y_prob = self.impute_df[prob_col].values
            y_pred = (y_prob > 0.5).astype(int)
            df = pd.DataFrame({'ID': ids, 'y_true': y_true, 'y_pred': y_pred, 'y_prob': y_prob})
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
        loci_for_type = self._get_loci_for_analysis_type(analysis_type)
        y_true = []
        imputed_loci = []
        true_loci = []
        for i in range(len(self.truth_df)):
            truth_locus_map = self._extract_loci_from_glstring(self.truth_df.loc[i, 'GLString'])
            if loci_for_type:
                truth_gl = '^'.join([truth_locus_map[locus] for locus in loci_for_type if locus in truth_locus_map])
            else:
                truth_gl = self._normalize_glstring(self.truth_df.loc[i, 'GLString'])
            impute_gl = self._normalize_glstring(self.impute_df.loc[i, gl_col])
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

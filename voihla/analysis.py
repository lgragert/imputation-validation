"""Analysis classes with preprocessing support."""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from .core import BaseAnalysis
from .preprocessing import ImputationPreprocessor
from .utils import GLStringParser


class SingleLocusAnalysis(BaseAnalysis):
    """Single locus analysis with preprocessing support."""

    def __init__(self, truth_file: str, impute_file: str,
                 raw_imputation: bool = False, high_res: bool = True):
        """Initialize with optional raw imputation preprocessing."""
        if raw_imputation:
            # Preprocess raw imputation file
            preprocessor = ImputationPreprocessor()
            processed_data = preprocessor.process_raw_imputation(impute_file)
            # Use processed data instead of file
            self.truth_data = pd.read_csv(truth_file, dtype={"ID": str})
            self.impute_data = processed_data
            self._align_data()
        else:
            super().__init__(truth_file, impute_file)

        self.high_res = high_res
        self.loci = self._get_available_loci()
        self.results = {}

    def _get_available_loci(self) -> List[str]:
        """Get list of available loci from the data."""
        all_loci = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
        available_loci = []

        # Check which loci are available in the data
        for locus in all_loci:
            prob_col = f'GENO_{locus}_Prob'
            if prob_col in self.impute_data.columns:
                available_loci.append(locus)

        # Filter by resolution if needed
        if not self.high_res:
            # Keep only basic loci for low resolution
            basic_loci = ['A', 'C', 'B', 'DRB1', 'DQB1']
            available_loci = [locus for locus in available_loci if locus in basic_loci]

        return available_loci

    def run_analysis(self):
        """Run single locus analysis for all available loci."""
        if not self.loci:
            raise ValueError("No loci available for analysis")

        # Parse GLStrings
        parser = GLStringParser()
        truth_parsed = parser.parse_glstring(self.truth_data, 'GLString')
        impute_parsed = parser.parse_glstring(self.impute_data, 'SLUG_GLString')

        for locus in self.loci:
            # Check if locus data exists in parsed data
            if (f'{locus}_1' in truth_parsed.columns and
                f'{locus}_1' in impute_parsed.columns):
                self.results[locus] = self._analyze_locus(
                    truth_parsed, impute_parsed, locus
                )

        return self.results

    def _analyze_locus(self, truth_data: pd.DataFrame,
                       impute_data: pd.DataFrame, locus: str) -> dict:
        """Analyze a single locus."""
        y_true = []
        y_prob = self.impute_data[f'GENO_{locus}_Prob'].values

        for i in range(len(truth_data)):
            truth_alleles = sorted([
                truth_data.loc[i, f'{locus}_1'],
                truth_data.loc[i, f'{locus}_2']
            ])
            impute_alleles = sorted([
                impute_data.loc[i, f'{locus}_1'],
                impute_data.loc[i, f'{locus}_2']
            ])

            # Perfect match = 1, otherwise 0
            y_true.append(1 if truth_alleles == impute_alleles else 0)

        y_true = np.array(y_true)
        y_pred = (y_prob > 0.5).astype(int)

        return self.calculate_metrics(y_true, y_pred, y_prob)


class MultiLocusAnalysis(BaseAnalysis):
    """Multi-locus unphased genotype analysis."""

    def __init__(self, truth_file: str, impute_file: str, raw_imputation: bool = False):
        """Initialize multi-locus analysis."""
        if raw_imputation:
            preprocessor = ImputationPreprocessor()
            processed_data = preprocessor.process_raw_imputation(impute_file)
            self.truth_data = pd.read_csv(truth_file, dtype={"ID": str})
            self.impute_data = processed_data
            self._align_data()
        else:
            super().__init__(truth_file, impute_file)

        self.analysis_types = self._get_available_analysis_types()
        self.results = {}

    def _get_available_analysis_types(self) -> List[str]:
        """Get available analysis types based on data columns."""
        available_types = []

        # Check what columns exist in impute_data
        columns = set(self.impute_data.columns)

        # Define required columns for each analysis type
        analysis_requirements = {
            '9loc': ['9loc_GLString', 'HapPair_Prob'],
            '7loc': ['7loc_GLString', '7loc_Prob'],
            'ClassI': ['ClassI_GLString', 'ClassI_Prob'],
            'DRDQ': ['DRDQ_GLString', 'DRDQ_Prob'],
            'DR': ['DR_GLString', 'DR_Prob'],
            'DQ': ['DQ_GLString', 'DQ_Prob']
        }

        # Check which analysis types have required columns
        for analysis_type, required_cols in analysis_requirements.items():
            if all(col in columns for col in required_cols):
                available_types.append(analysis_type)

        # Ensure at least 2 loci for multi-locus analysis
        if not available_types:
            # Fallback: check if we can create basic multi-locus from available loci
            available_types = self._create_fallback_analysis_types()

        return available_types

    def _create_fallback_analysis_types(self) -> List[str]:
        """Create fallback analysis types from available single loci."""
        available_loci = []
        all_loci = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

        for locus in all_loci:
            if f'GENO_{locus}_Prob' in self.impute_data.columns:
                available_loci.append(locus)

        fallback_types = []

        # Create analysis types based on available loci
        if len(available_loci) >= 9:
            fallback_types.append('9loc')
        if len(available_loci) >= 7:
            fallback_types.append('7loc')
        if set(['A', 'C', 'B']).issubset(set(available_loci)):
            fallback_types.append('ClassI')
        if any(locus in available_loci for locus in ['DRB345', 'DRB1', 'DQA1', 'DQB1']):
            if len([l for l in ['DRB345', 'DRB1', 'DQA1', 'DQB1'] if l in available_loci]) >= 2:
                fallback_types.append('DRDQ')

        return fallback_types

    def run_analysis(self):
        """Run multi-locus analysis for all available types."""
        if not self.analysis_types:
            raise ValueError("No multi-locus analysis types available (need at least 2 loci)")

        for analysis_type in self.analysis_types:
            try:
                self.results[analysis_type] = self._analyze_multilocus(analysis_type)
            except KeyError as e:
                print(f"Warning: Skipping {analysis_type} analysis due to missing data: {e}")
                continue

        return self.results

    def _analyze_multilocus(self, analysis_type: str) -> dict:
        """Analyze multi-locus genotypes."""
        # Get appropriate GLString column and probability
        if analysis_type == '9loc':
            if '9loc_GLString' in self.impute_data.columns:
                glstring_col = '9loc_GLString'
                prob_col = 'HapPair_Prob'
            else:
                glstring_col = 'SLUG_GLString'
                prob_col = 'HapPair_Prob'
        else:
            glstring_col = f'{analysis_type}_GLString'
            prob_col = f'{analysis_type}_Prob'

        # Check if required columns exist
        if glstring_col not in self.impute_data.columns:
            raise KeyError(f"Missing column: {glstring_col}")
        if prob_col not in self.impute_data.columns:
            raise KeyError(f"Missing column: {prob_col}")

        # Prepare truth data based on analysis type
        truth_glstrings = self._prepare_truth_glstrings(analysis_type)

        y_true = []
        y_prob = self.impute_data[prob_col].values

        for i in range(len(self.truth_data)):
            truth_glstring = truth_glstrings.iloc[i]
            impute_glstring = self.impute_data.loc[i, glstring_col]

            y_true.append(1 if truth_glstring == impute_glstring else 0)

        y_true = np.array(y_true)
        y_pred = (y_prob > 0.5).astype(int)

        return self.calculate_metrics(y_true, y_pred, y_prob)

    def _prepare_truth_glstrings(self, analysis_type: str) -> pd.Series:
        """Prepare truth GLStrings for different analysis types."""
        parser = GLStringParser()

        if analysis_type == '9loc':
            return self.truth_data['GLString']
        elif analysis_type == '7loc':
            return parser.create_subset_glstring(
                self.truth_data, ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1']
            )
        elif analysis_type == 'ClassI':
            return parser.create_subset_glstring(
                self.truth_data, ['A', 'C', 'B']
            )
        elif analysis_type == 'DRDQ':
            return parser.create_subset_glstring(
                self.truth_data, ['DRB345', 'DRB1', 'DQA1', 'DQB1']
            )
        elif analysis_type == 'DR':
            return parser.create_subset_glstring(
                self.truth_data, ['DRB345', 'DRB1']
            )
        elif analysis_type == 'DQ':
            return parser.create_subset_glstring(
                self.truth_data, ['DQA1', 'DQB1']
            )
        else:
            return self.truth_data['GLString']
"""Eplet analysis using Monte Carlo sampling."""

import pandas as pd
import numpy as np
import requests
import time
from typing import Dict, List, Tuple
from .core import BaseAnalysis
from .plotting import CalibrationPlotter


class MonteCarloEpletAnalysis:
    """Monte Carlo eplet mismatch analysis."""

    def __init__(self, api_key: str, delay: float = 1.0):
        """Initialize with API key and request delay."""
        self.api_key = api_key
        self.delay = delay
        self.base_url = "https://api.epregistry.com.br/eplet_mismatches"

    def create_random_pairs(self, truth_file: str, n_pairs: int,
                           random_state: int = 42) -> pd.DataFrame:
        """Create random donor-recipient pairs."""
        truth_data = pd.read_csv(truth_file)

        np.random.seed(random_state)
        pairs = []
        while len(pairs) < n_pairs:
            don_idx = np.random.choice(len(truth_data))
            rec_idx = np.random.choice(len(truth_data))
            # Allow self-pairing if needed to guarantee n_pairs
            pairs.append({
                'DON_ID': truth_data.iloc[don_idx]['ID'],
                'DON_GLString': truth_data.iloc[don_idx]['GLString'],
                'REC_ID': truth_data.iloc[rec_idx]['ID'],
                'REC_GLString': truth_data.iloc[rec_idx]['GLString']
            })

        return pd.DataFrame(pairs)

    def query_eplet_api(self, donor_alleles: str,
                       recipient_alleles: str) -> Dict:
        """Query the eplet registry API."""
        url = f"{self.base_url}?from={self.api_key}"
        url += f"&immunizer_alleles={donor_alleles}"
        url += f"&patient_alleles={recipient_alleles}"

        try:
            response = requests.get(url, headers={'Accept': 'application/json'})
            response.raise_for_status()
            time.sleep(self.delay)
            return response.json()
        except requests.RequestException as e:
            print(f"API request failed: {e}")
            return {}

    def analyze_eplet_mismatches(self, pairs_df: pd.DataFrame,
                               analysis_type: str = 'DRDQ') -> pd.DataFrame:
        """Analyze eplet mismatches for donor-recipient pairs."""
        results = []

        for idx, row in pairs_df.iterrows():
            if idx % 100 == 0:  # API rate limiting
                time.sleep(60.0)

            donor_alleles = row['DON_GLString']
            recipient_alleles = row['REC_GLString']

            api_result = self.query_eplet_api(donor_alleles, recipient_alleles)

            # Extract relevant eplet information
            eplet_data = self._extract_eplet_data(api_result, analysis_type)

            result_row = {
                'DON_ID': row['DON_ID'],
                'REC_ID': row['REC_ID'],
                'DON_GLString': donor_alleles,
                'REC_GLString': recipient_alleles,
                **eplet_data
            }

            results.append(result_row)
            time.sleep(self.delay)

        return pd.DataFrame(results)

    def _extract_eplet_data(self, api_result: Dict,
                          analysis_type: str) -> Dict:
        """Extract eplet data from API response."""
        if analysis_type == 'DRDQ':
            key = 'ALL'
        elif analysis_type == 'DR':
            key = 'DRB'
        elif analysis_type == 'DQ':
            key = 'DQ'
        else:
            key = 'ALL'

        eplet_info = api_result.get(key, {})

        return {
            'eplet_count': eplet_info.get('quantity', 0),
            'eplet_details': str(eplet_info.get('details', []))
        }


class EpletAnalysis:
    """Analysis of eplet-level imputation performance."""

    def __init__(self, truth_file: str, impute_file: str):
        """Initialize eplet analysis."""
        self.truth_data = pd.read_csv(truth_file)
        self.impute_data = pd.read_csv(impute_file)

    def run_calibration_analysis(self, n_bins: int = 10) -> Dict:
        """Run calibration analysis on eplet predictions."""
        plotter = CalibrationPlotter(n_bins)
        results = {}

        # Analyze eplet count predictions
        if 'eplet_count' in self.impute_data.columns:
            truth_counts = self.truth_data['eplet_count'].values
            impute_counts = self.impute_data['eplet_count'].values

            # Convert to binary classification (0 mismatches vs >0)
            y_true = (truth_counts > 0).astype(int)
            y_prob = self.impute_data.get('pair_probability',
                                         np.ones(len(y_true)) * 0.5)

            results['calibration_plot'] = plotter.create_calibration_plot(
                y_true, y_prob, 'Eplet Mismatch Calibration'
            )

        return results
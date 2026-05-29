"""Eplet analysis using Monte Carlo sampling and the EpRegistry API."""

import pandas as pd
import requests
import time
from typing import Dict, Optional
from .plotting import CalibrationPlotter


# Maps which_impute label to (API key prefix, quantity col, details col)
_LOCUS_CONFIG = {
    'Total': ('ALL', 'ALL_quantity', 'ALL_details'),
    'DR':   ('DRB', 'DRB_quantity', 'DRB_details'),
    'DQ':   ('DQ',  'DQ_quantity',  'DQ_details'),
    'ABC':  ('ABC', 'ABC_quantity', 'ABC_details'),
    'DP':   ('DP',  'DP_quantity',  'DP_details'),
    'MICA': ('MICA','MICA_quantity','MICA_details'),
}

_SKIP_API_KEYS = {'version'}


class MonteCarloEpletAnalysis:
    """Query the EpRegistry API for truth and imputed eplet mismatches.

    Mirrors the workflow in montecarlo_pairings.py as a reusable class:

    1. Sample random donor-recipient pairs from a truth file.
    2. Query the API with high-resolution truth genotypes -> truth eplet table.
    3. Filter matching rows from an imputation file.
    4. Query the API with low-resolution imputed genotypes -> imputed eplet table.

    Usage::

        mc = MonteCarloEpletAnalysis.from_key_file('api.key')

        truth_pairs = mc.sample_pairs('DRDQ_pairs_truth.csv', n_pairs=100,
                                      which_impute='DRDQ',
                                      save_path='DRDQ_pairs_truth_100.csv')

        truth_eplets = mc.build_truth_eplet_table(truth_pairs, which_impute='DRDQ',
                                                  save_path='DRDQ_eplet_truth_table100.csv')

        impute_rows = mc.get_imputation_rows('DRDQ_pairs_imputation.csv',
                                             truth_pairs, which_impute='DRDQ',
                                             save_path='DRDQ_pairs_imputation_100.csv')

        impute_eplets = mc.build_impute_eplet_table(impute_rows, which_impute='DRDQ',
                                                    save_path='DRDQ_eplet_lowres_impute100.csv')
    """

    def __init__(self, api_key: str, delay: float = 1.0):
        """
        Parameters
        ----------
        api_key : str
            EpRegistry API key.
        delay : float
            Seconds to sleep between API calls (default 1.0).
        """
        self.api_key = api_key
        self.delay = delay
        self._base_url = 'https://api.epregistry.com.br/eplet_mismatches'

    @classmethod
    def from_key_file(cls, key_file: str, delay: float = 1.0) -> 'MonteCarloEpletAnalysis':
        """Instantiate by reading the API key from a plain-text file."""
        with open(key_file) as f:
            api_key = f.readline().strip()
        return cls(api_key, delay)

    # ------------------------------------------------------------------
    # Step 1: sample pairs
    # ------------------------------------------------------------------

    def sample_pairs(self, truth_file: str, n_pairs: int,
                     random_state: int = 2272024,
                     save_path: Optional[str] = None) -> pd.DataFrame:
        """Sample n_pairs random donor-recipient rows from truth_file.

        Parameters
        ----------
        truth_file : str
            CSV with at minimum: DON_ID, REC_ID, DON_GLString, REC_GLString.
        n_pairs : int
            Number of pairs to sample.
        random_state : int
            Random seed for reproducibility (default matches montecarlo_pairings.py).
        save_path : str, optional
            If provided, saves the sampled DataFrame as a CSV.

        Returns
        -------
        pd.DataFrame
        """
        truth_data = pd.read_csv(truth_file, header=0)
        sampled = truth_data.sample(n=n_pairs, random_state=random_state).reset_index(drop=True)
        if save_path:
            sampled.to_csv(save_path, index=False, header=True)
        return sampled

    # ------------------------------------------------------------------
    # Step 2: truth eplet table
    # ------------------------------------------------------------------

    def build_truth_eplet_table(self, truth_pairs_df: pd.DataFrame,
                                which_impute: str,
                                save_path: Optional[str] = None) -> pd.DataFrame:
        """Query the EpRegistry API with high-res truth genotypes.

        Parameters
        ----------
        truth_pairs_df : pd.DataFrame
            Must have columns: DON_ID, REC_ID, DON_GLString, REC_GLString.
        which_impute : str
            One of 'DRDQ', 'DR', 'DQ'. Determines which eplet columns are kept.
        save_path : str, optional
            If provided, saves the result as a CSV.

        Returns
        -------
        pd.DataFrame
            Columns: ID, {eplet}_quantity, {eplet}_details
        """
        if which_impute not in _LOCUS_CONFIG:
            raise ValueError(f"which_impute must be one of {list(_LOCUS_CONFIG)}")
        _, qty_col, det_col = _LOCUS_CONFIG[which_impute]

        rows = []
        truth_pairs_df = truth_pairs_df.reset_index(drop=True)
        for line in range(len(truth_pairs_df)):
            donor_id   = str(truth_pairs_df.loc[line, 'DON_ID'])
            recip_id   = str(truth_pairs_df.loc[line, 'REC_ID'])
            donor_geno = str(truth_pairs_df.loc[line, 'DON_GLString'])
            recip_geno = str(truth_pairs_df.loc[line, 'REC_GLString'])
            id_pair    = f'{donor_id}+{recip_id}'

            if line % 100 == 0:
                time.sleep(60.0)

            parsed = self._parse_api_response(
                self._query_api(donor_geno, recip_geno), id_pair
            )
            rows.append(parsed)
            time.sleep(self.delay)

        result = pd.DataFrame(rows)
        if save_path:
            result.to_csv(save_path, index=False, header=True)
        return result

    # ------------------------------------------------------------------
    # Step 3: filter imputation rows for sampled pairs
    # ------------------------------------------------------------------

    def get_imputation_rows(self, impute_file: str,
                            truth_pairs_df: pd.DataFrame,
                            which_impute: str,
                            save_path: Optional[str] = None) -> pd.DataFrame:
        """Return all imputation rows corresponding to pairs in truth_pairs_df.

        Parameters
        ----------
        impute_file : str
            CSV with columns: DON_ID, REC_ID, PairProb_{which_impute},
            DON_{which_impute}, REC_{which_impute}.
        truth_pairs_df : pd.DataFrame
            Sampled truth pairs (output of sample_pairs).
        which_impute : str
            One of 'DRDQ', 'DR', 'DQ'.
        save_path : str, optional
            If provided, saves the result as a CSV.

        Returns
        -------
        pd.DataFrame
            All imputation rows for the sampled pairs, sorted by PairProb descending.
        """
        impute_data = pd.read_csv(impute_file, header=0)
        prob_col = f'PairProb_{which_impute}'

        frames = []
        for _, row in truth_pairs_df.iterrows():
            don, rec = row['DON_ID'], row['REC_ID']
            match = impute_data[
                (impute_data['DON_ID'] == don) & (impute_data['REC_ID'] == rec)
            ]
            if match.empty:
                match = impute_data[
                    (impute_data['DON_ID'] == rec) & (impute_data['REC_ID'] == don)
                ]
            match = match.sort_values(by=prob_col, ascending=False)
            frames.append(match)

        result = pd.concat(frames).reset_index(drop=True)
        if save_path:
            result.to_csv(save_path, index=False, header=True)
        return result

    # ------------------------------------------------------------------
    # Step 4: imputed eplet table
    # ------------------------------------------------------------------

    def build_impute_eplet_table(self, impute_pairs_df: pd.DataFrame,
                                 which_impute: str,
                                 save_path: Optional[str] = None) -> pd.DataFrame:
        """Query the EpRegistry API with low-res imputed genotypes.

        Parameters
        ----------
        impute_pairs_df : pd.DataFrame
            Output of get_imputation_rows. Must have: DON_ID, REC_ID,
            PairProb_{which_impute}, DON_{which_impute}, REC_{which_impute}.
        which_impute : str
            One of 'DRDQ', 'DR', 'DQ'.
        save_path : str, optional
            If provided, saves the result as a CSV.

        Returns
        -------
        pd.DataFrame
            Columns: DON_ID, REC_ID, PairProb_{which_impute},
            DON_{which_impute}, REC_{which_impute}, {eplet}_quantity, {eplet}_details
        """
        if which_impute not in _LOCUS_CONFIG:
            raise ValueError(f"which_impute must be one of {list(_LOCUS_CONFIG)}")
        _, qty_col, det_col = _LOCUS_CONFIG[which_impute]
        prob_col = f'PairProb_{which_impute}'
        don_col  = f'DON_{which_impute}'
        rec_col  = f'REC_{which_impute}'

        rows = []
        impute_pairs_df = impute_pairs_df.reset_index(drop=True)
        for line in range(len(impute_pairs_df)):
            donor_id   = str(impute_pairs_df.loc[line, 'DON_ID'])
            recip_id   = str(impute_pairs_df.loc[line, 'REC_ID'])
            donor_geno = str(impute_pairs_df.loc[line, don_col])
            recip_geno = str(impute_pairs_df.loc[line, rec_col])
            pair_prob  = impute_pairs_df.loc[line, prob_col]
            id_pair    = f'{donor_id}+{recip_id}'

            if line % 100 == 0:
                time.sleep(60.0)

            parsed = self._parse_api_response(
                self._query_api(donor_geno, recip_geno), id_pair
            )
            parsed['DON_ID'] = donor_id
            parsed['REC_ID'] = recip_id
            parsed[don_col]  = donor_geno
            parsed[rec_col]  = recip_geno
            parsed[prob_col] = pair_prob
            rows.append(parsed)
            time.sleep(self.delay)

        eplet_cols = [c for c in pd.DataFrame(rows).columns
                      if c not in {'ID', 'DON_ID', 'REC_ID', don_col, rec_col, prob_col}]
        result = pd.DataFrame(rows)[['DON_ID', 'REC_ID', prob_col, don_col, rec_col] + eplet_cols]
        if save_path:
            result.to_csv(save_path, index=False, header=True)
        return result

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _query_api(self, immunizer_alleles: str, patient_alleles: str) -> Dict:
        """Hit the EpRegistry API and return the raw JSON dict."""
        url = (f'{self._base_url}?from={self.api_key}'
               f'&immunizer_alleles={immunizer_alleles}'
               f'&patient_alleles={patient_alleles}')
        try:
            r = requests.get(url, headers={'Accept': 'application/json'})
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            print(f'API request failed for {immunizer_alleles} vs {patient_alleles}: {e}')
            return {}

    def _parse_api_response(self, api_result: Dict, id_pair: str) -> Dict:
        """Parse and clean one API response into a flat dict.

        Skips ABC, DP, MICA, and version keys (same as montecarlo_pairings.py).
        Converts detail lists from '[eplet1, eplet2]' to 'eplet1_eplet2'.
        """
        row: Dict = {'ID': id_pair}
        for key, value in api_result.items():
            if key in _SKIP_API_KEYS:
                continue
            for label, extended in value.items():
                cell = str(extended)
                if label == 'details':
                    cell = (cell
                            .replace(', ', '_')
                            .replace("'", '')
                            .replace('[', '')
                            .replace(']', ''))
                row[f'{key}_{label}'] = cell
        return row


class EpletAnalysis:
    """Calibration analysis comparing predicted vs true eplet counts.

    Accepts two file formats for the impute file:

    - **With PairProb** (output of MonteCarloEpletAnalysis.build_impute_eplet_table):
      DON_ID, REC_ID, PairProb_{locus}, DON_{locus}, REC_{locus},
      {eplet}_quantity, {eplet}_details.
      Probabilities are aggregated by summing PairProb per count value.

    - **Without PairProb** (e.g. EpRegistry_lowres_impute.csv with pair ID as index):
      probabilities are derived as empirical frequency across MC samples.

    truth_file must have columns: ID, {eplet}_quantity.

    Usage::

        ea = EpletAnalysis('DRDQ_eplet_truth_table100.csv',
                           'DRDQ_eplet_lowres_impute100.csv')
        plots = ea.run_calibration_analysis(n_bins=10)
        plots['calibration_plot_DR'].savefig('Calibration_DR.png', bbox_inches='tight')
    """

    def __init__(self, truth_file: str, impute_file: str):
        self.truth_data = pd.read_csv(truth_file, header=0)
        raw = pd.read_csv(impute_file, header=0)

        # Normalise to a unified index: DON_ID+REC_ID
        if 'DON_ID' in raw.columns and 'REC_ID' in raw.columns:
            raw.index = raw['DON_ID'] + '+' + raw['REC_ID']
            raw = raw.drop(columns=['DON_ID', 'REC_ID'])
        else:
            # EpRegistry_lowres_impute format: unnamed first column is the pair ID
            raw = raw.set_index(raw.columns[0])

        self.impute_data = raw
        self._prob_cols = {col.replace('PairProb_', ''): col
                          for col in raw.columns if col.startswith('PairProb_')}

        if not self._prob_cols:
            raise ValueError(
                f"No PairProb_* column found in '{impute_file}'. "
                "The impute file must contain at least one probability column "
                "(e.g. PairProb_DRDQ, PairProb_DR, PairProb_DQ). "
                "Use the output of MonteCarloEpletAnalysis.build_impute_eplet_table(), "
                "not the raw EpRegistry output."
            )

        # Populated by run_calibration_analysis()
        # pair_results:  {locus_label: DataFrame(ID, y_true, y_pred, y_prob)}
        # summary_results: {locus_label: {'brier': float, 'roc_auc': float,
        #                                  'confusion_matrix': DataFrame}}
        self.pair_results: Dict[str, pd.DataFrame] = {}
        self.summary_results: Dict[str, Dict] = {}

    def _aggregate_by_pair(self, locus_col: str,
                           which_impute: Optional[str] = None) -> pd.DataFrame:
        """Per pair: most probable predicted count and its probability."""
        prob_col = self._prob_cols.get(which_impute) if which_impute else None
        records = []
        for pair_id, group in self.impute_data.groupby(level=0):
            if prob_col and prob_col in group.columns:
                count_probs = group.groupby(locus_col)[prob_col].sum()
                top_count = count_probs.idxmax()
                top_prob  = count_probs.max()
            else:
                dist = group[locus_col].value_counts(normalize=True)
                top_count = dist.idxmax()
                top_prob  = dist.max()
            records.append({
                'ID': pair_id,
                f'pred_{locus_col}': top_count,
                f'prob_{locus_col}': top_prob,
            })
        return pd.DataFrame(records)

    def run_calibration_analysis(self, n_bins: int = 10) -> Dict:
        """Run calibration analysis for all loci present in both files.

        Results are also stored on the instance after calling this method:

        - ``self.pair_results[locus]``: DataFrame with columns ID, y_true, y_pred, y_prob
        - ``self.summary_results[locus]``: dict with keys 'brier', 'roc_auc', 'confusion_matrix'

        Returns
        -------
        dict
            Keys like 'calibration_plot_DR', values are matplotlib Figures.
        """
        plotter = CalibrationPlotter(n_bins)
        plots = {}
        self.pair_results = {}
        self.summary_results = {}

        locus_map = {
            'ALL_quantity': ('Total', 'DRDQ'),
            'DRB_quantity': ('DR',    'DR'),
            'DQ_quantity':  ('DQ',    'DQ'),
            'ABC_quantity': ('ABC',   'ABC'),
            'DP_quantity':  ('DP',    'DP'),
            'MICA_quantity':('MICA',  'MICA'),
        }

        truth_sorted = self.truth_data.sort_values('ID').reset_index(drop=True)

        for col, (locus_label, which_impute) in locus_map.items():
            if col not in self.impute_data.columns or col not in truth_sorted.columns:
                continue

            pred_df = self._aggregate_by_pair(col, which_impute)
            merged  = pd.merge(truth_sorted[['ID', col]], pred_df, on='ID', how='inner')

            y_true = (merged[col].astype(int) == merged[f'pred_{col}'].astype(int)).astype(int)
            y_prob = merged[f'prob_{col}'].to_numpy(dtype=float)
            y_pred = (y_prob >= 0.9).astype(int)

            analysis_results = pd.DataFrame({
                'ID':     merged['ID'],
                'y_true': y_true,
                'y_pred': y_pred,
                'y_prob': y_prob,
            })

            fig, confusion_mat, brier, roc_auc = plotter.calibration_plot(
                analysis_results, locus=locus_label,
                title=f'{locus_label} Eplet Count Calibration'
            )

            self.pair_results[locus_label] = analysis_results.reset_index(drop=True)
            self.summary_results[locus_label] = {
                'brier':            brier,
                'roc_auc':          roc_auc,
                'confusion_matrix': confusion_mat,
            }
            plots[f'calibration_plot_{locus_label}'] = fig

        return plots

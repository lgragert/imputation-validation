"""Utility functions for HLA data parsing and formatting."""

import pandas as pd
import numpy as np
from typing import List, Union


class GLStringParser:
    """Parse and manipulate GLString formatted HLA data."""

    def __init__(self):
        """Initialize parser."""
        self.loci_order = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

    def parse_glstring(self, data: pd.DataFrame,
                      glstring_col: str) -> pd.DataFrame:
        """Parse GLString into individual loci."""
        result = data.copy()

        # Split GLString by '^'
        loci_data = result[glstring_col].str.split('^', expand=True)
        loci_data.columns = self.loci_order[:loci_data.shape[1]]

        # Split each locus by '+'
        for locus in loci_data.columns:
            allele_data = loci_data[locus].str.split('+', expand=True)
            result[f'{locus}_1'] = allele_data[0]
            result[f'{locus}_2'] = allele_data[1] if allele_data.shape[1] > 1 else allele_data[0]

        return result

    def create_subset_glstring(self, data: pd.DataFrame,
                              loci_subset: List[str]) -> pd.Series:
        """Create GLString for subset of loci."""
        parsed = self.parse_glstring(data, 'GLString')

        glstrings = []
        for idx in range(len(parsed)):
            locus_strings = []
            for locus in loci_subset:
                allele1 = parsed.loc[idx, f'{locus}_1']
                allele2 = parsed.loc[idx, f'{locus}_2']
                locus_strings.append(f'{allele1}+{allele2}')
            glstrings.append('^'.join(locus_strings))

        return pd.Series(glstrings)

    def clean_eplet_string(self, eplet_series: pd.Series) -> pd.Series:
        """Clean eplet string formatting."""
        cleaned = eplet_series.str.replace(', ', '_')
        cleaned = cleaned.str.replace('\'', '')
        cleaned = cleaned.str.replace(r'[\[\]]', '', regex=True)
        return cleaned


class NGSFormatter:
    """Format NGS high-resolution data to GLString."""

    def __init__(self, ard_version: str = "3520"):
        """Initialize with ARD version."""
        try:
            import pyard
            self.ard = pyard.init(ard_version)
        except ImportError:
            self.ard = None
            warnings.warn("pyard not available, ARD reduction disabled")

    def format_ngs_data(self, ngs_data: pd.DataFrame) -> pd.DataFrame:
        """Format NGS data to GLString format."""
        result = ngs_data.copy()

        loci = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

        for locus in loci:
            result = self._add_locus_prefix(result, locus)

        # Create GLString
        result['GLString'] = (
            result['DON_A'] + '^' + result['DON_C'] + '^' + result['DON_B'] +
            '^' + result['DON_DRB345'] + '^' + result['DON_DRB1'] +
            '^' + result['DON_DQA1'] + '^' + result['DON_DQB1'] +
            '^' + result['DON_DPA1'] + '^' + result['DON_DPB1']
        )

        return result[['unos', 'GLString']].rename(columns={'unos': 'ID'})

    def _add_locus_prefix(self, data: pd.DataFrame, locus: str) -> pd.DataFrame:
        """Add locus prefix to typing data."""
        allele1_col = f'NGS_{locus}i'
        allele2_col = f'NGS_{locus}ii'

        if allele1_col not in data.columns or allele2_col not in data.columns:
            return data

        processed_alleles = []

        for idx in range(len(data)):
            allele1 = str(data.loc[idx, allele1_col])
            allele2 = str(data.loc[idx, allele2_col])

            # Handle special cases
            if locus == 'DRB345':
                allele1 = self._process_drb345(allele1)
                allele2 = self._process_drb345(allele2)
            else:
                allele1 = f'{locus}*{allele1}' if allele1 != 'nan' else f'{locus}*NNNN'
                allele2 = f'{locus}*{allele2}' if allele2 != 'nan' else f'{locus}*NNNN'

            # ARD reduction if available
            if self.ard and not allele1.endswith('NEW'):
                try:
                    allele1 = self.ard.redux(allele1, 'lgx')
                except:
                    pass

            if self.ard and not allele2.endswith('NEW'):
                try:
                    allele2 = self.ard.redux(allele2, 'lgx')
                except:
                    pass

            # Sort and combine
            sorted_alleles = sorted([allele1, allele2])
            processed_alleles.append('+'.join(sorted_alleles))

        data[f'DON_{locus}'] = processed_alleles
        return data

    def _process_drb345(self, allele: str) -> str:
        """Process DRB345 alleles with special handling."""
        if allele == 'nan':
            return 'DRBX*NNNN'
        elif allele == '5*01:08N':
            return 'DRB5*01:02'
        else:
            return f'DRB{allele}'
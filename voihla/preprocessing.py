"""Preprocessing functions to convert raw imputation data to analysis format."""

import pandas as pd
import numpy as np
from typing import List, Dict, Set
from .utils import GLStringParser


class ImputationPreprocessor:
    """Convert raw imputation data with flexible loci to analysis-ready format."""

    def __init__(self):
        """Initialize preprocessor with flexible loci detection."""
        self.standard_loci_order = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
        self.detected_loci = []

    def detect_loci_from_data(self, raw_file: str) -> List[str]:
        """Detect available loci from the raw data."""
        raw_data = pd.read_csv(raw_file)

        # Get a sample haplotype to detect loci
        sample_hap = raw_data.iloc[0]['Hap1']
        detected_loci = set()

        parts = sample_hap.split('~')
        for part in parts:
            if '*' in part:
                locus = part.split('*')[0]
                if locus in ['DRB3', 'DRB4', 'DRB5']:
                    detected_loci.add('DRB345')
                else:
                    detected_loci.add(locus)

        # Order loci according to standard order
        ordered_loci = [locus for locus in self.standard_loci_order if locus in detected_loci]

        return ordered_loci

    def process_raw_imputation(self, raw_file: str, output_file: str = None) -> pd.DataFrame:
        """Process raw imputation file with flexible loci to analysis format."""
        # Detect available loci first
        self.detected_loci = self.detect_loci_from_data(raw_file)

        raw_data = pd.read_csv(raw_file)
        top_predictions = raw_data.groupby('ID').first().reset_index()

        processed_data = self._convert_haplotypes_to_glstring(top_predictions)

        if output_file:
            processed_data.to_csv(output_file, index=False)

        return processed_data

    def _parse_haplotype(self, hap_string: str) -> Dict[str, str]:
        """Parse haplotype string with flexible loci."""
        alleles = {}
        parts = hap_string.split('~')

        for part in parts:
            if '*' in part:
                locus = part.split('*')[0]
                # Handle DRB345 loci (DRB3, DRB4, DRB5)
                if locus in ['DRB3', 'DRB4', 'DRB5']:
                    alleles['DRB345'] = part
                else:
                    alleles[locus] = part

        return alleles

    def _create_unphased_genotype(self, hap1: Dict, hap2: Dict) -> Dict[str, List[str]]:
        """Create unphased genotype from two haplotypes with flexible loci."""
        genotype = {}

        # Only process detected loci
        for locus in self.detected_loci:
            allele1 = hap1.get(locus, f'{locus}*NNNN')
            allele2 = hap2.get(locus, f'{locus}*NNNN')

            # Special handling for DRB345
            if locus == 'DRB345':
                if allele1.startswith('DRB345*NNNN') and not allele2.startswith('DRB345*NNNN'):
                    allele1 = 'DRB345*NNNN'
                elif allele2.startswith('DRB345*NNNN') and not allele1.startswith('DRB345*NNNN'):
                    allele2 = 'DRB345*NNNN'

            genotype[locus] = sorted([allele1, allele2])

        return genotype

    def _build_glstring(self, genotype: Dict[str, List[str]]) -> str:
        """Build GLString from genotype dictionary with flexible loci."""
        locus_strings = []

        # Build GLString only for detected loci
        for locus in self.detected_loci:
            if locus in genotype:
                locus_strings.append('+'.join(genotype[locus]))
            else:
                locus_strings.append(f'{locus}*NNNN+{locus}*NNNN')

        return '^'.join(locus_strings)

    def _create_all_subsets(self, genotype: Dict[str, List[str]],
                           probability: float) -> Dict[str, any]:
        """Create all possible subset GLStrings based on available loci."""
        subsets = {}

        # Define subset combinations - only include if loci are available
        subset_definitions = {
            'ClassI': ['A', 'C', 'B'],
            'DRDQ': ['DRB345', 'DRB1', 'DQA1', 'DQB1'],
            'DR': ['DRB345', 'DRB1'],
            'DQ': ['DQA1', 'DQB1']
        }

        # Add dynamic subsets based on available loci
        if len(self.detected_loci) >= 7:
            subset_definitions['7loc'] = [l for l in ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1']
                                         if l in self.detected_loci]

        if len(self.detected_loci) >= 9:
            subset_definitions['9loc'] = self.detected_loci

        for subset_name, loci in subset_definitions.items():
            # Check if we have enough loci for this subset
            available_subset_loci = [locus for locus in loci if locus in self.detected_loci]

            if len(available_subset_loci) >= 2 or subset_name in ['9loc', '7loc']:  # Multi-locus needs at least 2
                subset_strings = []
                for locus in available_subset_loci:
                    if locus in genotype:
                        subset_strings.append('+'.join(genotype[locus]))
                    else:
                        subset_strings.append(f'{locus}*NNNN+{locus}*NNNN')

                if subset_strings:  # Only add if we have data
                    subsets[f'{subset_name}_GLString'] = '^'.join(subset_strings)
                    subsets[f'{subset_name}_Prob'] = probability

        return subsets

    def _convert_haplotypes_to_glstring(self, data: pd.DataFrame) -> pd.DataFrame:
        """Convert haplotypes to GLString with flexible loci."""
        results = []

        for idx, row in data.iterrows():
            # Parse haplotype strings
            hap1_alleles = self._parse_haplotype(row['Hap1'])
            hap2_alleles = self._parse_haplotype(row['Hap2'])

            # Create unphased genotype
            genotype = self._create_unphased_genotype(hap1_alleles, hap2_alleles)

            # Build GLString
            glstring = self._build_glstring(genotype)

            result_row = {
                'ID': row['ID'],
                'SLUG_GLString': glstring,
                'HapPair_Prob': row['HapPair_Prob']
            }

            # Add full GLString if we have enough loci
            if len(self.detected_loci) >= 9:
                result_row['9loc_GLString'] = glstring

            # Add individual locus probabilities for detected loci
            for locus in self.detected_loci:
                result_row[f'GENO_{locus}_Prob'] = row['HapPair_Prob']

            # Add subset GLStrings and probabilities
            subsets = self._create_all_subsets(genotype, row['HapPair_Prob'])
            result_row.update(subsets)

            results.append(result_row)

        return pd.DataFrame(results)
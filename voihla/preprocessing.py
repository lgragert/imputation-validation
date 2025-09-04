"""Preprocessing functions to convert raw imputation data to analysis format."""

import pandas as pd
from collections import defaultdict
import gzip
import os


class ImputationPreprocessor:
    """Convert raw imputation data with flexible loci to analysis-ready format."""

    ALL_LOCI = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

    def __init__(self):
        self.gf_dicts = {locus: defaultdict(dict) for locus in self.ALL_LOCI}
        self.multiloc_freq = defaultdict(dict)
        self.other_freqs = defaultdict(lambda: defaultdict(dict))
        self.all_loci_found = []

    def open_impute_file(self, filepath):
        if filepath.endswith('.gz'):
            return gzip.open(filepath, "rt")
        else:
            return open(filepath, "r")

    def process_files(self, filename_list):
        """Process multiple imputation files to extract GLString and probabilities."""
        for filename in filename_list:
            impute_outfile = self.open_impute_file(filename)
            # Detect loci from first data line (order preserved)
            for line in impute_outfile:
                if line.strip() == "" or line.startswith("ID"):
                    continue
                (_, _, hap1, _, _) = line.split(',', 4)
                loci_in_file = [h.split('*')[0] for h in hap1.strip().split('~')]
                for locus in loci_in_file:
                    if locus not in self.all_loci_found:
                        self.all_loci_found.append(locus)
                break
            impute_outfile.close()

            # 1st pass: sum total freq per subject
            happair_id_total = {}
            impute_outfile = self.open_impute_file(filename)
            for line in impute_outfile:
                if line.strip() == "" or line.startswith("ID"):
                    continue
                (subject_id, _, _, _, freq) = line.split(',', 4)
                happair_id_total.setdefault(subject_id, 0)
                happair_id_total[subject_id] += float(freq)
            impute_outfile.close()

            # 2nd pass: process haplotype pairs
            impute_outfile = self.open_impute_file(filename)
            for line in impute_outfile:
                if line.strip() == "" or line.startswith("ID"):
                    continue
                (subject_id, rank, hap1, hap2, freq) = line.strip().split(',', 4)
                hap1_alleles = hap1.split('~')
                hap2_alleles = hap2.split('~')
                happair_freq = float(freq)
                prob = happair_freq / happair_id_total[subject_id]

                geno_locus = {}
                for idx, locus in enumerate(loci_in_file):
                    alleles = sorted([hap1_alleles[idx], hap2_alleles[idx]])
                    geno_locus[locus] = '+'.join(alleles)

                multiloc_glstring = '^'.join([geno_locus[locus] for locus in loci_in_file])
                self.multiloc_freq[subject_id][multiloc_glstring] = self.multiloc_freq[subject_id].get(multiloc_glstring, 0) + happair_freq

                for locus in loci_in_file:
                    self.gf_dicts[locus][subject_id][geno_locus[locus]] = self.gf_dicts[locus][subject_id].get(geno_locus[locus], 0) + prob

                # Build other combinations only if loci are present
                sevenloc = [l for l in ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1'] if l in loci_in_file]
                if sevenloc:
                    sevenloc_gl = '^'.join([geno_locus[l] for l in sevenloc])
                    self.other_freqs['sevenloc'][subject_id][sevenloc_gl] = self.other_freqs['sevenloc'][subject_id].get(sevenloc_gl, 0) + prob
                classI = [l for l in ['A', 'C', 'B'] if l in loci_in_file]
                if classI:
                    classI_gl = '^'.join([geno_locus[l] for l in classI])
                    self.other_freqs['classI'][subject_id][classI_gl] = self.other_freqs['classI'][subject_id].get(classI_gl, 0) + prob
                drdq = [l for l in ['DRB345', 'DRB1', 'DQA1', 'DQB1'] if l in loci_in_file]
                if drdq:
                    drdq_gl = '^'.join([geno_locus[l] for l in drdq])
                    self.other_freqs['DRDQ'][subject_id][drdq_gl] = self.other_freqs['DRDQ'][subject_id].get(drdq_gl, 0) + prob
                dr = [l for l in ['DRB345', 'DRB1'] if l in loci_in_file]
                if dr:
                    dr_gl = '^'.join([geno_locus[l] for l in dr])
                    self.other_freqs['DR'][subject_id][dr_gl] = self.other_freqs['DR'][subject_id].get(dr_gl, 0) + prob
                dq = [l for l in ['DQA1', 'DQB1'] if l in loci_in_file]
                if dq:
                    dq_gl = '^'.join([geno_locus[l] for l in dq])
                    self.other_freqs['DQ'][subject_id][dq_gl] = self.other_freqs['DQ'][subject_id].get(dq_gl, 0) + prob
            impute_outfile.close()

        # Helper to get top impute for any dictionary
        def top_impute_df(top_df, geno_dict, locus, which_impute):
            top_singleloc = pd.DataFrame()
            for id, values in geno_dict.items():
                top_genotype = max(values, key=values.get)
                top_freq = values[top_genotype]
                if which_impute == 'singleloc':
                    line = pd.DataFrame({'GENO_' + locus: top_genotype, 'GENO_' + locus + '_Prob': top_freq}, index=[id])
                elif which_impute == 'multiloc':
                    line = pd.DataFrame({'GLString': top_genotype, 'HapPair_Prob': top_freq}, index=[id])
                else:
                    line = pd.DataFrame({which_impute + '_GLString': top_genotype, which_impute + '_Prob': top_freq}, index=[id])
                top_singleloc = pd.concat([top_singleloc, line])
            top_df = pd.concat([top_df, top_singleloc], axis=1)
            return top_df

        # Build top impute DataFrame for multilocus
        top_multiloc_impute = pd.DataFrame()
        top_multiloc_impute = top_impute_df(top_multiloc_impute, self.multiloc_freq, '', 'multiloc')

        # Build top impute DataFrame for single loci present
        top_singleloc_impute = pd.DataFrame()
        for locus in self.all_loci_found:
            if any(self.gf_dicts[locus].values()):
                top_singleloc_impute = top_impute_df(top_singleloc_impute, self.gf_dicts[locus], locus, 'singleloc')

        # Build GLString for single locus genotype (SLUG) for present loci only (preserve order)
        if not top_singleloc_impute.empty:
            slug_cols = [f'GENO_{locus}' for locus in self.all_loci_found if f'GENO_{locus}' in top_singleloc_impute.columns]
            top_singleloc_impute['SLUG_GLString'] = top_singleloc_impute[slug_cols].agg('^'.join, axis=1)
            prob_cols = [f'GENO_{locus}_Prob' for locus in self.all_loci_found if f'GENO_{locus}_Prob' in top_singleloc_impute.columns]
            top_singleloc_impute = top_singleloc_impute[['SLUG_GLString'] + prob_cols]

        # Build top impute DataFrames for other combinations present
        other_top_imputes = []
        for which in ['sevenloc', 'classI', 'DRDQ', 'DR', 'DQ']:
            if self.other_freqs[which]:
                df = pd.DataFrame()
                df = top_impute_df(df, self.other_freqs[which], '', which)
                other_top_imputes.append(df)

        # Concatenate all top impute DataFrames
        dfs_to_concat = [df for df in [top_multiloc_impute, top_singleloc_impute] + other_top_imputes if not df.empty]
        if dfs_to_concat:
            top_impute = pd.concat(dfs_to_concat, axis=1)
            top_impute = top_impute.reset_index(names=['ID'])
            return top_impute
        else:
            return pd.DataFrame()

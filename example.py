
from voihla import ImputationPreprocessor, SingleLocusAnalysis, MultiLocusAnalysis, CalibrationPlotter, MonteCarloEpletAnalysis, EpletAnalysis
import pandas as pd

# Example script to demonstrate usage of the voihla package

# Import and preprocess imputation data
preprocessor = ImputationPreprocessor()
imputation_filename = 'imputation.csv'
imputation_data = preprocessor.process_files([imputation_filename])

# Load in truth data
truth_filename = 'truth_table.csv'
truth_data = pd.read_csv(truth_filename)

# Single Locus analysis
slug_analysis = SingleLocusAnalysis(truth_data, imputation_data)
slug_results = slug_analysis.get_results_df()
slug_results_A = slug_results['A']

# Single Locus Calibration Plot
plotter = CalibrationPlotter(n_bins=2)
plotter.calibration_plot(slug_results_A, locus='A', title='Single Locus A Calibration', save_path='Calibration_A.png')


# Multi Locus analysis
ml_analysis = MultiLocusAnalysis(truth_data, imputation_data)
ml_results = ml_analysis.get_results_df()
ml_results_AB = ml_results['multiloc']
# Calibration Plot for MUG
plotter.calibration_plot(ml_results_AB, locus='MUG', title='Multi Locus Calibration for HLA-A and HLA-B Imputation', save_path='Calibration_AB.png')

# Eplet Analysis
EpReg_api_key_file = 'api.key'
for line in open(EpReg_api_key_file):
    epReg_api_key = line
eplet_analysis = MonteCarloEpletAnalysis(api_key=epReg_api_key, delay=1.0) # This will use the API key when making eplets
random_pairs = eplet_analysis.create_random_pairs('truth_table.csv', 2, 42)

# Multiple Files for Imputation
preprocessor = ImputationPreprocessor()
imputation_files = ['../Imputation-Validation/impute.nyulowres.AFA.csv.gz',
                    '../Imputation-Validation/impute.nyulowres.API.csv.gz',
                    '../Imputation-Validation/impute.nyulowres.CAU.csv.gz',
                    '../Imputation-Validation/impute.nyulowres.HIS.csv.gz',
                    '../Imputation-Validation/impute.nyulowres.NAM.csv.gz']
imputation_data = preprocessor.process_files(imputation_files)

truth_filename = '../Imputation-Validation/genotype_truth_table.csv'
truth_data = pd.read_csv(truth_filename)
slug_analysis = SingleLocusAnalysis(truth_data, imputation_data)
slug_results = slug_analysis.get_results_df()

all_loci = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
plotter = CalibrationPlotter(n_bins=10)
for locus in all_loci:
    locus_results = slug_results[locus]
    plotter.calibration_plot(locus_results, locus=locus, title=f'Single Locus {locus} Calibration', save_path=f'Calibration_{locus}.png')

ml_analysis = MultiLocusAnalysis(truth_data, imputation_data)
ml_results = ml_analysis.get_results_df()
ml_results_AB = ml_results['multiloc']
# Calibration Plot for MUG
analysis_types = ['9-loci', '7-loci', 'HLA-Class I', 'DRDQ', 'DR', 'DQ']
plotter.calibration_plot(ml_results['multiloc'], locus='MUG', title=f'Multi Locus Calibration for 9-loci Imputation', save_path=f'Calibration_9-loci.png')
plotter.calibration_plot(ml_results['sevenloc'], locus='MUG', title=f'Multi Locus Calibration for 7-loci Imputation', save_path=f'Calibration_7-loci.png')
plotter.calibration_plot(ml_results['classI'], locus='MUG', title=f'Multi Locus Calibration for HLA-Class I Imputation', save_path=f'Calibration_HLA-Class_I.png')
plotter.calibration_plot(ml_results['DRDQ'], locus='MUG', title=f'Multi Locus Calibration for DRDQ Imputation', save_path=f'Calibration_DRDQ.png')
plotter.calibration_plot(ml_results['DR'], locus='MUG', title=f'Multi Locus Calibration for DR Imputation', save_path=f'Calibration_DR.png')
plotter.calibration_plot(ml_results['DQ'], locus='MUG', title=f'Multi Locus Calibration for DQ Imputation', save_path=f'Calibration_DQ.png')

check_pred = ml_analysis.pred_incorrect_high_prob('DR')
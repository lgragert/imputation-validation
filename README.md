# VOIHLA - Validation of Imputed HLA

Compute high resolution HLA imputation validation metrics using the `voihla` Python package and scikit-learn model evaluation statistics.

## Overview

The `voihla` package provides tools to preprocess, analyze, and visualize HLA imputation results. It supports single-locus, multilocus, and eplet-level analyses using standard metrics and calibration plots.

## Installation

Install dependencies using pip:
```
pip install -r requirements.txt
```

## Example Input Files

- `imputation.csv`: Imputation output with predicted HLA haplotype pairs and probabilities.
- `truth_table.csv`: High resolution genotype truth table in GLString format.

Example contents of files that are input for the package:

imputation.csv
```
ID,Rank,Hap1,Hap2,HapPair_Prob
D3505,1,A*30:02~B*14:02,A*32:01~B*39:10,0.3150459288416418
D3505,2,A*30:02~B*14:02,A*32:01~B*39:01,0.2673517305598033
D3505,3,A*30:02~B*39:10,A*32:01~B*14:02,0.09971243338882652
D3505,4,A*30:02~B*14:02,A*32:01~B*39:06,0.09552497014201682
D3505,5,A*30:01~B*14:02,A*32:01~B*39:10,0.0787155933156964
D3505,6,A*30:01~B*14:02,A*32:01~B*39:01,0.06679899077690125
D3505,7,A*30:01~B*39:10,A*32:01~B*14:02,0.02920649369703246
D3505,8,A*30:01~B*14:02,A*32:01~B*39:06,0.02386732857916924
D3505,9,A*30:02~B*14:02,A*32:01~B*39:24,0.009453181190469357
D3505,10,A*30:01~B*14:02,A*32:01~B*39:24,0.002361918368107546
D3505,11,A*30:02~B*39:01,A*32:01~B*14:02,0.002257484294942207
D13880,1,A*30:02~B*07:02,A*34:02~B*53:01,0.40269177048888066
D13880,2,A*30:01~B*07:02,A*34:02~B*53:01,0.20308144129546576
D13880,3,A*30:02~B*53:01,A*34:02~B*07:02,0.13918038201193
D13880,4,A*30:01~B*53:01,A*34:02~B*07:02,0.11366198792610206
D13880,5,A*30:02~B*07:05,A*34:02~B*53:01,0.04839353857099623
D13880,6,A*30:01~B*07:05,A*34:02~B*53:01,0.04136109507629823
D13880,7,A*30:02~B*53:01,A*34:02~B*07:09,0.011395011799743957
D13880,8,A*30:02~B*53:01,A*34:02~B*07:05,0.010829180412477053
D13880,9,A*30:01~B*53:01,A*34:02~B*07:09,0.009305763318635456
D13880,10,A*30:01~B*53:01,A*34:02~B*07:05,0.008843675778868338
D13880,11,A*30:04~B*53:01,A*34:02~B*07:02,0.0027297420425591353
```

`truth_table.csv`
```
ID,GLString
D3505,A*30:02+A*32:01^B*14:02+B*39:01
D13880,A*30:02+A*34:02+B*07:05+B*53:01
```

## Usage

All main modules are in the `voihla` folder.

### Preprocessing

Convert raw imputation files to analysis-ready format:

```python
from voihla.preprocessing import ImputationPreprocessor

preprocessor = ImputationPreprocessor()
top_impute = preprocessor.process_files(['imputation.csv'])
top_impute.to_csv('lowres_topprob_impute.csv', index=False)
```
### Single-Locus Analysis
```
import pandas as pd
from voihla.analysis import SingleLocusAnalysis
from voihla.preprocessing import ImputationPreprocessor

preprocessor = ImputationPreprocessor()
impute_df = preprocessor.process_files(['imputation.csv'])
truth_df = pd.read_csv('truth_table.csv')
analysis = SingleLocusAnalysis(truth_df, impute_df)
results = analysis.run()
print(results)
```

### Multilocus Analysis
```
from voihla.analysis import MultiLocusAnalysis

analysis = MultiLocusAnalysis(truth_df, impute_df)
results = analysis.run()
print(results)
```

### Calibration Plots
```
from voihla.plotting import CalibrationPlotter

plotter = CalibrationPlotter(n_bins=4)
locus = 'A'
df = analysis.get_results_df()[locus]
fig = plotter.calibration_plot(df['y_true'], df['y_prob'], f'Calibration {locus}')
fig.savefig(f'Calibration_{locus}.png')
```

### Eplet-Level Analysis
```
from voihla.eplet import MonteCarloEpletAnalysis

eplet_analysis = MonteCarloEpletAnalysis(api_key='YOUR_API_KEY')
pairs_df = eplet_analysis.create_random_pairs('truth_table.csv', n_pairs=100)
results_df = eplet_analysis.analyze_eplet_mismatches(pairs_df)
results_df.to_csv('DRDQ_eplet_lowres_impute100.csv', index=False)
```

### Output
- Calibration plots saved as PNG files
- ROC curves saved as PNG files
- Classification reports
- Summary CSV files

### API Reference
Please go to the Eplet Registry for an API key.


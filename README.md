# imputation-validation
Compute high resolution HLA imputation validation metrics using Python scikit-learn model evaluation statistics


## Overview

The objective is to write a Python script/module that computes imputation performance metrics using built-in functions the scikit-learn package.

### Two Input Files:

1) High Resolution Imputation Output - 9-locus HLA haplotype pairs and probabilities 
2) True High Resolution Genotype

### Procedure:

1) Decompose the predicted HLA and true HLA into several different categories 
	1) Start with single-locus HLA genotype, by locus. Start with just A, B, C, DRB1, DQB1.
	2) After you have program working, add tests for DRB3/4/5, DQA1, DPA1, DPB1, the whole multi-locus unphased genotype and then by amino acid position, then by eplet.
2) Assign the true positive, false positive, true negative, false negative cases for the category based on the most probable multilocus unphased genotype. Use a probability threshold for positive vs negative predictions of 50%. The truth is binary - either the high resolution genotype was correct or it was not - all alleles have to be correct. 
3) Compute prediction performance metrics using scikit-learn functions.


### Imputation output format:

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

### High Resolution Genotype truth table format:

```
ID,GLString
D3505,A*30:02+A*32:01^B*14:02+B*39:01
D13880,A*30:02+A*34:02+B*07:05+B*53:01
```

Genotype List String: a grammar for describing HLA and KIR genotyping results in a text string
https://pubmed.ncbi.nlm.nih.gov/23849068/

### Formatting for all scripts to run
All scripts mentioned below use the Genotype List String (GLString) format. Creating the truth table from NGS data is done in `NGS_glstring_format.py`. Formatting imputation data to GLString is done in `impute_glstring_format.py` and creates multiple files necessary for different levels of analysis.

Command line prompt: `python3 NGS_glstring_format.py NGS_file.csv`

The input it expects:
1. NGS columns needed for script:
```
unos,NGS_Ai,NGS_Aii,NGS_Bi,NGS_Bii,NGS_Ci,NGS_Cii,NGS_DRB1i,NGS_DRB1ii,NGS_DRB345i,NGS_DRB345ii,NGS_DQA1i,NGS_DQA1ii,NGS_DQB1i,NGS_DQB1ii,NGS_DPA1i,NGS_DPA1ii,NGS_DPB1i,NGS_DPB1ii

(where unos is the ID of the subject)
```

Output for NGS_glstring_format.py:
1. `genotype_truth_table.csv`

This script will be used as the truth table for every single analysis.

An example using above IDs.
```
ID,GLString
D3505,A*30:02+A*32:01^B*14:02+B*39:01
D13880,A*30:02+A*34:02^B*07:05+B*53:01
```

Command line prompt: `python3 impute_glstring_format.py imputefile1.csv.gz imputefile2.csv.gz imputefile3.csv.gz`

The input it expects:
1. Imputation format before GLString:
   The imputation files are gzip CSV files and the script can take several files as arguments at once. The above prompt is an example, but you can have more than 3 files or only 1 file. The expected information between each comma is below:
   
   Header:
```
ID,Rank,Hap1,Hap2,HapPair_Prob
```

Output for impute_glstring_format.py:
1. `lowres_topprob_impute.csv`

This takes the top probable imputation from each recipient for either SLUG or MUG analyses, but not for eplet analysis. It will have separate GLStrings for each analysis with the probability promptly after it.

Header:
```
ID,9loc_GLString,HapPair_Prob,SLUG_GLString,GENO_A_Prob,GENO_B_Prob,GENO_C_Prob,GENO_DRB345_Prob,GENO_DRB1_Prob,GENO_DQA1_Prob,GENO_DQB1_Prob,GENO_DPA1_Prob,GENO_DPB1_Prob,7loc_GLString,7loc_Prob,ClassI_GLString,ClassI_Prob,DRDQ_GLString,DRDQ_Prob,DR_GLString,DR_Prob,DQ_GLString,DQ_Prob
```
2. `lowres_*_impute_csv`, where * is DRDQ, DR, or DQ. For Class II eplet-level only and keeps all probable imputations rather than the top probable.

Header for each file:
```
ID,DRDQ_GLString,DRDQ_freq
ID,DR_GLString,DR_freq
ID,DQ_GLString,DQ_freq
```

### Single Locus Unphased Genotype (SLUG) level analysis:

Multiple rows have the pair of alleles for a single locus.

We make predictions on DPA1 and DPB1, but those loci can be ignored depending on if we have high resolution for them and want to evaluate those predictions.

Command line prompt:
```
python3 singleloc_analysis.py genotype_truth_table.csv lowres_topprob_impute.csv quant_num

(where quant_num=an integer for the number of quantiles you want to make for the calibration plots)
```

Output: 
```
Calibration_locus.png
ROC_AUC_9loc.png
Classification report print statements

(where locus=A,B,C,DRB345,DRB1,DQA1,DQB1,DPA1,DPB1)
```
Here are two example plots outputted. The commandline prompt would be: `python3 singleloc_analysis.py genotype_truth_table.csv lowres_topprob_impute.csv 4`

Calibration Plot for DQB1 locus:

![Calibration_DQB1](https://github.com/lgragert/imputation-validation/assets/114109419/f14c87b1-81cb-4b49-8f3b-62c2e2187a3f)

ROC Curve of DQB1 locus:

![ROC_DQB1](https://github.com/lgragert/imputation-validation/assets/114109419/4d2a86cc-3fe0-4f33-b93e-0d705746b774)


### Multilocus unphased genotype (MUG) level analysis:

Multiple rows have the same multilocus unphased genotype, but different phasing/haplotype arrangements. If you want to measure the ability to get the entire multilocus unphased genotype correct, you'd want to sum across those. We rarely know the haplotypes experimentally, unless we have family pedigree data.

This scripts makes plots for different combinations of loci such as 9-loci, 7-loci, Class I only (HLA-A,-C,-B), and DR-DQ only (HLA-DRB345, -DRB1, -DQA1, -DQB1)

Command line prompt:
```
python3 multiloc_analysis.py genotype_truth_table.csv lowres_topprob_impute.csv quant_num

(where quant_num=an integer for the number of quantiles you want to make for the calibration plots)
```

Output: 
```
Calibration_*.png

(where *=9loc, 7loc, ClassI, DRDQ)
```


### Eplet-level analysis:

Finally we could test quality of predictions by eplets, by converting the allele names using HLAGenie.

Eplet-Level Analysis:

We use the [Eplet Registry API](https://www.epregistry.com.br/), which can only handle 100 pairings at a time, so be cautious of that when creating these plots. Right now, scripts only analyzes Class II eplet mismatches.

**Notice**: Please have your own API key ready in a file to use these scripts. API key not provided in scripts. Please get in touch with the Eplet Registry for one.

If there are no donor-recipient pairings, then you can resort to Monte-Carlo Pairings simulation. 

Pipeline for the Monte-Carlo Pairings:
1. Create simulated donor-recipient pairings for both the truth table file and the imputation file.

Command line prompt:  
```
python3 DRDQ_pair_simulation.py genotype_truth_table.csv lowres_*_impute.csv

(where *=DRDQ, DR, DQ, depending on which analyses you choose to do)
```

Output: 
```
*_pairs_truth.csv
*_pairs_imputation.csv

(where *=DRDQ, DR, DQ)
```

2. Random sample the simulated dataset and run only those pairs in the eplet API. `API_key.key` is a file with the eplet registry secret API key.

   If you do not have an API key and cannot use the Eplet Registry API, then you may have to skip this step and continue to step 3 by doing the calculations yourself on the website.

Command line prompt: 
```
python3 montecarlo_pairings.py *_pairs_truth.csv *_pairs_imputation.csv * n_pairs API_key.key

(where *=DRDQ, DR, DQ, and n_pairs is an integer of how many pairings you want to keep)
```

Output:

   If you cannot use the API Registry, the format of the output should be simple if you can only use the website and need to calculate it by hand. 
   1. The first output is a CSV called `*_pairs_truth#.csv` where * = DR,DQ, or DRDQ and # = the number of pairs.
      ```
      HEADER FORMAT:
      DON_ID,DON_GLString,REC_ID,REC_GLString
      ```
   2. The second output is a CSV called `*_eplet_lowres_impute#.csv` with the * and # does the exact same as the first output.
      ```
      HEADER FORMAT:
      DON_ID,REC_ID,PairProb_DRDQ,DON_DRDQ,REC_DRDQ,ALL_quantity,ALL_details
      Where ALL_quantity and ALL_details are lists of the total unique eplets and eplet counts from the eplet calculation. 
      ```

4. Create calibration plots based on results from the eplet API.

Command line prompt: 
```
python3 eplet_MC_analysis.py *_pairs_truth#.csv *_eplet_lowres_impute#.csv quant_num * #

(where *=DR,DQ, or DRDQ and #=number of pairs, typically 100)
(where quant_num=number of quantiles you want on the calibration plots)
```

An example would be: `python3 eplet_MC_analysis.py DRDQ_pairs_truth100.csv DRDQ_eplet_lowres_impute100.csv 4 DRDQ 100`

The plot below would be comparing the DRDQ eplet mismatch counts for 100 pairings.

![Calibration_DRDQ_counts_100](https://github.com/lgragert/imputation-validation/assets/114109419/a54b54f9-5594-41f0-a24d-b653dded03ab)

Output: 
```
Calibration_*_counts_#.png
Calibration_*_eplets_#.png

(where *=DRDQ, DR, DQ, and #=number of pairs)
```

4. Validation of the eplet-level analysis to make sure calibration is accurate. Creates histogram plots and a CSV file.

Script: `eplet_validation.py`

Input:
```
DRDQ_eplet_lowres_impute#.csv
DR_eplet_lowres_impute#.csv
DQ_eplet_lowres_impute#.csv

(where #=number of pairs)
```

Output:
```
truth_impute_eplet_validation.csv
Histo_DRDQ_absdiff.png
Histo_DR_absdiff.png
Histo_DQ_absdiff.png
```

Here is an example histogram of comparing the DRDQ eplet mismatch counts for the 100 pairs we had previously.

![Histo_DRDQ_absdiff](https://github.com/lgragert/imputation-validation/assets/114109419/3995c7c7-6a63-44a2-b5b8-c5317cf75969)

To test exon 1 positions, we need to gather some typing that is true two-field resolution where exon 1 was typed, we can test to see if the probabilities are well calibrated, even position by position.  If it turns out as we expect that the predictions are poorly calibrated for the exon 1 positions, then we’ll know the TRS results for those positions are also flawed.  My hope is that the ARD positions have well-calibrated predictions even when imputing to two-field, then we won’t need to reimpute the SRTR file.


## The Average Typing Resolution Score (TRS) by HLA-Amino Acid Position

Create a table of amino acid-level TRS per each position by population.

Command line prompt:
```
typing_res_score_9loc_AA.py *

(Where *=one population group)
```

Input:
```
impute.*.csv.gz

(Where *=population group)
```

Output: CSV file which goes by each amino acid position per locus average for each population group.
```
HLA_AA_TRS_Average_*.csv

Where * is the population group.
Heading:
Subject_Type, Ethnicity, Locus, AA_Position, TRS_Average
```

### Manhattan Plots of the Average TRS for AA Positions by Population and Subject Type

Creates a visualization of the TRS by position in a Manhattan Plot.

Script: `AA_TRS_Manhattan_Plot.R`

Input: `SRTR_HLA_AA_TRS_Average_*.csv`, it takes all the CSV files at once.

Output:
```
AA_TRS_Manhattan_*subj.jpg
AA_TRS_Manhattan_*recips.jpg
AA_TRS_Manhattan_*donors.jpg

Where * = population group.
```

An example of all subject types for the CAU population group for all 9-loci.

![AA_TRS_Manhattan_CAUsubj](https://github.com/lgragert/imputation-validation/assets/114109419/3a7f79a3-5adb-4e0b-9749-a37d3de7195b)


If you would like to concatenate all Manhattan plots created, use ImageMagick.

In the command line, use [ImageMagick](https://imagemagick.org/script/download.php#google_vignette) [montage](https://imagemagick.org/script/montage.php) command to put the output together vertically so you can compare across population groups. 

Command line prompt:
```
magick montage AA_*#.jpg -tile 1x6 -geometry 5400x1200 AA_TRS_Manhat_%.jpg

* is just a wild card in this scenario to get all the jpg images.
Where # = (subj, donors, recips)
Where % = (SUBS, DONORS, RECIPS)
```


# Scikit-Learn Resources

The first step is set a probability threshold for positive vs negative predictions (50%) and then determine which cases are true positives, false positives, true negatives, false negatives to form a confusion matrix.

From there, scikit-learn can quickly generate all the statistics we could possibly need.

https://scikit-learn.org/stable/modules/model_evaluation.html
https://scikit-learn.org/stable/modules/classes.html#module-sklearn.metrics

See sidebar - Terminology and derivations from a Confusion matrix:

https://en.wikipedia.org/wiki/Receiver_operating_characteristic

## Probability Calibration curves


https://scikit-learn.org/stable/auto_examples/calibration/plot_calibration_curve.html

A lot of the probabilities for HLA imputation from antigen-level data fall in the midrange, so that means the most probable prediction was completely wrong quite a lot of the time. Even when that is the case, the predictions were extremely well calibrated, so this is what we aim for.

A lot of the standard model performance stats that are derived from the standard confusion matrix, like ROC, are much more dependent on the typing resolution of the input data than the calibration.

Having better typing resolution on the input will push more of the predictions to either end on the x-axis (towards 0% or 100% probability), but in our case we’re stuck with antigen-level typing. 

If you had 100 individuals with the same antigen-level HLA typing and ethnicity, and you went on and did the actual high resolution typing, the fact is you wouldn’t get the same high resolution genotype for all 100.

That means given the data we have there’s really no possible way for the top prediction from imputation to be the correct one every time, because the underlying HLA allelic diversity really is present in the population.

Many other published studies have measured imputation “accuracy” by computing the fraction of the time the top prediction was correct, and I argue it’s absolutely the wrong way to measure it.

But it is possible for imputation to provide high resolution genotype probabilities, many with mid-range values like 40%, 50%, 60%, that are in good calibration with this reality, and that’s what we hope to show.


See example of a well-calibrated DRB1 allele match prediction from Figure 1 Eberhard et al. 2010. Estimating unbiased haplotype frequencies from stem cell donor samples typed at heterogeneous resolutions: a practical study based on over 1 million German donors

https://pubmed.ncbi.nlm.nih.gov/20604895/

See also NMDP studies Madbouly et al. 2014 and Dehn et al 2016, and Sajulga et al. 2021 where calibration was not as good as the German study. The US population is more complex.

https://pubmed.ncbi.nlm.nih.gov/25040134/
https://pubmed.ncbi.nlm.nih.gov/27496216/
https://pubmed.ncbi.nlm.nih.gov/34362573/

## Brier score

Metric for evaluating the quality of the returned probabilities.

Mean squared difference between the predicted probability and the actual outcome.

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.brier_score_loss.html#sklearn.metrics.brier_score_loss

The value ranges between 0 and 1. The smaller the Brier score loss, the better, hence the naming with “loss”. 

https://en.wikipedia.org/wiki/Brier_score

An alternative option we've used is city block distance, weighted by size of the probability decile.

## F1 score

harmonic mean of the precision and recall. 1 is best and 0 is worst.

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html


## ROC Curve

Receiver-operating characteristic

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_curve.html

Plots true positive rate (TPR) on the Y axis and false positive rate (FPR) on the X axis.

Statistical power as a function of the Type I Error. Sensitivity or recall as a function of false positive rate.

A larger area under the curve (AUC) score (AKA c-statistic) is generally better.

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html

The “steepness” of ROC curves is also important, since it is ideal to maximize the TPR while minimizing the FPR.

One trick though is determining the threshold for calling "positive" vs "negative".

## Multiple Imputation

If we use Rubin’s Rules for multiple imputation, the association of HLA mismatch with transplant outcomes that we compute downstream will be unbiased, because we can propagate the uncertainty in HLA assignments throughout the statistical models.

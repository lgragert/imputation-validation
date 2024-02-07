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


### Single Locus Unphased Genotype (SLUG) level analysis:

Multiple rows have the pair of alleles for a single locus.

We make predictions on DPA1 and DPB1, but those loci can be ignored depending on if we have high resolution for them and want to evaluate those predictions.

### Multilocus unphased genotype (MUG) level analysis:

Multiple rows have the same multilocus unphased genotype, but different phasing/haplotype arrangements. If you want to measure the ability to get the entire multilocus unphased genotype correct, you'd want to sum across those. We rarely know the haplotypes experimentally, unless we have family pedigree data.
 
### Amino acid level or eplet-level analysis:

Finally we could test quality of predictions by single amino acid positions or eplets / amino acid motifs, by converting the allele names using HLAGenie.


To test exon 1 positions, we need to gather some typing that is true two-field resolution where exon 1 was typed, we can test to see if the probabilities are well calibrated, even position by position.  NYU and Penn have some data that is true two-field typing where exon 1 was covered.  If it turns out as we expect that the predictions are poorly calibrated for the exon 1 positions, then we’ll know the TRS results for those positions are also flawed.  My hope is that the ARD positions have well-calibrated predictions even when imputing to two-field, then we won’t need to reimpute the SRTR file.

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


# Datasets

## Deceased donor typing data from NYU

Input low resolution typing data:

```
/lustre/project/lgragert/imputation-validation/
n217_lowres.csv
```

Imputation output :

```
/lustre/project/lgragert/imputation-validation/
impute.nyulowres.*.csv.gz
```

True high resolution typings:

```
/lustre/project/lgragert/imputation-validation/
n217_withNGS.csv
```

This dataset has the advantage that it is real and in the solid organ setting, even though it is small.


## SRTR simulated dataset

For the SRTR dataset the input typing was antigen-level and high resolution genotypes are unknown. However we can simulate truth using weighted random choice to validate the imputation algorithm itself.

The quality of the reference population haplotype frequency estimates and alignment of the data we want to impute to the reference population data and Hardy-Weinberg equilibrium model can't be measured using this method.

To create the input files needed for the Python module:

1) Use the `impute.srtr.*.csv.gz` files for imputation output.
2) Do a weighted choice draw from the haplotype pair probability distribution to select a “true” high resolution multilocus unphased genotype.


## National Kidney Registry (NKR) Dataset

50K high resolution typings in kidney transplant setting.

Roll back HLA typing to antigen level to simulate how data would appear in SRTR then reimpute.

TODO - Get LG onto IRB to gain access to NKR data, then run imputation on NYU cluster.


## NMDP confirmatory typing validation datasets

Mostly stem cell donors with their initial recruitment typing and their true high resolution HLA typing, typically a customized typing confirmatory typing ordered on behalf of a searching patient. There are also some random prospective high resolution typings. 

TODO - Need a DUA to request dataset from Martin Maiers at NMDP. This is a large dataset (>>100,000 test cases) in stem cell transplant setting.


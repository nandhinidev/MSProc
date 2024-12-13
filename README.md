# Metabolomics/Proteomics Data Processing for Statistical Analysis
## Summary
This code processes mass spectrometry data (MS) using R package ['MALDIquant'](https://cran.r-project.org/web/packages/MALDIquant/index.html) and generates features for machine learning. The code executes the following steps in order:
* Pre-processing: Individual spectrum files are baselined, smoothed and internally recalibrated.
* Processing: Peaks are detected in each technical replicate spectrum and binned with a mass tolerance.
* Outlier detection: Average pairwise cosine similarity values are calculated for technical replicates of each sample and replicates that did not meet the reproducibility threshold (cosine simularity < 0.9) were removed.
* Average replicates: Average spectra of technical replicates for each sample.
* Post-processing: Peaks are detected in the averaged spectrum and binned with a mass tolerance.
* Imputation: Missing values are imputed and data is normalized using total ion current for each sample.
## R packages used
MALDIquant, MALDIquantForeign, tidyverse, DescTools, reshape2, grDevices, matrixStats, ggplot2, coop
## Input
## How to run the code?
## Output
A directory named 'output' is created inside 'proj' directory. Inside a sub-directory, the following files are saved:
* Directory with averaged spectra (data_avg) for each sample 
* Directories with individual peak lists for replicates (data_cent) and averaged spectra (data_avg_cent)
* List of outlier replicates (outliers.txt)
* Cosine similarity values plotted as a bar chart before and after removing outlier replicates
* Feature tables with m/z features in rows and samples in columns. 
  + All replicates (preprocessed_data.csv)
  + Averaged replicates (preprocessed_data_avg.csv)
  + Averaged replicates with imputation and normalization (preprocessed_data_avg_imp_norm.csv)


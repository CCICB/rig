
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rig

<!-- badges: start -->
<!-- badges: end -->

The goal of rig is to identify Radiation Induced Glioma (RIG) signatures
by analyzing gene expression data and looking for overexpression of RIG
genes within a given cohort.

rig bases classification on the proportion of RIG genes that are
overexpressed (based on robust Zscore approaches) so will only work if
your RIG samples make up \<50% of the cohort you supply.

If your RNAseq pipeline doesn’t compute expression for some of the genes
(e.g. the long noncoding RNAs), rig automatically adjusts the metrics
and should still perform sensibly.

## Installation

You can install the development version of rig from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("selkamand/rig")
```

## Example

Use the rig package to identify RIG signatures in your dataset:

``` r
library(rig)

# Load example data
data <- rig_example_data()

# Print the example data
print(head(data))
#>   sample_id   gene expression
#> 1   sample1   MYT1      16.49
#> 2   sample1  PCSK2      12.33
#> 3   sample1 KCNJ10      15.23
#> 4   sample1   PLP1      18.00
#> 5   sample1  OLIG2      21.00
#> 6   sample1  CROC4      20.00
```

``` r

# Identify RIG signatures in the example data
rig_identification <- rig_predict(data, col_samples = "sample_id", col_genes = "gene", col_expression = "expression")
#> ℹ Successfully tailored 49/59 gene names in the signature to your dataset
```

``` r

# Print the identification results
print(rig_identification)
#>            sample n_overexpressed_rig_genes prop_overexpressed_rig_genes
#> sample1   sample1                         7                          0.7
#> sample10 sample10                         0                          0.0
#> sample2   sample2                         9                          0.9
#> sample3   sample3                         1                          0.1
#> sample4   sample4                         0                          0.0
#> sample5   sample5                         0                          0.0
#> sample6   sample6                         0                          0.0
#> sample7   sample7                         0                          0.0
#> sample8   sample8                         0                          0.0
#> sample9   sample9                         0                          0.0
#>          median_rig_zscore predicted_to_be_rig
#> sample1          2.3642141                TRUE
#> sample10        -0.9702285               FALSE
#> sample2          5.3105009                TRUE
#> sample3          0.8554214               FALSE
#> sample4         -0.4189385               FALSE
#> sample5          0.1748680               FALSE
#> sample6         -0.2620046               FALSE
#> sample7         -0.5586053               FALSE
#> sample8          0.2695987               FALSE
#> sample9         -0.6090791               FALSE
```

# Plot Results

Plots coming soon

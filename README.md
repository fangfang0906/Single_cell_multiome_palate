# Single-cell multi-omics decodes regulatory programs during development of mouse secondary palate
The abnormal perturbation in gene regulation during palatogenesis may lead to cleft palate, a major congenital birth defect in humans and mice. However, **a comprehensive multi-omic map** of the developing secondary palate at single-cell resolution is **lacking**. In this study, we performed **single-cell multiome sequencing (10x Multiome)** and profiled chromatin accessibility (snATAC-seq) and gene expression (snRNA-seq) simultaneously within the same cells (n = 36,154) isolated from mouse secondary palate across embryonic days (E) 12.5, E13.5, E14.0, and E14.5. For more info please see our [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2022.11.02.514609v1.abstract).

<p align="center">
<img src="Figure 1_workflow.png">
</p>

<p align="center">
<img src="Summary.gif">
</p>

## Reproducibility
To reproduce the analysis and figures presented in our manuscript please see the [*Reproducibility*](https://github.com/lkmklsmn/empty_nn/tree/master/Reproducibility) folder.

## GEO dataset
Check out our jupyter notebook (in R environment) tutorial at [*EmptyNN - Cell Hashing Dataset Tutorial*](https://github.com/lkmklsmn/empty_nn/blob/master/tutorial/EmptyNN%20-%20Cell%20Hashing%20Dataset%20Tutorial.ipynb).

## Installation
EmptyNN is implemented in R and depends on the **keras** and **Matrix** R packages.

#### Option 1
```
$ git clone http://github.com/lkmklsmn/empty_nn
$ cd empty_nn

## enter R and install packages
$ R

> install.packages("EmptyNN_1.0.tar.gz", repos = NULL, type = "source")
```
#### Option 2
```
> install.packages("devtools")
> library(devtools)
> install_github("lkmklsmn/empty_nn")
```


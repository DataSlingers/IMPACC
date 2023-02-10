# IMPACC

`IMPACC` is a tool for unsupervised clustering and feature importance discovery. This document provides a tutorial of how to use `IMPACC`.

## Brief description of `IMPACC`
IMPACC (Interpretable Minipatch Adaptive Consensus Clustering) is a powerful methodology for consensus clustering using minipatch learning with adaptive feature and observation sampling schemes. IMPACC offers interpretable results by discovering features that differentiate clusters. This method is particularly applicable to sparse, high-dimensional data sets common in bioinformatics. MPCC (MiniPatch Consensus Clustering) provides consensus clustering by subsampling a tiny fraction of both observations and features at each iteration.

## Installation 
Install the package with the following code:
```{r}
library(devtools)
install_github("DataSlingers/IMPACC")
```

This package provide implementment of the `IMPACC` and  `MPCC` methods in R.

Gan, Luqin, and Genevera I. Allen. "Fast and Interpretable Consensus Clustering via Minipatch Learning." arXiv preprint arXiv:2110.02388 (2021).

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



## Tutorial 

###  Data Input
The input data format needs to be a matrix where columns are observations, rows are features, and cells are numerical values. If your input data is single-cell RNA-seq data, it needs to be normalized to achieve optimal results. 
```{r, message=FALSE, warning=FALSE}
library(IMPACC)
data(yan)
yan$raw[1:3, 1:3]
yan$sc_cnt[1:3, 1:3]
head(yan$sc_label)
K = length(unique(yan$sc_label))
```
The `yan` list contains three components: `raw` is gene expression count matrix before log2 transformation; `sc_cnt` is the log2 transformed gene expression matrix; and `sc_label` correspond to the cell labels provided by authors of the original publication. 

### Run IMPACC
```{r}
impacc = IMPACC(d=yan$sc_cnt,K = K,reps = 100,verbose=FALSE)
```

`IMPACC` returns a list containing ConsensusMatrix (numerical matrix),  labels (consensus class asssignments), feature_importance (feature), and nIter (stopping point).

Users can run MPACC (Minipatch Adaptive Consensus Clustering) by simply setting `adaptiveFeature` as FALSE. 

###  Construct different clustering results by passing the ConsensusMatrix argument to IMPACC_cluster() function 
```{r}
clus = IMPACC_cluster(ConsensusMatrix = impacc$ConsensusMatrix,K = 3)
```

###  Run MPACC with adaptive observation subsampling and random feature subsampling. 
```{r}
mpacc = IMPACC(d=yan$sc_cnt,K = K,adaptiveFeature = FALSE,verbose=FALSE)
```


###  Run MPCC with random minipatch subsampling. 
```{r}
mpcc = MPCC(d=yan$sc_cnt,K = K,verbose=FALSE)
```

###  Run IMPACC with multinomial feature evaluation
```{r}
impacc_multinomial = IMPACC(d=yan$sc_cnt,K = K,reps = 1,feature_evaluation = 'multinomial', verbose=FALSE)
```

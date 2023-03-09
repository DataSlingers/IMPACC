\name{IMPACC}
\alias{IMPACC}
\alias{MPCC}
\alias{IMPACC_cluster}

\title{run IMPACC}
\description{
  IMPACC function for consensus matrix and feature importance scores, using adaptive observation and/or feature subsampling;
  MPCC function for consensus matrix, using random feature and subsampling;
  IMPACC_cluster function for calculating clustering results from consensus matrix;
}

\usage{
IMPACC(d=NULL,
                  K=NULL,
                  adaptiveFeature = TRUE,
                  reps=300,
                  pItem=0.25,
                  pFeature=0.1,
                  innerLinkage="ward.D",
                  distance="manhattan",
                  h=0.95,
                  E= 3,
                  qI=0.95,
                  qF=1,
                  alpha_I=0.5,
                  alpha_F=0.5,
                  pp=0.05,
                  finalAlgorithm='hclust',
                  finalLinkage='ward.D',
                  early_stop=TRUE,
                  num_unchange = 5,
                  eps = 0.00001,
                  feature_evaluation = 'ANOVA',
                  seed=NULL,
                  verbose=TRUE)


MPCC(d=NULL,
                K = NULL,
                reps=300,
                pItem=0.25,
                pFeature=0.1,
                innerLinkage="ward.D",
                distance="manhattan",
                h=0.95,
                finalAlgorithm='hclust',
                finalLinkage='ward.D',
                early_stop=TRUE,
                num_unchange = 5,
                eps = 0.00001,
                seed=NULL,
              verbose=TRUE)


IMPACC_cluster(ConsensusMatrix=NULL,
              K=NULL,
              finalAlgorithm='hclust',
              finalLinkage='ward.D')
}

\arguments{
  \item{d}{data to be clustered. A data matrix where columns are observation and rows are features.}
  \item{K}{integer value. Number of clusters}
  \item{adaptiveFeature}{boolean. If TRUE, algorithm will adaptively select features. }
  \item{reps}{integer value. Maximum number of minipatches.  }
  \item{pItem}{numerical value. Proportion of items to sample.}
  \item{pFeature}{numerical value. Proportion of features to sample.  }
  \item{innerLinkage}{Heirarchical linkage method for minipatch.}
  \item{distance}{character value. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.   }
  \item{h}{numerical value. quantile cutoff of dendrogram height.  }
  \item{E}{interger value. the least number of times a feature is subsampled in burn-in stage}
  \item{qI}{numerical value. high uncertainty cutoff for observations}
  \item{qF}{numerical value. high importance cutoff for features}
  \item{alpha_I}{numerical value. learning rate for observation weight updates}
  \item{alpha_F}{numerical value. learning rate for feature weight updates}
  \item{pp}{numerical value. percentile cutoff for p-value of ANOVA test}
  \item{early_stop}{boolean. If TRUE, the algorithm will stop when the consensus matrix is stable.}
  \item{eps}{numeric value. Measure the change of stability of consensus matrix}
  \item{num_unchange}{integer value. Number of continous stable minipatches needed to stop the function. }
  \item{feature_evaluation}{character value. 'ANOVA': evaluate feature importance by ANOVA test, 'rankANOVA' evaluate feature importance by Kruskal-Wallis test, 'multinomial': evaluate feature importance by multinomial regression test. }
  \item{seed}{optional numerical value.  sets random seed for reproducible results.}
  \item{verbose}{ boolean. If TRUE, print messages to the screen to indicate progress.  This is useful for large datasets.}

  \item{ConsensusMatrix}{result of IMPACC/MPCC.}
  \item{finalAlgorithm}{character value. cluster algorithm.  'hc' heirarchical (hclust), 'spec' for spectral clustering,'pam' for paritioning around medoids, 'km' for k-means upon data matrix, 'ap' for affinity propagation clustering, 'convex' for convex clustering or a function that returns a clustering.  See example and vignette for more details. }
  \item{finalLinkage}{heirarchical linkage method for consensus matrix. }

}




\details{
IMPACC implements the Interpretable Minipatch Adaptive Consensus Clustering of Gan, L., & Allen, G. I. (2021). Its utility is to provide stable and robust clustering membership, interpretability in terms of feature importance in a time efficient manne. IMPACC takes a numerical data matrix of obervations as columns and features as rows.  This function construct minipatches according to pItem, pFeature, and clusters the data by fitting a hierarchical clustering with ward.D linkage and manhanttan distance. And this function perform adaptive subsampling schemes according to E, qI, qF, alpha_I, alpha_F and pp.

MPCC implements the Minipatch Consensus Clustering of Gan, L., & Allen, G. I. (2021), which is based on random minipatches.

For a detailed description of usage, output and images, see the vignette by: openVignette().
}

\value{
IMPACC returns a list containing ConsensusMatrix (numerical matrix),  labels (consensus class asssignments), feature_importance (feature), and nIter (stopping point).

MPCC returns a list containing ConsensusMatrix (numerical matrix), labels (consensus class asssignments), and nIter (stopping point).

IMPACC_cluster returns a list of size N contaning clustering labels.
}

\author{ Luqin Gan luqin_gan@rice.edu
Genevera I. Allen gallen@rice.edu
 }

\references{
Gan, L., & Allen, G. I. (2021). Fast and Interpretable Consensus Clustering via Minipatch Learning. arXiv preprint arXiv:2110.02388.

}

\examples{
data(yan)
impacc = IMPACC(d=yan$sc_cnt,K = K,reps = 100,verbose=FALSE)
}

\keyword{ methods }

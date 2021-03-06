% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crules.R
\name{ProjectSimUMAP}
\alias{ProjectSimUMAP}
\title{Project a normalized simulated vector of cell expression
onto UMAP space given an normalized expression matrix.
If these vectors are large (large number of genes), PCA
may be performed on both `exp.mat` and `x.prime` to reduce
computational load.}
\usage{
ProjectSimUMAP(exp.mat, x.prime)
}
\arguments{
\item{exp.mat}{original expression data (E^0)}

\item{x.prime}{predicted expression data (E^t)}
}
\value{
a data frame with umap coordinates for the cells
         and simulated states (x1,y1) -> (x2,y2) where each
         row corresponds to a cell in the original dataset.
}
\description{
Project a normalized simulated vector of cell expression
onto UMAP space given an normalized expression matrix.
If these vectors are large (large number of genes), PCA
may be performed on both `exp.mat` and `x.prime` to reduce
computational load.
}

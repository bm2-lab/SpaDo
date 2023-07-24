\name{InitialClustering}
\alias{InitialClustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Initial clustering of spatial transcriptomics data for single slice.
}
\description{
    Initial clustering of spatial transcriptomics data for single slice.
}
\usage{
InitialClustering(expression_profile, user_offered = FALSE, sample_information_user_offered = NULL, nfeatures = 2000, resolution = (ifelse(user_offered == FALSE, 2, FALSE)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
    A data frame showing gene expression profile. The column is cell/spot names. The row is gene names.
}
  \item{user_offered}{
    Logical. If users offer cell type lables of each cell/spot. Default is FALSE.
}
  \item{sample_information_user_offered}{
    A vector showing cell type lables of each cell/spot. The names of vector should be corresponded to the colnames of 'expression_profile'.
}
  \item{nfeatures}{
    Numeric. The number of selected highly variable genes in Seurat package. Default is 2000.
}
  \item{resolution}{
    Numeric. The resolution for clustering in Seurat package. Default is 2.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    \item{sample_information }{Vector. The cell type lables of each cell/spot. The names of vector should be cell/spot names.}
    \item{high_var_genes }{Vector. The selected highly variable genes.}
    \item{sce_seurat }{Seurat object.}
}

\author{
    Duan Bin
}

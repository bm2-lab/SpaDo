\name{SpatialPlot}
\alias{SpatialPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Draw clustering result and spatial distribution of clusters.
}
\description{
    Drawing tsne graph with cluster lable and cluster distribution in spatial coordinates.
}
\usage{
SpatialPlot(initial_clustering_result, sample_information_coordinate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{initial_clustering_result}{
    A list. The output of function 'InitialClustering()'. The list containing three parts: 'sample_information','sce_seurat' and 'high_var_genes'.
}
  \item{sample_information_coordinate}{
    The spatial coordinates of each cell/spot. A data frame with columns are 'X' and 'Y' coordinates, rows are cell's names.
}
}

\author{
    Duan Bin
}

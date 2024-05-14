\name{DomainPlot}
\alias{DomainPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Draw three graphs of detected domains for single slice.
}
\description{
    For single slice analysis, draw three graphs of detected domains for single slice, that is, umap of domains, hierarchical clustering tree, and spatial distribution of detected domains.
}
\usage{
DomainPlot(domain_hclust, distribution_distance, sample_information_coordinate = NULL, k = ncol(domain_hclust[[3]]), size = 1, shape = 19, aspect_ratio = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{domain_hclust}{
    A list. The output of function 'DomainHclust()'.
}
  \item{distribution_distance}{
    A matrix. The output of function 'DistributionDistance()' or the output of function "SpatialCellTypeDistribution()".
}
  \item{sample_information_coordinate}{
    The spatial coordinates of each cell/spot. A data frame with columns are 'X' and 'Y' coordinates, rows are cell's names.
}
  \item{k}{
    Numeric. The number of selected domains.
}
  \item{size}{
    Numeric. The size of each point in the graph. Default is 1.
}
  \item{shape}{
    Numeric. The shape of each point in the graph. Default is 19 (solid circle).
}
  \item{aspect_ratio}{
    Numeric. Length to width ratio of graph. Default is 1.
}
}


\author{
    Duan Bin
}

\name{DomainPlot_multiple}
\alias{DomainPlot_multiple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Draw umap of detected domains and hierarchical clustering tree for multiple slices.
}
\description{
    Draw umap of detected domains with spatially proximate cell type distribution profile and corresponding hierarchical clustering tree for multiple slices.
}
\usage{
DomainPlot_multiple(domain_hclust, distribution_distance, k = ncol(domain_hclust[[3]]), size = 1, shape = 19, aspect_ratio = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{domain_hclust}{
    A list. The output of function 'DomainHclust()'.    
}
  \item{distribution_distance}{
    A matrix. The output of function 'DistributionDistance()'. When cell number is too large(>30k), "cell_type_distribution_multiple$cell_type_distribution_combine" from function "SpatialCellTypeDistribution_multiple()" is recommanded here to save memory and improve speed.
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

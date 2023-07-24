\name{DistributionDistance}
\alias{DistributionDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Calculating the distance between spatially proximate cell type distribution of each cell/spot.
}
\description{
    Calculating the distance between spatially proximate cell type distribution of each cell/spot.
}
\usage{
DistributionDistance(cell_type_distribution, distance = c("jensen-shannon", "manhattan"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell_type_distribution}{
    A data frame showing spatially proximate cell type distribution of each cell/spot. The row is cell/spot name. The column is cell type name. 
}
  \item{distance}{
    Character. The metric to calculate distance. 'jensen-shannon' and 'manhattan' are supported. 'jensen-shannon' will get better performance with lower speed. 'manhattan' is faster with a little worse performance.
}
}
\value{
    A matrix showing distance between every two cells/spots.
}
\author{
    Duan Bin
}

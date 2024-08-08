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
DistributionDistance(cell_type_distribution, distance = c("JSD", "manhattan"),no_cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell_type_distribution}{
    A data frame showing spatially proximate cell type distribution of each cell/spot. The row is cell/spot name. The column is cell type name. 
}
  \item{distance}{
    Character. The metric to calculate distance. The metric used to calculate distance supports both 'JSD' and 'manhattan' in this version. The 'jensen-shannon' metric typically delivers better performance but at a slower speed, making it recommended for datasets with fewer than 30,000 cells. For larger datasets (more than 30,000 cells), 'manhattan' is recommended due to its significantly faster computation, albeit with slightly reduced performance.
}
  \item{no_cores}{
    Numeric. The cpu number to calculate distance when choosing "JSD". 
}
}
\value{
    A matrix showing distance between every two cells/spots.
}
\author{
    Duan Bin
}

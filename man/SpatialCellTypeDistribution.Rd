\name{SpatialCellTypeDistribution}
\alias{SpatialCellTypeDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    The calculation of spatially proximate cell type distributions for single slice.
}
\description{
    The calculation of spatially proximate cell type distributions for single slice.
}
\usage{
SpatialCellTypeDistribution(sample_information_coordinate, sequence_resolution = c("single_cell", "spot"), sample_information_cellType = NULL, sample_information_decon = NULL, neighbour_search_method = (ifelse(sequence_resolution == "single_cell", "KNN", "radius")), k = (ifelse(neighbour_search_method == "KNN", 30, FALSE)), r = (ifelse(neighbour_search_method == "radius", 2, FALSE)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_information_coordinate}{
    A data frame showing coordinates of cells/spots. The rows are cell/spot names, and the columns are coordinates "X" and "Y".
}
  \item{sequence_resolution}{
    Character. The resolution of spatial transcriptomics. "single_cell" means single-cell resolution. "spot" means non-single-cell resolution.
}
  \item{sample_information_cellType}{
    A vector showing cell type lables of each cell/spot. The names of vector should be cell/spot names.
}
  \item{sample_information_decon}{
    A data frame showing deconvolution result of spots. The column is cell type name. The row is spot name.
}
  \item{neighbour_search_method}{
    Character. The method to search neighbours for each cell/spot. If 'sequence_resolution' is 'single_cell', this parameter is set to 'KNN'. If 'sequence_resolution' is 'spot', this parameter is set to 'radius'.
}
  \item{k}{
    Numeric. The number of neighbors. Default is 30. It's useful only when 'neighbour_search_method' is set to 'KNN'.
}
  \item{r}{
    Numeric. The value of radius. Default is 2. It's useful only when 'neighbour_search_method' is set to 'radius'.
}
}

\value{
    A data frame showing spatially proximate cell type distributions for single slice.
}

\author{
    Duan Bin
}

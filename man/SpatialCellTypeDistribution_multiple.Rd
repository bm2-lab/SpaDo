\name{SpatialCellTypeDistribution_multiple}
\alias{SpatialCellTypeDistribution_multiple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    The calculation of spatially proximate cell type distributions for multiple slices.
}
\description{
    The calculation of spatially proximate cell type distributions for multiple slices.
}
\usage{
SpatialCellTypeDistribution_multiple(sample_information_coordinate_list, sequence_resolution = c("single_cell", "spot"), sample_information_cellType_list = NULL, sample_information_decon_list = NULL,r=2,k=30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_information_coordinate_list}{
    A list of data frames, each of which shows coordinates of cells/spots for each slice. The rows of data frame are cell/spot names, and the columns are coordinates "X" and "Y".
}
  \item{sequence_resolution}{
    Character. The resolution of spatial transcriptomics. "single_cell" means single-cell resolution. "spot" means non-single-cell resolution.
}
  \item{sample_information_cellType_list}{
    A list of vectors, each of which shows cell type lables of each cell/spot. The names of each vector should be cell/spot names.
}
  \item{sample_information_decon_list}{
    A list of data frames, each of which shows deconvolution result of spots. The column is cell type name. The row is spot name.
}
  \item{r}{
  	Numeric. The value of radius. Default is 2. 
}
  \item{k}{
    Numeric. The number of neighbors. Default is 30.
}
}

\value{
	\item{cell_type_distribution_combine }{A data frame showing spatially proximate cell type distributions for multiple slices. The column is cell type names. The row is cell/spot names.}
	\item{datasets_lable }{A vector showing slice lables of each cell/spot name corresponding to rownames of 'cell_type_distribution_combine'.}
}

\author{
    Duan Bin
}
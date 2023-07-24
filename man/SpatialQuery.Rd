\name{SpatialQuery}
\alias{SpatialQuery}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Spatial Domain annotation for new slices.
}
\description{
	Spatial Domain annotation for new slices with spatial reference.
}
\usage{
SpatialQuery(spatialReference_result, cellType_distri_query)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spatialReference_result}{
  	Alist. The output of SpatialReference() function.
  }
  \item{cellType_distri_query}{
  	Spatially adjacent cell type distribution of new slices.
  }
}
\value{
	A data frame showing the query results
}

\author{
	Bin Duan
}

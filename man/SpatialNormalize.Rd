\name{SpatialNormalize}
\alias{SpatialNormalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Normalization for Spatial Transcriptomics.
}
\description{
    Different normalization strategies for different spatial transcriptomics platforms.
}
\usage{
SpatialNormalize(expression_profile, ST_method = c("osmFISH", "merFISH", "STARmap", "Others"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_profile}{
    A data frame showing gene expression profile. The columns are cells, and rows are genes.
}
  \item{ST_method}{
    Spatial transcriptomics platforms. In theory, SpatialAlign supports all spatial transcriptomics platforms. For 'osmFISH','merFISH', and 'STARmap', we adopt normalization methods used in their original papers. For 'Others' platforms, we adopt default normalization in Seuat.
}
}

\value{
    A data frame showing normalized gene expression profile. The columns are cells, and rows are genes.
}

\author{
    Duan Bin
}

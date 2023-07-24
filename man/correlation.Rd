\name{Correlation}
\alias{Correlation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Calculating the correlation of each column of a matrix.
}
\description{
    Calculating the correlation of each column of a matrix.
}
\usage{
Correlation(matrix, matrix2 = NULL, method = c("pearson", "spearman", "cosine", "euclidean"), parallel_run = FALSE, cpu_num = ifelse(parallel_run == TRUE, 4, FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{
  A matrix with column is sample, row is feature.
}
  \item{matrix2}{
  Another matrix with column is sample, row is feature.
}
  \item{method}{
  The metric to calculate the correlation. Four metrics("pearson", "spearman", "cosine", "euclidean") are supported.
}
  \item{parallel_run}{
    Logical. Default is FALSE. If TRUE, run this function in parallel.
}
  \item{cpu_num}{
    Numeric. Default is 8. It works only when 'parallel_run' is TRUE.
}
}
\value{
    A matrix showing the correlation between each sample.
}

\author{
    Duan Bin
}
\name{DomainDistribution}
\alias{DomainDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculating domain composition for each slice.
}
\description{
  Calculating domain composition for each slice.
}
\usage{
DomainDistribution(hclust_result_df, datasets_lable_combine, domain_num = ncol(hclust_result_df))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hclust_result_df}{
  A data frame showing hierarchical clustering result. The first output of function 'DomainHclust()'.
}
  \item{datasets_lable_combine}{
  A vector showing slice lable for each cell/spot.
}
  \item{domain_num}{
  Numeric. Selecting the number of domain.
}
}
\value{
  \item{domain_distribution }{A data frame showing domain composition of all slices.}
  \item{domain_result_each_list }{A list showing domain composition of each slice.}

}

\author{
  Duan Bin
}

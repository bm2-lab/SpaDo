\name{DomainHclust}
\alias{DomainHclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Hierarchical clustering of spatially proximate cell type distribution.
}
\description{
    Hierarchical clustering of spatially proximate cell type distribution.
}
\usage{
DomainHclust(distribution_distance, autoselection = TRUE, auto_resolution = c(0, 1, 2, 3, 4), domain_num = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distribution_distance}{
    A matrix showing distance between every two cells/spots.
}
  \item{autoselection}{
    Logical. Default is TRUE, which means the function will select optimal domain number automatically.
}
  \item{auto_resolution}{
    Numeric. Different resolutions to select optimal domain number automatically which are consistent to parameter 'deepSplit' of cutreeDynamic() finction from R package dynamicTreeCut. Default is 0. This parameter is used only when 'autoselection' is TRUE.
}
  \item{domain_num}{
    Numeric. Users decide domain number. Default is 10. This parameter is used only when 'autoselection' is FALSE.
}
}

\value{

  \item{hclust_result_df }{A data frame of character showing hierarchical clustering result. Each character of data frame is domain name.}
  \item{hclust_result_model }{A hclust object obtained from function 'hclust', showing hierarchical clustering result.}
  \item{hclust_result_order }{A data frame of numeric showing hierarchical clustering result. Each numeic of data frame is the order in hierachical clustering result.}

}
\author{
    Duan Bin
}

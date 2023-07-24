\name{DomainHclust}
\alias{DomainHclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DomainHclust(distribution_distance, autoselection = TRUE, auto_resolution = c(0, 1, 2, 3, 4), domain_num = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distribution_distance}{
%%     ~~Describe \code{distribution_distance} here~~
}
  \item{autoselection}{
%%     ~~Describe \code{autoselection} here~~
}
  \item{auto_resolution}{
%%     ~~Describe \code{auto_resolution} here~~
}
  \item{domain_num}{
%%     ~~Describe \code{domain_num} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (distribution_distance, autoselection = TRUE, auto_resolution = c(0, 
    1, 2, 3, 4), domain_num = 10) 
{
    result_hc <- hclust(d = as.dist(distribution_distance), method = "ward.D2")
    if (autoselection) {
        require(dynamicTreeCut)
        auto_resolution <- auto_resolution[1]
        cluster_hc_dy <- cutreeDynamic(result_hc, distM = distribution_distance, 
            deepSplit = auto_resolution)
        cluster_hc_dy_num <- length(table(cluster_hc_dy))
    }
    else {
        cluster_hc_dy_num = domain_num
    }
    cluster_hc_all = cutree(result_hc, k = 1:cluster_hc_dy_num)
    cluster_hc_all_order <- cluster_hc_all[result_hc$order, ]
    cluster_hc_all_order_rename <- as.data.frame(cluster_hc_all_order)
    colnames(cluster_hc_all_order_rename) <- paste("SO_level", 
        1:cluster_hc_dy_num, sep = "_")
    cluster_hc_all_order_rename_c <- apply(cluster_hc_all_order_rename, 
        2, function(x) {
            return(as.character(x))
        })
    row.names(cluster_hc_all_order_rename_c) <- row.names(cluster_hc_all_order_rename)
    cluster_hc_all_order_rename_c <- as.data.frame(cluster_hc_all_order_rename_c)
    for (i in 1:ncol(cluster_hc_all_order_rename_c)) {
        cluster_hc_all_order_rename_c[, i] <- paste("Term", i, 
            cluster_hc_all_order_rename_c[, i], sep = "_")
    }
    return(list(hclust_result_df = cluster_hc_all_order_rename_c, 
        hclust_result_model = result_hc, hclust_result_order = cluster_hc_all_order))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
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
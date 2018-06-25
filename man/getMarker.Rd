\name{getMarker}
\alias{getMarker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marker selection
}
\description{
Selecting markers if gene pairs are mutual nearest coexpressed in a number of cells.

}
\usage{
getMarker(filterres, MNN, MNNIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filterres}{
The list oject from function GeneFilter
}
  \item{MNN}{
The maximum number that how many nearest neighbors would be considered.
}
  \item{MNNIndex}{
The minimum number which is the coexpressed cells for each gene pairs
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\res{data}{The raw expresion data after filter widly expressed genes}
\res{geneSumm}{A data frame including gene names as, corresponding number of expressed cells and the peak number of probability density distribution function}
\res{cellSumm}{A data frame including cell ID and corresponding number of expressed genes}
\res{binadata}{The binary expression matrix corresponding to raw data}
\res{marker}{The filnal marker selected}

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
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
filterres=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2,cutoff=2)
filterres=GeneFilter(filterres=filterres,maxexp=dim(filterres$data)[2]*0.8)
filterres=getMarker(filterres=filterdata,MNN=200,MNNIndex=20)

## The function is currently defined as
function (filterres, MNN, MNNIndex) 
{
    data = filterres$data
    binadata = filterres$binadata
    genename = row.names(binadata)
    geneindex = rowSums(binadata)
    HamD = tcrossprod(binadata)
    diag(HamD) = 0
    marker = getMNN(HamD = HamD, genename = genename, MNN = MNN, 
        MNNIndex = MNNIndex)
    filterres$marker = marker
    return(filterres)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

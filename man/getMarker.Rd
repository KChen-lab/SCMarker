\name{getMarker}
\alias{getMarker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marker selection
}
\description{
Selecting markers if gene pairs are mutual nearest coexpressed or exclusivity in a number of cells.

}
\usage{
getMarker(obj, MNN, MNNIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
The list oject from function GeneFilter
}
  \item{MNN}{
The maximum number that how many nearest neighbors would be considered.
}
  \item{MNNIndex}{
The minimum number which is the coexpressed/mutually exclusivie cells for each gene pairs
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
{data}{The raw expresion data after filter widly expressed genes}
{geneSumm}{A data frame including gene names as, corresponding number of expressed cells and the peak number of probability density distribution function}
{cellSumm}{A data frame including cell ID and corresponding number of expressed genes}
{binadata}{The binary expression matrix corresponding to raw data}
{marker}{The filnal marker selected}
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
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
res=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2,cutoff=2)
res=GeneFilter(obj=res)
res=getMarker(obj=res,MNN=200,MNNIndex=20)

## The function is currently defined as
function (obj, MNN, MNNIndex)
{
    data = obj$newdata
    binadata = obj$binadata
    genename = row.names(binadata)
    geneindex = rowSums(binadata)
    HamD = tcrossprod(binadata)
    diag(HamD) = 0
    HamDD = tcrossprod((1 - binadata), binadata)
    MNNmarker = getMNN(HamD = HamD, genename = genename, MNN = MNN,
        MNNIndex = MNNIndex)
    MENmarker = getMEN(HamD = HamDD, genename = genename, MNN = MNN,
        MNNIndex = MNNIndex)
    obj$marker = union(MNNmarker, MENmarker)
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

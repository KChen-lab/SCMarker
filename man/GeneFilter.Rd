\name{GeneFilter}
\alias{GeneFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filtering gene
}
\description{
Filtering genes which are widely expressed in most of the cell population.
}
\usage{
GeneFilter(filterres, maxexp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filterres}{
The list object from function ModalFilter
}
  \item{maxexp}{
The maxmum number of expressed sampels. That means the gene would be filtered if it expressed in more then the maxexp number of cells.
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
ModalFilter
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
res=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2,cutoff=2)
res=GeneFilter(filterres=res,maxexp=dim(res$data)[2]*0.8)

## The function is currently defined as
function (filterres, maxexp) 
{
    data(excludeGene)
    excludeGene = as.character(excludeGene[, 1])
    geneSumm = filterres$geneSumm
    index = match(geneSumm$gene, excludeGene)
    geneSumm = geneSumm[is.na(index), ]
    data = filterres$data
    index = match(geneSumm$gene, row.names(data))
    data = data[index, ]
    binadata = filterres$binadata
    binadata = binadata[index, ]
    geneindex = rowSums(binadata)
    binadata = binadata[geneindex < maxexp, ]
    data = data[geneindex < maxexp, ]
    geneSumm = geneSumm[geneindex < maxexp, ]
    filterres$data = data
    filterres$geneSumm = geneSumm
    filterres$binadata = as.matrix(binadata)
    return(filterres)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

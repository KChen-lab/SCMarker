\name{ModalFilter}
\alias{ModalFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filtering genes with unimodal distribution
}
\description{
Filtering genes(cells) that expressed (non zero) distribution is similar with normal distribution.
}
\usage{
ModalFilter(data, geneK, cellK, width = 1, cutoff = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The raw continuous expression matrix
}
  \item{geneK}{
The least number of genes which expressed in a cell
}
  \item{cellK}{
The least number of cells which have expression for a gene
}
  \item{width}{
The numeric width which is used to estimate probability density distribution function, default 1.
}
  \item{cutoff}{
The numeric width which is used to estimate probability density distribution function, default 2.
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
\res{data}{The raw expresion data after filter genes or cells which less than cellK or geneK}
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
peak, genepeak
}
\examples{
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
res=ModalFilter(data=melanoma1,geneK=10,cellK=10)

## The function is currently defined as
function (data, geneK, cellK, width = 1, cutoff = 2) 
{
    rawdata = data
    cellSumm = data.frame(cell = colnames(data))
    cellSumm$count <- sapply(1:dim(data)[2], genecount, data = data, 
        index = "col")
    data = data[, cellSumm$count > cellK]
    binadata = Bina(data, cutoff = cutoff)
    geneSumm = data.frame(gene = row.names(data))
    geneSumm$count <- rowSums(binadata)
    data = data[geneSumm$count > geneK, ]
    binadata = binadata[geneSumm$count > geneK, ]
    geneSumm = geneSumm[geneSumm$count > geneK, ]
    cellSumm = cellSumm[cellSumm$count > cellK, ]
    geneSumm$exppeak = apply(data, 1, genepeak, width = width)
    data = data[geneSumm$exppeak >= 2, ]
    binadata = binadata[geneSumm$exppeak >= 2, ]
    geneSumm = geneSumm[geneSumm$exppeak >= 2, ]
    obj = list(rawdata = rawdata, newdata = data, geneSumm = geneSumm, 
        cellSumm = cellSumm, binadata = binadata)
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

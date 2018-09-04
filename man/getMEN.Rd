\name{getMEN}
\alias{getMEN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mutually nearest exclusive gene pairs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMEN(HamDD, genename, MNN, MNNIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HamDD}{
A matrix from multiplicate of binary data tcrossprod((1 - binadata), binadata).
}
  \item{genename}{
The rownames and colnames of matrix HamDD
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A list including the MNN nearest neighbors which mutually exclusive in more than MNNindex samples.
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
data = res$newdata
binadata = res$binadata
genename = row.names(binadata)
HamDD = tcrossprod((1 - binadata), binadata)
marker = getMEN(HamDD = HamDD, genename = genename, 
        MNN = 200, MNNIndex = 20)

## The function is currently defined as
function (HamDD, genename, MNN, MNNIndex) 
{
    MNNgene = lapply(1:dim(HamDD)[1], RankGene, MNN = MNN, HamD = HamDD, 
        geneName = genename, MNNIndex = MNNIndex)
    genePair = unique(as.character(do.call(cbind, lapply(1:length(MNNgene), 
        MNNpair, MNNgene = MNNgene, geneName = genename))))
    return(genePair)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

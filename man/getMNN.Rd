\name{getMNN}
\alias{getMNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mutual nearest coexpressed gene pairs
}
\description{
Idenitify mutual nearest expressed gene pairs based on binary expression data
}
\usage{
getMNN(HamD, genename, MNN, MNNIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HamD}{
A matrix from multiplicate of binary data tcrossprod(binadata) and diag are zero.
}
  \item{genename}{
The rownames and colnames of matrix HamD
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
A list including the MNN nearest neighbors which co-expressed in more than MNNindex samples.
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
getMarker
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
filterres=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2,cutoff=2)
filterres=GeneFilter(filterres=filterres,maxexp=dim(filterres$data)[1]*0.8)
data = filterres$data
binadata = filterres$binadata
genename = row.names(binadata)
HamD = tcrossprod(binadata)
diag(HamD) = 0
marker = getMNN(HamD = HamD, genename = genename, 
        MNN = 200, MNNIndex = 20)

## The function is currently defined as
function (HamD, genename, MNN, MNNIndex) 
{
    RankGene <- function(k, MNN, HamD, geneName) {
        x = HamD[k, ]
        xrank = rank(-x)
        MNNgene = geneName[xrank <= MNN & x > MNNIndex]
        return(MNNgene)
    }
    MNNgene = lapply(1:dim(HamD)[1], RankGene, MNN = MNN, HamD = HamD, 
        geneName = genename)
    MNNpair <- function(k, MNNgene, geneName) {
        subgene = MNNgene[[k]]
        index = match(subgene, geneName)
        PP <- function(i, index, MNNgene, k, geneName) {
            if (geneName[k] \%in\% MNNgene[[index[i]]]) 
                return(geneName[index[i]])
        }
        if (length(index[!is.na(index)]) > 0) {
            a = do.call(cbind, lapply(1:length(index), PP, index = index, 
                MNNgene = MNNgene, k = k, geneName = geneName))
            if (!is.null(a)) {
                return(a)
            }
        }
    }
    genePair = unique(as.character(do.call(cbind, lapply(1:length(MNNgene), 
        MNNpair, MNNgene = MNNgene, geneName = genename))))
    return(genePair)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

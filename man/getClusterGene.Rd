\name{getClusterGene}
\alias{getClusterGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying cluster specific genes under clutering from Seurat or dbscan
}
\description{
Identifying cluster specific genes under clutering from Seurat or dbscan
}
\usage{
getClusterGene(obj, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
The list objective from the function SCcluster
}
  \item{method}{
clutering results from "Seurat" or "dbscan"
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
\res{clustergene}{A data frame including cluster specfic genes}
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
res=ModalFilter(data=melanoma1,geneK=10,cellK=10)
res=GeneFilter(obj=res)
res=getMarker(obj=res,MNN=200,MNNIndex=20)
res=SCcluster(obj=res)
res=getClusterGene(obj=res,method="Seurat")

## The function is currently defined as
function (obj, method) 
{
    marker = obj$marker
    if (method == "Seurat") {
        seuratCluster = obj$seuratCluster
        gene = get_marker_genes(obj$rawdata[marker, names(seuratCluster)], 
            seuratCluster)
        gene$gene = marker
        obj$clustergene = gene
        obj$method = method
    }
    else if (method == "dbscan") {
        dbscanCluster = obj$dbscanCluster
        gene = get_marker_genes(obj$rawdata[marker, names(dbscanCluster)], 
            dbscanCluster)
        gene$gene = marker
        obj$clustergene = gene
        obj$method = method
    }
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

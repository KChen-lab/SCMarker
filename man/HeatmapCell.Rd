\name{HeatmapCell}
\alias{HeatmapCell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap of single cell based on top cluster specific genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HeatmapCell(obj, top)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
The list objective from the function getClusterGene
}
  \item{top}{
The number of gene which showed top cluster specificity
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
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
res=ModalFilter(data=melanoma1,geneK=10,cellK=10)
res=GeneFilter(obj=res)
res=getMarker(obj=res,MNN=200,MNNIndex=20)
res=SCcluster(obj=res)
res=getClusterGene(obj=res,method="Seurat")
HeatmapCell(obj=res,10)
## The function is currently defined as
function (obj, top) 
{
    feature = obj$clustergene
    topmarker = as.data.frame(feature \%>\% group_by(clusts) \%>\% 
        top_n(top, auroc))
    topmarker = topmarker[order(topmarker$clusts), ]
    method = obj$method
    if (method == "Seurat") {
        cluster = obj$seuratCluster
        cluster = cluster[order(cluster)]
    }
    else if (method == "dbscan") {
        cluster = obj$dbscanCluster
        cluster = cluster[order(cluster)]
    }
    clustertype = unique(cluster)
    data = obj$rawdata[topmarker$gene, names(cluster)]
    clustercount = table(cluster)
    gapIndex = c()
    for (i in 2:(length(clustertype) - 1)) {
        gapIndex = c(gapIndex, sum(clustercount[1:i]))
    }
    pheatmap(data, cluster_cols = F, cluster_rows = F, show_colnames = F, 
        gaps_col = gapIndex)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{HeatmapCluster}
\alias{HeatmapCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap of cluster based on top cluster specific genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HeatmapCluster(obj, top)
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
HeatmapCluster(obj=res,top=10)


## The function is currently defined as
function(obj,top,scale="none"){
    feature=obj$clustergene
    topmarker=as.data.frame(feature %>% group_by(clusts) %>% top_n(top, auroc))
    topmarker=topmarker[order(topmarker$clusts),]
    method=obj$method
    if (method=="Seurat"){
        cluster=obj$seuratCluster
    }else if (method=="dbscan"){
        cluster=obj$dbscanCluster
    }
    clustertype=levels(factor(cluster))
    genemean=c()
    for (i in clustertype){
        cell=names(cluster[cluster==i])
        subdata=obj$rawdata[topmarker$gene,cell]
        genemean=cbind(genemean,apply(subdata,1,mean))
    }
    row.names(genemean)=topmarker$gene
    colnames(genemean)=clustertype

    mat_col <- data.frame(cluster=levels(factor(cluster)))
    row.names(mat_col) <- levels(factor(cluster))
    mat_colors <- list(cluster=rainbow(length(clustertype)))
    names(mat_colors$cluster)=levels(factor(cluster))
    pheatmap(genemean,cluster_cols = F,scale=scale,cluster_rows = F,show_colnames=F,annotation_col = mat_col,annotation_colors = mat_colors)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

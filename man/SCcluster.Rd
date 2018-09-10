\name{SCcluster}
\alias{SCcluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Application in Seurat cluster using markers selected from SCMarker
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCcluster(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
The list objective from the function getMarker
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
{tsne}{tSNE embeding from Seurat}
{seuratCluster}{Cluster results from Seurat}
{dbscanCluster}{Cluster results from dbscan}
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
## The function is currently defined as
function (obj)
{
    data = obj$rawdata
    marker = obj$marker
    pbmc <- CreateSeuratObject(raw.data = data, min.cells = 3,
        min.genes = 200, project = "project")
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",
        scale.factor = 10000)
    pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean,
        dispersion.function = LogVMR, x.low.cutoff = 0.0125,
        x.high.cutoff = 3, y.cutoff = 0.5)
    pbmc <- ScaleData(object = pbmc)
    pbmc <- RunPCA(object = pbmc, pc.genes = gene, do.print = TRUE,
        pcs.print = 1:5, genes.print = 5)
    pbmc <- FindClusters(object = pbmc, reduction.type = "pca",
        dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE,
        force.recalc = TRUE)
    pbmc <- RunTSNE(object = pbmc, dims.use = 1:15, do.fast = TRUE)
    TSNE = pbmc@dr$tsne@cell.embeddings
    seuratCluster = pbmc@ident
    dbscanCluster = dbscan(TSNE, eps = 1.2, minPts = 15)$cluster
    names(dbscanCluster) = row.names(TSNE)
    obj$tsne = TSNE
    obj$seuratCluster = seuratCluster
    obj$dbscanCluster = dbscanCluster
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{SCMarker-package}
\Rdversion{2.0}
\alias{SCMarker-package}
\alias{SCMarker}
\docType{package}
\title{
\packageTitle{SCMarker}
}
\description{
\packageDescription{SCMarker}
}
\details{

\tabular{ll}{
Package: \tab SCMarker\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2018-09-04\cr
License: \tab GPL(>=2)\cr
}
}
\author{
\packageAuthor{SCMarker}

Maintainer: \packageMaintainer{SCMarker}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
res=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2)
res=GeneFilter(obj=res)
res=getMarker(obj=res,MNN=300,MNNIndex=30)
library(Seurat)
library(SC3)
library(dbscan)
library(dplyr)
res=SCcluster(obj=res)
res=getClusterGene(obj=res,method="Seurat")
HeatmapCluster(obj=res,top=10)
HeatmapCell(obj=res,5)
}

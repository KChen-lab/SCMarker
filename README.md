# SCMarker

SCMarker performs cell-type-specific marker selection from single cell RNA sequenceing data. 
It provides users a tool for selecting features from tens of thousands of genes for further cell-type clustering analysis. 

SCMarker is done based on two hypotheses: 
First, the expression of a gene should follow bi/multi-modal distribution in a mixed cell population if it is a marker of a specific cell-type. 
Second, marker genes of a cell type express synergistically in a subset of cells.

Developer
------------
Fang Wang (fwang9@mdanderson.org)


Marker selection
---------------------
The three main functions for this package are `ModalFilter()`, `GeneFilter()` and `getMarker()`. 
ModalFilter() performs the initial filter based on the least expressed number of genes(cells) and whether the gene has unimodal distribution. 
GeneFilter() takes the output of ModalFilter() and filters out genes that have unimodal distributed expressions and are expressed in more than maxexp cells.
getMarker() takes the output of GeneFilter() and selects the final markers based on synergistically (co- or mutual-exclusively) expressed gene pairs. 



Installation
----------------------
Download SCMarker_2.0.tar.gz
```R
install.packages("SCMarker_2.0.tar.gz",repos=NULL,type="source")
```
or install through GitHub
```R
library(devtools)
install_github("KChen-lab/SCMarker")
```


Usage
----------------------

```R
library(SCMarker)
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
res=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2)# default width = 1 for UMI data, width =2 for TPM data.
res=GeneFilter(obj=res)
res=getMarker(obj=res,k=300,n=30)
head(res$marker)


##Integrating with other analyses
library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(pheatmap)
library(Seurat)
library(dbscan)
res=SCcluster(obj=res)
res=getClusterGene(obj=res,method="Seurat")
HeatmapCluster(obj=res,top=10)
HeatmapCell(obj=res,top=10)
```

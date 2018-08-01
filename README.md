# SCMarker

SCMarker is a R package that performs ab initio marker selection from single cell RNA sequenceing data to enhance downstream cell-type clustering, trajectory inference and cell-type specific biological analysis.

SCMarker tests two properties of a scRNAseq dataset:

1. A gene expression level should follow bi/multi-modal distribution in a mixed cell population if it is a marker of a specific cell type.

2. Marker genes of the same cell type tend to be co-expressed in the same cells and mutually exclusive with genes which are markers of different cell type.

SCMarker is written by Fang Wang (fwang9@mdanderson.org) in Ken Chen' lab at the Department of Bioinformatics and Computational Biology in the MD Anderson Cancer Center


Marker selection
---------------------
SCMarker contains three main functions:

1. `ModalFilter()`,
    performs the initial filtering based on the least expressed number of genes (cells) and whether the gene is unimodal distribution.

2. `GeneFilter()`
    takes the output of `ModalFilter()` and filtered the genes which are unimodal distribution and express in more than maxexp cells.

3. `getMarker()`.
    takes the output of `GeneFilter()` and selects the final markers based on mutually coexpressed gene pairs.



Usage
----------------------
You can install SCMarker from GitHub
```
library(devtools)
install_github("fang0828/SCMarker")
```
Alternatively, you can download from GitHub
```
#https://github.com/Fang0828/SCMarker/releases
#Download SCMarker_1.1.tar.gz
install.packages("SCMarker_1.1.tar.gz",repos=NULL,type="source")
```

```R
library(SCMarker)
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
filterres=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2,cutoff=2)
filterres=GeneFilter(filterres=filterres,maxexp=dim(melanoma1))
filterres=getMarker(filterres=filterres,MNN=200,MNNindex=20)
filterres$marker
```

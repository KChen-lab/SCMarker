# SCMarker

SCMarker performs cell-type-specific marker selection from single cell RNA sequencing data.
It provides users a tool for selecting features from tens of thousands of genes for further cell-type clustering analysis.

SCMarker is done based on two hypotheses:
1) The expression of a gene should follow bi/multi-modal distribution in a mixed cell population if it is a marker of a specific cell-type.
2) Marker genes of a cell type express synergistically in a subset of cells.

Developer
------------
Fang Wang (fwang9@outlook.com)


Marker selection
---------------------
The three main functions for this package are `ModalFilter()`, `GeneFilter()` and `getMarker()`.

`ModalFilter()` performs the initial filter based on the least expressed number of genes(cells) and whether the gene has unimodal distribution.

`GeneFilter()` takes the output of ModalFilter() and filters out genes that have unimodal distributed expressions and are expressed in more than maxexp cells.

`getMarker()` takes the output of GeneFilter() and selects the final markers based on synergistically (co- or mutual-exclusively) expressed gene pairs.



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

```

An example to show how SCMarker [improve identification of NK cell in GBM data.](https://github.com/KChen-lab/SCMarker/blob/master/test/NK%20cell%20identification%20from%20GBM%20data.pdf)

Publication
-----------------------
Wang, Fang, et al. "SCMarker: ab initio marker selection for single cell transcriptome profiling." PLoS computational biology 15.10 (2019): e1007445.

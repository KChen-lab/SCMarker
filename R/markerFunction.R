# TODO: Add comment
#
# Author: FWang9
###############################################################################


####preprocess data
##count expression cell or gene
genecount<-function(k,data,index){
	if (index=="row"){
		x=data[k,]
	}else{
		x=data[,k]
	}
	return(sum(x!=0))
}
####peak of density
peak<-function(pp){
	y=pp$y
	x=pp$x
	Hpeak=c()
	for (i in 2:(length(y)-1)){
		if (y[i]>=y[i-1]&y[i]>y[i+1]){
			Hpeak=rbind(Hpeak,c(x[i],y[i]))
		}
	}
	if (length(Hpeak)==0){
		Hpeak=rbind(Hpeak,c(0,0))
	}
	Hpeak=as.data.frame(Hpeak)
	names(Hpeak)=c("H","Density")
	return(Hpeak)
}
###number of peak
genepeak<-function(x,width=0.05){
	if (length(x[x!=0])>1){
		pp=density(x[x!=0],bw=width)
		Hpeak=peak(pp)
		return(dim(Hpeak)[1])
	}else{
		return(0)
	}
}

#####filter data
ModalFilter<-function(data,geneK,cellK,width = 0.05,cutoff = 2){
	rawdata=data
	cellSumm=data.frame(cell=colnames(data))
	cellSumm$count<-sapply(1:dim(data)[2],genecount,data=data,index="col")
	data=data[,cellSumm$count>cellK]
	binadata=Bina(data,cutoff=cutoff)
	geneSumm=data.frame(gene=row.names(data))
	geneSumm$count<-rowSums(binadata)
	data=data[geneSumm$count>geneK,]
	binadata=binadata[geneSumm$count>geneK,]
	geneSumm=geneSumm[geneSumm$count>geneK,]
	cellSumm=cellSumm[cellSumm$count>cellK,]
	geneSumm$exppeak=apply(data,1,genepeak,width=width)
	data=data[geneSumm$exppeak>=2,]
	binadata=binadata[geneSumm$exppeak>=2,]
	geneSumm=geneSumm[geneSumm$exppeak>=2,]
	obj=list(rawdata=rawdata,newdata=data,geneSumm=geneSumm,cellSumm=cellSumm,binadata=binadata)
	return(obj)
}
#####Intial marker
GeneFilter<-function(obj){
	maxexp=dim(obj$newdata)[2]*0.8
	data(excludeGene)
	excludeGene=as.character(excludeGene[,1])
	geneSumm=obj$geneSumm
	index=match(geneSumm$gene,excludeGene)
	geneSumm=geneSumm[is.na(index),]
	data=obj$newdata
	index=match(geneSumm$gene,row.names(data))
	data=data[index,]
	binadata=obj$binadata
	binadata=binadata[index,]
	geneindex=rowSums(binadata)
	binadata=binadata[geneindex<maxexp,]
	data=data[geneindex<maxexp,]
	geneSumm=geneSumm[geneindex<maxexp,]
	obj$newdata=data
	obj$geneSumm=geneSumm
	obj$binadata=as.matrix(binadata)
	return(obj)
}
######Inial cluster
Bina<-function(data,cutoff=2){
	data[data<cutoff]=0
	data[data>=cutoff]=1
	return(data)
}
RankGene<-function(kk,k,HamD,geneName,n){
	x=HamD[kk,]
	xrank=rank(-x)
	MNNgene=geneName[xrank<=k&x>n]
	return(MNNgene)
}
MNNpair<-function(k,MNNgene,geneName){
	subgene=MNNgene[[k]]
	index=match(subgene,geneName)
	PP<-function(i,index,MNNgene,k,geneName){
		if (geneName[k] %in% MNNgene[[index[i]]])
			return(geneName[index[i]])
	}
	if(length(index[!is.na(index)])>0){
		a=do.call(cbind,lapply(1:length(index),PP,index=index,MNNgene=MNNgene,k=k,geneName=geneName))
		if (!is.null(a)){
			return(a)
		}
	}
}
getMNN<-function(HamD,genename,k,n){
	MNNgene=lapply(1:dim(HamD)[1],RankGene,k=k,HamD=HamD,geneName=genename,n=n)
	genePair=unique(as.character(do.call(cbind,lapply(1:length(MNNgene),MNNpair,MNNgene=MNNgene,geneName=genename))))
	return(genePair)
}
getMEN<-function(HamDD,genename,k,n){
	MNNgene=lapply(1:dim(HamDD)[1],RankGene,k=k,HamD=HamDD,geneName=genename,n=n)
	genePair=unique(as.character(do.call(cbind,lapply(1:length(MNNgene),MNNpair,MNNgene=MNNgene,geneName=genename))))
	return(genePair)
}
getMarker<-function(obj,k=300,n=30){
	data=obj$newdata
	binadata=obj$binadata
	genename=row.names(binadata)
	geneindex=rowSums(binadata)
	HamD=tcrossprod(binadata)
	diag(HamD)=0
	HamDD=tcrossprod((1-binadata),binadata)
	MNNmarker=getMNN(HamD=HamD,genename=genename,k=k,n=n)
	MENmarker=getMEN(HamD=HamDD,genename=genename,k=k,n=n)
	obj$marker=union(MNNmarker,MENmarker)
	return(obj)
}

######cluster
SCcluster<-function(obj){
	data=obj$rawdata
	gene=row.names(data)
	gene=unique(gene)
	index=match(gene,row.names(data))
	data=data[index,]
	marker=obj$marker
	pbmc <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200,
	                           project = "project")
	pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",
	                      scale.factor = 10000)
	pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
	                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
	pbmc <- ScaleData(object = pbmc)
	pbmc <- RunPCA(object = pbmc, pc.genes = marker, do.print = TRUE, pcs.print = 1:5,
	                genes.print = 5)
	pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:8,
	                      resolution = 0.6, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
	pbmc <- RunTSNE(object = pbmc, dims.use = 1:15, do.fast = TRUE)
	TSNE=pbmc@dr$tsne@cell.embeddings
	seuratCluster=pbmc@ident
	dbscanCluster=dbscan(TSNE,eps=1.2,minPts=15)$cluster
	names(dbscanCluster)=row.names(TSNE)
	obj$tsne=TSNE
	obj$seuratCluster=seuratCluster
	obj$dbscanCluster=dbscanCluster
	return(obj)
}



#####cluster-specific marker
getClusterGene<-function(obj,method){
	marker=obj$marker
	if (method=="Seurat"){
		seuratCluster=obj$seuratCluster
		gene=get_marker_genes(obj$rawdata[marker,names(seuratCluster)], seuratCluster)
		gene$gene=marker
		obj$clustergene=gene
		obj$method=method
	}else if (method == "dbscan"){
		dbscanCluster=obj$dbscanCluster
		gene=get_marker_genes(obj$rawdata[marker,names(dbscanCluster)], dbscanCluster)
		gene$gene=marker
		obj$clustergene=gene
		obj$method=method
	}
	return(obj)
}

HeatmapCluster<-function(obj,top,scale="none"){
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
HeatmapCell<-function(obj,top){
	feature=obj$clustergene
	topmarker=as.data.frame(feature %>% group_by(clusts) %>% top_n(top, auroc))
	topmarker=topmarker[order(topmarker$clusts),]
	method=obj$method
	if (method=="Seurat"){
		cluster=obj$seuratCluster
		cluster=cluster[order(cluster)]
	}else if (method=="dbscan"){
		cluster=obj$dbscanCluster
		cluster=cluster[order(cluster)]
	}
	clustertype=unique(cluster)
	data=obj$rawdata[topmarker$gene,names(cluster)]
	clustercount=table(cluster)
	gapIndex=c()
	for (i in 2:(length(clustertype)-1)){
		gapIndex=c(gapIndex,sum(clustercount[1:i]))
	}
	pheatmap(data,cluster_cols = F,cluster_rows = F,show_colnames=F,gaps_col=gapIndex)
}

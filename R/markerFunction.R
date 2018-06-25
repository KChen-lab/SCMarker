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
genepeak<-function(x,width){
	if (length(x[x!=0])>1){
		pp=density(x[x!=0],width=width)
		Hpeak=peak(pp)
		return(dim(Hpeak)[1])
	}else{
		return(0)
	}
}

#####filter data
ModalFilter<-function(data,geneK,cellK,width,cutoff){
	cellSumm=data.frame(cell=colnames(data))
	cellSumm$count<-sapply(1:dim(data)[2],genecount,data=data,index="col")
	data=data[,cellSumm$count>cellK]
	binadata=Bina(data,cutoff=2)
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
	filterres=list(data=data,geneSumm=geneSumm,cellSumm=cellSumm,binadata=binadata)
	return(filterres)
}
#####Intial marker
GeneFilter<-function(filterres,maxexp){
	data(excludeGene)
	excludeGene=as.character(excludeGene[,1])
	geneSumm=filterres$geneSumm
	index=match(geneSumm$gene,excludeGene)
	geneSumm=geneSumm[is.na(index),]
	data=filterres$data
	index=match(geneSumm$gene,row.names(data))
	data=data[index,]
	binadata=filterres$binadata
	binadata=binadata[index,]
	geneindex=rowSums(binadata)
	binadata=binadata[geneindex<maxexp,]
	data=data[geneindex<maxexp,]
	geneSumm=geneSumm[geneindex<maxexp,]
	filterres$data=data
	filterres$geneSumm=geneSumm
	filterres$binadata=as.matrix(binadata)
	return(filterres)
}
######Inial cluster
Bina<-function(data,cutoff){
	data[data<cutoff]=0
	data[data>=cutoff]=1
	return(data)
}
getMNN<-function(HamD,genename,MNN,MNNIndex){
	RankGene<-function(k,MNN,HamD,geneName){
		x=HamD[k,]
		xrank=rank(-x)
		MNNgene=geneName[xrank<=MNN&x>MNNIndex]
		return(MNNgene)
	}
	MNNgene=lapply(1:dim(HamD)[1],RankGene,MNN=MNN,HamD=HamD,geneName=genename)
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
	genePair=unique(as.character(do.call(cbind,lapply(1:length(MNNgene),MNNpair,MNNgene=MNNgene,geneName=genename))))
	return(genePair)
}
getMarker<-function(filterres,MNN,MNNIndex){
	data=filterres$data
	binadata=filterres$binadata
	genename=row.names(binadata)
	geneindex=rowSums(binadata)
	HamD=tcrossprod(binadata)
	diag(HamD)=0
	marker=getMNN(HamD=HamD,genename=genename,MNN=MNN,MNNIndex=MNNIndex)
	filterres$marker=marker
	return(filterres)
}



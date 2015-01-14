#Import an MAF file as a matrix
importMAF<-function(fn, col_loci="Hugo_Symbol",col_samples="Tumor_Sample_Barcode"){

    maf<-read.table(fn,sep="\t", header=TRUE, comment.char="")
    maf<-maf[,c(col_loci,col_samples)]

    mutationList<-as.character(unique(maf[,1]))
    sampleList<-as.character(unique(maf[,2]))

    mutMat<-matrix(0,length(mutationList),length(sampleList))
    rownames(mutMat)<-mutationList
    colnames(mutMat)<-sampleList

    for(i in 1:nrow(maf))
        mutMat[as.character(maf[i,1]),as.character(maf[i,2])]=1

    return(mutMat)
}


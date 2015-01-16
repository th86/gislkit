#Import data from the TCGA compressed package as a matrix
importTCGA<-function(path_customized="./", pattern_customized=".rsem.genes.normalized_results", verbose=TRUE){
    fileList=dir(path=path_customized,pattern=pattern_customized)
    example<-read.table(paste(path_customized, fileList[1], sep=""), header=TRUE, sep="\t"   )

    ge<-matrix(0,nrow(example),length(fileList))
    colnames(ge)<-fileList
    rownames(ge)<-as.character(example[,1])
    rm(example)
    for( fileName in fileList){
        e<-read.table(paste(path_customized, fileName, sep=""), header=TRUE, sep="\t"   )
        ge[,fileName]<-as.numeric(e[,2])
        if(verbose==TRUE)
            cat(which(fileList==fileName),"\n")
    }
    rm(e)
    return(ge)
}

#Converts UUIDs in a TCGA data matrix to TCGA barcodes. 
convertUUID2barcode<-function(ge,path ="./"){
    map<-read.table(paste(path, "FILE_SAMPLE_MAP.txt",sep=""),sep="\t",header=T)
    rownames(map)<-map[,1]
    cn<-map[colnames(ge),2]
    colnames(ge)<-cn
    return(ge)
}
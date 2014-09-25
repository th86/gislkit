#Extract gene symbols
getGeneSymbols = function(innames){
    outnames = sapply(innames, function(x){
    
        if(regexpr("\\?", x) > 0){
            o = strsplit(x, "\\|")[[1]][2]
        }else{
            o = strsplit(x, "\\|")[[1]][1]
        }
        return (o)

    }
    )
}

#Create attractor metagene list
attractorList.create<-function(){
attractorList<-list(LYM=c("SASH3","CD53","NCKAP1L","LCP2","IL10RA","PTPRC","EVI2B","BIN2","WAS","HAVCR2"), CIN=c("TPX2","KIF4A","KIFC1","NCAPG","BUB1","NCAPH","CDCA5","KIF2C","PLK1","CENPA"),
                    MES=c("COL3A1","COL5A2","COL1A2","THBS2","COL5A1","VCAN","COL6A3","SPARC","AEBP1","FBN1"),
                    CIN=c("TPX2","KIF4A","KIFC1","NCAPG","BUB1","NCAPH","CDCA5","KIF2C","PLK1","CENPA"),
                    END=c("CDH5","ROBO4","CXorf36","CD34","CLEC14A","ARHGEF15","CD93","LDB2","ELTD1","MYCT1"),
                    FS=c("FGD3","SUSD3"),
                    AHSA2=c("AHSA2","LOC91316","PILRB","ZNF767","TTLL3","CCNL2","PABC1L","LENG8","CHKB-CPT1B","SEC31B"),
                    IFIT=c("IFIT3","MX1","OAS2","RSAD2","CMPK2","IFIT1","IFIT44L","IFI44","IFI6","OAS1"),
                    WDR38=c("WDR38","YSK4","ROPN1L","C1orf194","MORN5","WDR16","RSPH4A","FAM183A","ZMYND10","DNAI1"),
                    MHC=c("HLA-DPA1","HLA-DRA","HLA-DPB1","HLA-DRB1","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DQA1","HLA-DRB5","HLA-DQB1"),
                    GIMAP=c("GIMAP4","GIMAP7","GIMAP6","GIMAP5","GIMAP8","GIMAP1","GIMAP2","TMEM176B","TMEM176A","ZNF777"),
                    chr8q24.3=c("SHARPIN","HSF1","TIGD5","GPR172A","ZC3H3","EXOSC4","SCRIB","CYHR1","MAF1","PUF60") )
}

#Create metagenes
meta.create<-function( ge, attractorList  ){
    meta=matrix(NA, length(attractorList), ncol(ge))
    rownames(meta)=names(attractorList)
    colnames(meta)=colnames(ge)
    for(i in 1:length(attractorList)){
        attractor_members=intersect( rownames(ge), attractorList[[i]])
        meta[i,]=colSums( ge[ attractor_members ,] )/length(attractor_members)
    }

    return(meta)
}

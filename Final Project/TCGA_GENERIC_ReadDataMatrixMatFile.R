TCGA_GENERIC_ReadDataMatrixMatFile <- function(Filename) {
     DataList=readMat(Filename)
     
     MATdata=as.matrix(DataList$RawData)
     rownames(MATdata)=DataList[[2]]
     colnames(MATdata)=DataList[[3]]
     
     return(MATdata)
}
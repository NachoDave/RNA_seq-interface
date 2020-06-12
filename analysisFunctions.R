## 
library(DESeq2)
source("RNA_SeqSaveClass.R")

deseq2CntNrm <- function(ob, expSmpNm, designFacs, ddsNm){
  
  # takes in a rna_seq analysis object and the name of the experimental sample table
  # and creates the count matrix and colData to create a DESeq2 object and normalizes
  # designFacs contains the names of the design factors
  
  browser()
  
  colDt <- ob@ExpSmpTab[[expSmpNm]]
  # make count matrix
 for (dx in 1:nrow(colDt)){
   
   tGnCnt <- as.data.frame(ob@GeneCntTables[colDt[[dx, 1]] == ob@GeneMeta[[1]]])
   
   if (dx == 1){
     
     gnCntMat <- data.frame(x = tGnCnt[, 4])
     rownames(gnCntMat) <- tGnCnt[, 1]
     
   } else
   {
     gnCntMat$x <- tGnCnt[, 4]
     
   }
   colnames(gnCntMat)[dx] <- colDt[[dx, 1]]
   
 } 
  
  desFacFrm <- "~"
  # make design formula
  cnt <- 0
  for (dx in designFacs$Factors){
    print(dx)
    cnt <- cnt + 1 
    if (cnt > 1){
      desFacFrm <- paste0(desFacFrm, " + ")
    }
    
    desFacFrm <- paste0(desFacFrm, dx)
  }
  
  # Create dds object and normalize
  
  dd <- DESeqDataSetFromMatrix(countData = gnCntMat, colData = colDt, design =  as.formula(desFacFrm))
  dd <- estimateSizeFactors(dd)
  
  ob <- addDESeqDS(ob, dd, ddsNm, expSmpNm)
  
  browser()
  return(ob)
}
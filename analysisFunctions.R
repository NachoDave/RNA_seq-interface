## 
library(DESeq2)
source("RNA_SeqSaveClass.R")
source("nrmCntResults.R")

# Function to normalize counts using DESeq2 
deseq2CntNrm <- function(ob, expSmpNm, designFacs, ddsNm, rmCnt = 10){
  
  # takes in a rna_seq analysis object and the name of the experimental sample table
  # and creates the count matrix and colData to create a DESeq2 object and normalizes
  # designFacs contains the names of the design factors
  
  #browser()
  
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
  for (dx in rev(designFacs$Factors)){
    print(dx)
    cnt <- cnt + 1 
    if (cnt > 1){
      desFacFrm <- paste0(desFacFrm, " + ")
    }
    
    desFacFrm <- paste0(desFacFrm, dx)
  }
  
  # Remove genes with low counts
  
  # Create dds object and normalize
  
  dd <- DESeqDataSetFromMatrix(countData = gnCntMat, colData = colDt, design =  as.formula(desFacFrm))
  
  # Remove genes with low counts
  #browser()
  keep <- rowSums(counts(dd)) >= rmCnt
  dd <- dd[keep,]
  
  
  # Estimate size factors
  dd <- estimateSizeFactors(dd)
  
  # Create a nrmCnt save object and add results
  nrmCntOb <- new("NrmCntResults")
  nrmCntOb <- addDESeqNrmCnts(nrmCntOb, dd, expSmpNm, paste(designFacs$Factors, collapse = "+"))
  #browser()
  ob <- addNrmCntsDS(ob, nrmCntOb, ddsNm, expSmpNm)
  
  #browser()
  return(ob)
}

# Function to perform differential analysis on DESeq2 object

deseq2DA <- function(ob, nrmCnts){
  browser()
  ob@NrmCnts[[nrmCnts]]@dds <- DESeq(ob@NrmCnts[[nrmCnts]]@dds)
  return(ob)
}

# Make DESeq2 results

deseq2Res <- function(ob, nrmCnts, fac, cntrst){
  
  # ob - analysis object
  # nrmCnts - name of the normed counts
  # fac - the analysis factor
  # Cntrst - the control level
  browser()
  
  levs <- levels(ob@NrmCnts[[nrmCnts]]@dds@colData@listData[[fac]])
  cntrstDESeq <- c(fac, cntrst, levs[!levs == cntrst])
  
  ob@NrmCnts[[nrmCnts]]@dds <- results(ob@NrmCnts[[nrmCnts]]@dds, contrast = cntrstDESeq)
  
}
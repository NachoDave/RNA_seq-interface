## 
library(DESeq2)
source("RNA_SeqSaveClass.R")
source("nrmCntResults.R")
source("pathWaySave.R")
# Function to normalize counts using DESeq2 
deseq2CntNrm <- function(ob, expSmpNm, designFacs, ddsNm, rmCnt = 10, selCol = 2){
  
  # takes in a rna_seq analysis object and the name of the experimental sample table
  # and creates the count matrix and colData to create a DESeq2 object and normalizes
  # designFacs contains the names of the design factors
  
  #browser()
  
  colDt <- ob@ExpSmpTab[[expSmpNm]]
  # make count matrix
 for (dx in 1:nrow(colDt)){
   
   tGnCnt <- as.data.frame(ob@GeneCntTables[colDt[[dx, 1]] == ob@GeneMeta[[1]]])
   
   if (dx == 1){
     
     gnCntMat <- data.frame(x = tGnCnt[, selCol])
     rownames(gnCntMat) <- tGnCnt[, 1]
     
   } else
   {
     gnCntMat$x <- tGnCnt[, selCol]
     
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
# 
# deseq2DA <- function(ob, nrmCnts){
#   browser()
#   ob@NrmCnts[[nrmCnts]]@dds <- DESeq(ob@NrmCnts[[nrmCnts]]@dds)
#   return(ob)
# }

# Make DESeq2 results

deseq2DA <- function(ob, nrmCnts, cntrstRef, cntrstCnd){
  #browser()
  
  ob@NrmCnts[[nrmCnts]]@dds[[cntrstRef]] <- relevel(ob@NrmCnts[[nrmCnts]]@dds[[cntrstRef]], ref = cntrstCnd)
  ob@NrmCnts[[nrmCnts]]@dds <- DESeq(ob@NrmCnts[[nrmCnts]]@dds)
  return(ob)
}

# Make DESeq2 results

deseq2Res <- function(ob, nrmCnts, fac, cntrst, LFC = FALSE, lfcMethod = "apeglm"){
  
  # ob - analysis object
  # nrmCnts - name of the normed counts
  # fac - the analysis factor
  # Cntrst - the control level
  #browser()
  
  # When labelling the comparison conditions DESeq2 puts the control variable 2nd on the contrast and coef
  levs <- levels(ob@NrmCnts[[nrmCnts]]@dds@colData@listData[[fac]])
  cntrstDESeq <- c(fac, levs[!levs == cntrst], cntrst  ) # now using the coef to allow the type of LFC to be set to apeglm
  coefDESeq <- paste(fac, levs[!levs == cntrst], 'vs', cntrst,   sep = "_")
  
  tryCatch(
    {ob@NrmCnts[[nrmCnts]] <- addDESeqRes(ob@NrmCnts[[nrmCnts]], results(ob@NrmCnts[[nrmCnts]]@dds, contrast = cntrstDESeq)) 
    #{ob@NrmCnts[[nrmCnts]] <- addDESeqRes(ob@NrmCnts[[nrmCnts]], results(ob@NrmCnts[[nrmCnts]]@dds, coef = coefDESeq)) 
    if (LFC){
      
      if (lfcMethod =="normal") # use the contrast for normal LFC/ use the coef for the others
      {
        ob@NrmCnts[[nrmCnts]] <- addDESeqResLFC(ob@NrmCnts[[nrmCnts]], lfcShrink(ob@NrmCnts[[nrmCnts]]@dds, contrast = cntrstDESeq, type = lfcMethod))
      }
      else
      {
        ob@NrmCnts[[nrmCnts]] <- addDESeqResLFC(ob@NrmCnts[[nrmCnts]], lfcShrink(ob@NrmCnts[[nrmCnts]]@dds, coef = coefDESeq, type = lfcMethod))
        
      }
    }
    return(ob)
    },
    error = function(e){return(e)}
  )
  
  
}

# Use webGesalt for pathway analysis
runWebGestaltR <- function(ob, db, enrichMeth, pth, nm, gns, idType, output, gnTb){
  tryCatch(
    {
      #browser()

     for (dbdx in db){
       for (mthdx in enrichMeth){
        
        #browser()
          
         refSet = "genome"
         minN = 5
         maxN = 2000
         topThr = 25 # number of results to report
         

         if (mthdx == "ORA"){
           tGns <- as.vector(gns[, 1])
           topThr = 25
         } else {
           tGns <- gns
           topThr = 20
         }

         # "GO Molecular Function", "GO Cellular Component", "GO Biological Process", # GO 
         # "KEGG Pathway", "Panther Pathway", 
         # "Reactome Pathway"
         
         if (dbdx == "GO Molecular Function"){
           enDb<- "geneontology_Molecular_Function"
         }
         else if(dbdx == "GO Cellular Component"){
           enDb<- "geneontology_Cellular_Component"
         }
         else if (dbdx == "GO Biological Process"){
           enDb<- "geneontology_Biological_Process"
         }
         else if (dbdx == "KEGG Pathway"){
           enDb<- "pathway_KEGG"
         }
         else if (dbdx == "Panther Pathway"){
           enDb<- "pathway_Panther"
         }
         else if (dbdx == "Reactome Pathway"){
           enDb<- "pathway_Reactome"
         }
         
         prjNm <- paste(nm, enDb, mthdx, sep = "_")
          # Overlap tissue with cell =====================================================================================#
         #browser() 
         x <- WebGestaltR(enrichMethod=mthdx, 
                                    organism="hsapiens",
                                    enrichDatabase=enDb, 
                                    interestGene = tGns,
                                     interestGeneType=idType, # gene id type 
                                     referenceSet = refSet,
                                     outputDirectory="/tmp/", 
                                    projectName=prjNm,
                                    isOutput=output, 
                                    sigMethod='top',
                                     minNum=minN,
                                    topThr = topThr
          ) # overlapping tissue, biological process
        assign(prjNm, x)
        #browser()
        xx <- new("webGStalt_Analysis")
        xx <- newPthWyAnl(xx, x, gnTb, mthdx, enDb)
        
        ob <- addPthWyAnl(ob, xx, prjNm)
        
        
      } # loop 1
    } # loop 2

      return(ob)
    }, # trycatch
  error = function(e){
    #browser()
    e$message <- paste("WebGStalt error: ", e$message, "for", prjNm, ". Try rerunning and removing this analysis", sep = " ")
    return(e)}
  )
}

# Compare gene lists ======================================================================================================= #
cmpGnList <- function(lst1, lst2, setOp){
  browser()
  # get stuff from the lists
  
  if (setOp == "In both"){
    
    lst <- intersect(rownames(lst1), rownames(lst2))
    
    
    if (any(colnames(lst1) == "Gene Symbol")){

      lst <- data.frame(ID = lst, `Gene Symbol` = lst1[lst, "Gene Symbol"], 
                        `log2FoldChange_1` = lst1[lst, "log2FoldChange"], 
                        `pvalue_1` = lst1[lst, "pvalue"],
                        `padj_1` = lst1[lst, "padj"],
                        log2FoldChange_2 = lst2[lst, "log2FoldChange"],
                        pvalue_2 = lst2[lst, "pvalue"],
                        padj_2 = lst2[lst, "padj"],
                        log2FoldChange_Mean = (lst1[lst, "log2FoldChange"] + lst2[lst, "log2FoldChange"])/2)

    } else if (any(colnames(lst1) == "Gene.Symbol")){
      lst <- data.frame(ID = lst, `Gene Symbol` = lst1[lst, "Gene.Symbol"])
      
    }
    else {
      
      lst <- data.frame(ID = lst, 
                        `log2FoldChange_1` = lst1[lst, "log2FoldChange"], 
                        `pvalue_1` = lst1[lst, "pvalue"],
                        `padj_1` = lst1[lst, "padj"],
                        log2FoldChange_2 = lst2[lst, "log2FoldChange"],
                        pvalue_2 = lst2[lst, "pvalue"],
                        padj_2 = lst2[lst, "padj"],
                        log2FoldChange_Mean = (lst1[lst, "log2FoldChange"] + lst2[lst, "log2FoldChange"])/2)
      
    }

  } else if (setOp == "In list 1 not 2"){
    
    lst <-setdiff(rownames(lst1), rownames(lst2))
    
    if (any(colnames(lst1) == "Gene Symbol")){
      
      lst <- data.frame(ID = lst, `Gene Symbol` = lst1[lst, "Gene Symbol"], 
                        `log2FoldChange` = lst1[lst, "log2FoldChange"], 
                        `pvalue` = lst1[lst, "pvalue"],
                        `padj` = lst1[lst, "padj"]
          )
      
    }
    else {
      
      lst <- data.frame(ID = lst, 
                        `log2FoldChange` = lst1[lst, "log2FoldChange"], 
                        `pvalue` = lst1[lst, "pvalue"],
                        `padj` = lst1[lst, "padj"]
      )
      
    }
    
  } else  { # in list 2 not 1
    
    lst <-setdiff(rownames(lst2), rownames(lst1))
    
    
    if (any(colnames(lst1) == "Gene Symbol")){
      
      lst <- data.frame(ID = lst, `Gene Symbol` = lst2[lst, "Gene Symbol"], 
                        `log2FoldChange` = lst2[lst, "log2FoldChange"], 
                        `pvalue` = lst2[lst, "pvalue"],
                        `padj` = lst2[lst, "padj"]
      )
      
    }
    else {
      
      lst <- data.frame(ID = lst, 
                        `log2FoldChange` = lst1[lst, "log2FoldChange"], 
                        `pvalue` = lst2[lst, "pvalue"],
                        `padj` = lst2[lst, "padj"]
      )
      
    }
    
  }
  
  rownames(lst) <- lst$ID
  
  return(lst)
  
}


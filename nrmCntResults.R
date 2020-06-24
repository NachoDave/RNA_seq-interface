library(DESeq2)
library(dplyr)
# class to save the nrm results
setClass("NrmCntResults", slots = list(dds = "DESeqDataSet", ExpSampNm = "character", NrmCnts = "matrix", NrmMethod = "character", design = "character"))

# dds - a deseq2 dataset (this is optional)
# expSampNm - name of experimental sample table used to generate counts matrix
# nrmCnts - matrix of normalzied counts
# nrmMethod - method used to normalize the data


# Add DESeq2 data

setGeneric(name = "addDESeqNrmCnts", 
           def=function(object, ds,  smExpNm, desFac)
           {
             standardGeneric("addDESeqNrmCnts")
             
           }
)

setMethod(f="addDESeqNrmCnts", 
          signature = "NrmCntResults",
          definition=function(object, ds, smExpNm, desFac){
            # ds is the DESeq object
            # smExpNm refers to the sample experiment list used to generate the normed object
            # 
            #browser()
            object@dds <- ds
            object@NrmCnts <- counts(ds, normalized = TRUE)
            object@ExpSampNm <- smExpNm
            object@NrmMethod <- "DESeq2"
            object@design <- desFac
            return(object) 
          }
)


# Remove DESeq2 data




# Add other methods
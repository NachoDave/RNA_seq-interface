# setClass("student", slots=list(name="character", age="numeric", GPA="numeric"))
# setMethod("show",
#           "student",
#           function(object) {
#             cat(object@name, "\n")
#             cat(object@age, "years old\n")
#             cat("GPA:", object@GPA, "\n")
#           }
# )
# 
# 
# x <- new("student", age = c(1, 23, 45))
# show(x)
library(dplyr)

setClass("RNASeqAnalysis", slots = list(GeneCntTables="list", GeneMeta = "data.frame", factorsTab = "list", ExpSmpTab = "list", ExpSmpFactorsTab = "list",
                                        DESeqNrmCnts = "list", DESeqNrmCntsExpSmp = "list"))
# GeneCntTables contains the gene count tables for an experiment

# Methods

# Add new gene counts table
setGeneric(name="newGeneCnts",
           def=function(object, geneCntTab, geneMeta)
           {
             standardGeneric("newGeneCnts")
           }
)

setMethod(f="newGeneCnts", # fntion name
          signature="RNASeqAnalysis",
          definition=function(object, geneCntTab, geneMeta){
            
            # add the gene count data frameto list
            
            object@GeneCntTables[[length(object@GeneCntTables) + 1]] <- geneCntTab
            
            # Add the meta Data
            
            
            if (ncol(object@GeneMeta) == 0){
              geneMeta$ID <- 1  
              object@GeneMeta <- (geneMeta)
              colnames(object@GeneMeta) <- c( "Description", "FileName", "ID")
              
            } else{
              geneMeta$ID <- max(object@GeneMeta$ID) + 1
              object@GeneMeta[nrow(object@GeneMeta) + 1,] <- geneMeta  
            }

            return(object)
          }
)

# Remove gene counts table
setGeneric(name="rmGeneCnts",
           def=function(object, rwDx)
           {
             standardGeneric("rmGeneCnts")
           }
)

setMethod(f="rmGeneCnts",
          signature="RNASeqAnalysis",
          definition=function(object, rwDx){
            
            object@GeneCntTables <- object@GeneCntTables[!(1:length(object@GeneCntTables) %in% rwDx)] # remove data
            
            object@GeneMeta <- object@GeneMeta[!(1:nrow(object@GeneMeta)) %in% rwDx,] # remove meta data
            
            
            return(object)
          }
          )

# Add new factors table
setGeneric(name = "addFactorsTab", 
           def=function(object, facTab, dx)
           {
             standardGeneric("addFactorsTab")
             
           }
           )

setMethod(f="addFactorsTab", 
          signature = "RNASeqAnalysis",
          definition=function(object, facTab, dx){
          #browser()
            object@factorsTab[[dx]] <- facTab 
           
            return(object) 
          }
          )

# Remove factors table

setGeneric(name = "rmFactorsTab", 
           def=function(object, facTab, rwDx)
           {
             standardGeneric("rmFactorsTab")
             
           }
)

setMethod(f="rmFactorsTab", 
          signature = "RNASeqAnalysis",
          definition=function(object, facTab, rwDx){
            #browser()
            object@factorsTab <- object@factorsTab[!(1:length(object@factorsTab) %in% rwDx)] # remove data 
            
            return(object) 
          }
)


# Add Experimental sample table

setGeneric(name = "addExpSmpTab", 
           def=function(object, smpExpTab, nm, facDx)
           {
             standardGeneric("addExpSmpTab")
             
           }
)

setMethod(f="addExpSmpTab", 
          signature = "RNASeqAnalysis",
          definition=function(object, smpExpTab, nm, facDx){
            #browser()
            object@ExpSmpTab[[nm]] <- smpExpTab 
            object@ExpSmpFactorsTab[[nm]] <- facDx
            return(object) 
          }
)

# Remove Experimental sample table

setGeneric(name = "rmExpSmpTab", 
           def=function(object,  nm)
           {
             standardGeneric("rmExpSmpTab")
             
           }
)

setMethod(f="rmExpSmpTab", 
          signature = "RNASeqAnalysis",
          definition=function(object, nm){
            #browser()
            object@ExpSmpTab[[nm]] <- NULL 
            object@ExpSmpFactorsTab[[nm]] <- NULL
            return(object) 
          }
)

# Add deseq data set
setGeneric(name = "addDESeqDS", 
           def=function(object, ds, nm, smExpNm)
           {
             standardGeneric("addDESeqDS")
             
           }
)

setMethod(f="addDESeqDS", 
          signature = "RNASeqAnalysis",
          definition=function(object, ds, nm, smExpNm){
            #browser()
            object@DESeqNrmCnts[[nm]] <- ds
            object@DESeqNrmCntsExpSmp[[nm]] <- smExpNm
            return(object) 
          }
)

# Remove deseq2 data

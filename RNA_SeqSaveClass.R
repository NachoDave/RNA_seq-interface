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

setClass("RNASeqAnalysis", slots = list(GeneCntTables="list", GeneMeta = "data.frame"))
# GeneCntTables contains the gene count tables for an experiment

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
              
              object@GeneMeta <- (geneMeta)
              colnames(object@GeneMeta) <- c( "Description", "FileName")
              
            } else{
              object@GeneMeta[nrow(object@GeneMeta) + 1,] <- geneMeta  
            }

            return(object)
          }
)

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




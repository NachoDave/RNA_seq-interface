# Create a class for storing webGstalt results

setClass("webGStalt_Analysis", slots = c(pthWyAnl = "data.frame", geneTable = "character", 
                                         method = "character", db = "character"))

setGeneric(name = "newPthWyAnl", {
  
  def = function(object, pwa, gnTb, mth, db){
    
    standardGeneric("newPthWyAnl")
    
  }
  

})

setMethod(f = "newPthWyAnl", 
          signature = "webGStalt_Analysis",
          definition = function(object, pwa, gnTb, mth, db){

            object@pthWyAnl <- pwa
            object@geneTable <- gnTb
            object@method <- mth
            object@db <- db
            
            return(object)
            
          })
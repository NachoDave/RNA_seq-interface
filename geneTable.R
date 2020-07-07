
# set a class to store gene tables from rna seq app

setClass("geneTableDA", slots = c(gnTbl = "data.frame", nrmCnts = "character", filters = "character"))

# Genetable - gene table
# nrmCnts - name of normed counts object used to create table
# filter - filters used on table

# Method to add a egencount table
setGeneric(name = "addGeneTableDA", {
  
  def = function(object, gnTb, nrmCntNm, fil){
    
    standardGeneric("addGeneTableDA")
    
  }
    
  
})

setMethod(f = "addGeneTableDA", 
          signature = "geneTableDA",
          definition = function(object, gnTb, nrmCntNm, fil){
            
            object@gnTbl <- gnTb
            object@nrmCnts <- nrmCntNm
            object@filters <- fil
            
            return(object)
            
          })

library(ggplot2)
library(plotly)
library("pheatmap")
library("RColorBrewer")

maplot <- function(DEseqRes, tit = "", xlims = NULL, ylims = NULL){
  #browser()
  df <- as.data.frame(DEseqRes)
  df$`p < 0.05` <- df$padj < 0.05
  df$`p < 0.05`[is.na(df$`p < 0.05`)] = FALSE
  plt <- ggplot(df,aes(x = baseMean, y = log2FoldChange, colour = `p < 0.05`, text = paste("Gene:", gsub("\\.\\d+", "", rownames(df))))) +
    geom_point(size = 0.35)  +  ylab("log fold change") + xlab("Mean of normalized Counts") + 
    #theme(legend.position = "none") + s
    geom_hline(yintercept = 0, linetype = "dashed") + 
    #browser()
    ggtitle(paste0(gsub("log2\ fold\ change\ \\(.+\\):\ ", "", DEseqRes[[1]]@elementMetadata$description[[2]]), tit))  + 
    scale_color_manual(values = c("red", "blue"), name = "P < 0.05")
  
  # set axis limits if given
    if (!is.null(xlims)){
    
      plt <- plt + scale_x_continuous(trans='log10', limits = 10^xlims) 
      
    }
  else
  {
    plt <- plt + scale_x_continuous(trans='log10')
  } 
    if (!is.null(ylims)){
    
    plt <- plt + ylim(low = ylims[1], high = ylims[2])  
  }


  # 
  fig <- ggplotly(plt, tooltip = c("text", "baseMean", "log2FoldChange"), source = "A") %>% config(displaylogo = FALSE,
                                                                                     modeBarButtonsToRemove = list(
                                                                                       'sendDataToCloud',
                                                                                       'pan2d',
                                                                                       'autoScale2d',
                                                                                       #'resetScale2d',
                                                                                       'hoverClosestCartesian',
                                                                                       'hoverCompareCartesian', 
                                                                                       'select2d',
                                                                                       'lasso2d',
                                                                                       'drawline',
                                                                                       'toggleSpikelines ',
                                                                                       'zoomIn2d',
                                                                                       'zoomOut2d',
                                                                                       'toggleSpikelines'
                                                                                     )) %>% event_register("plotly_click")
  
  return(fig) 
  
}


volplot <- function(DEseqRes, tit = "", xlims = NULL, ylims = NULL){
  #browser()
  df <- as.data.frame(DEseqRes)
  # Remove padj == na
  df <- df[!is.na(df$padj) ,]
  
  # Get colours up/down
  df$upDwn <- "Not Sig"
  df[df$padj < 0.05 & df$log2FoldChange > 0,]$upDwn <- "Up"
  df[df$padj < 0.05 & df$log2FoldChange < 0,]$upDwn <- "Down"
  
  # df$`p < 0.05` <- df$padj < 0.05
  # df$`p < 0.05`[is.na(df$`p < 0.05`)] = FALSE
  df$log10P <- -log10(df$padj)
  plt <- ggplot(df,aes(x = log2FoldChange, y = log10P, colour = upDwn , text = paste("Gene:", gsub("\\.\\d+", "", rownames(df))))) +
    geom_point(size = 0.35) + ylab("-log10(P adjusted)") + xlab("log fold change") + 
    #theme(legend.position = "none") + #geom_hline(yintercept = 0, linetype = "dashed") + 
    ggtitle(paste0(gsub("log2\ fold\ change\ \\(.+\\):\ ", "", DEseqRes[[1]]@elementMetadata$description[[2]]), tit)) + 
    scale_color_manual(values = c("blue", "gray", "red"), name = "Up/down Regulated") 
  
  if (!is.null(xlims)){
    
    plt <- plt + xlim(low = xlims[1], high = xlims[2])  
    
  }
  
  if (!is.null(ylims)){
    
    plt <- plt + ylim(low = ylims[1], high = ylims[2])  
  }
  
  
  
  
  # 
  fig <- ggplotly(plt, tooltip = c("text", "baseMean", "log2FoldChange"),  source = "B") %>% config(displaylogo = FALSE,
                                                                                     modeBarButtonsToRemove = list(
                                                                                       'sendDataToCloud',
                                                                                       'pan2d',
                                                                                       'autoScale2d',
                                                                                       #'resetScale2d',
                                                                                       'hoverClosestCartesian',
                                                                                       'hoverCompareCartesian', 
                                                                                       'select2d',
                                                                                       'lasso2d',
                                                                                       'drawline',
                                                                                       'toggleSpikelines ',
                                                                                       'zoomIn2d',
                                                                                       'zoomOut2d',
                                                                                       'toggleSpikelines'
                                                                                     )) %>% event_register("plotly_click")
  
  return(fig)
  
}

plotPthWybar <- function(pthWy){
  
  #browser()
  df <- pthWy@pthWyAnl
  #df <- df[order(df$enrichmentRatio, decreasing = T), ]
  df$p0p05 <- "FDR > 0.05"
  df$p0p05[df$FDR <= 0.05] <- "FDR < 0.05"
  p <- ggplot(df, aes(x = description, y = enrichmentRatio,  text = paste("Enrichment:", enrichmentRatio, "P Val:", pValue,"Genes:", userId))) +
    geom_col(colour = "black", aes(fill = p0p05), width = 0.7 ) + coord_flip() + 
    labs(fill = "", x = "", y = "Enrichment Ratio" , title = paste("Gene Table:", pthWy@geneTable), subtitle = paste("DB:", pthWy@db, "Method:", pthWy@method))
             #ggtitle(paste("Gene Tab:", pthWy@geneTable, "\nDB:", pthWy@db, "Method:", pthWy@method))
  
  fig <- ggplotly(p, tooltip = c( "text"),  source = "C") %>% config(displaylogo = FALSE,
                                                                                                    modeBarButtonsToRemove = list(
                                                                                                      'sendDataToCloud',
                                                                                                      'pan2d',
                                                                                                      'autoScale2d',
                                                                                                      #'resetScale2d',
                                                                                                      'hoverClosestCartesian',
                                                                                                      'hoverCompareCartesian', 
                                                                                                      'select2d',
                                                                                                      'lasso2d',
                                                                                                      'drawline',
                                                                                                      'toggleSpikelines ',
                                                                                                      'zoomIn2d',
                                                                                                      'zoomOut2d',
                                                                                                      'toggleSpikelines'
                                                                                                    )) %>% event_register("plotly_click") %>%
    layout(title = list(text = paste0(paste("Gene Tab:", pthWy@geneTable),
                                      '<br>',
                                      '<sup>',
                                      paste("DB:", pthWy@db, "Method:", pthWy@method),
                                      '</sup>')))
  return(fig)
}


# Sample distance plot
pltSmpDist <- function(cntMat, trnsFrm){
  
  if (trnsFrm == "vst"){
    trnDt <- vst(cntMat@dds, blind = F)
  } 
  else if (trnsFrm == "rlog"){
    trnDt <- rlog(cntMat@dds, blind = F)
    
  }
  else {
    trnDt <- cntMat@dds
    
  }

  dst <- dist(t(assay(trnDt)))
    
  sampleDistMatrix <- as.matrix(dst)
  rownames(sampleDistMatrix) <- colnames(cntMat@dds)
  colnames(sampleDistMatrix) <- colnames(cntMat@dds)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  fig <- pheatmap(sampleDistMatrix,
           clustering_distance_rows=dst,
           clustering_distance_cols=dst,
           col=colors)
  
  return(fig)
  
}

pltHtMp <- function(cntMat, trnsFrm){
  
  if (trnsFrm == "vst"){
    trnDt <- vst(cntMat@dds, blind = F)
  } 
  else if (trnsFrm == "rlog"){
    trnDt <- rlog(cntMat@dds, blind = F)
    
  }
  else {
    trnDt <- normTransform(cntMat@dds)
    
  }
  
  #browser()
  select <- order(rowMeans(counts(cntMat@dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(cntMat@dds)[,names(colData(cntMat@dds))[3: length(names(colData(cntMat@dds)))- 1]])
  fig <- pheatmap(assay(trnDt)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  
  return(fig)
}

pltPCA <- function(cntMat, trnsFrm){
  
  if (trnsFrm == "vst"){
    trnDt <- vst(cntMat@dds, blind = F)
  } 
  else if (trnsFrm == "rlog"){
    trnDt <- rlog(cntMat@dds, blind = F)
    
  }
  else {
    trnDt <- normTransform(cntMat@dds)
    
  }
  #browser()
  pcaDt <- plotPCA(trnDt, intgroup=c(names(colData(cntMat@dds))[3: length(names(colData(cntMat@dds)))- 1]), returnData=TRUE)
  p <- ggplot(pcaDt, aes(x = PC1, y = PC2, color = name)) +
    geom_point(size=3) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) 
  # percentVar <- round(100 * attr(pcaData, "percentVar"))
  # ggplot(pcaData, aes(PC1, PC2, color=as.formula(names(colData(cntMat@dds))[2]), shape=as.formula(names(colData(cntMat@dds))[3]))) +
  #   geom_point(size=3) +
  #   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #   coord_fixed()
  fig <- ggplotly(p) %>% config(displaylogo = FALSE,
                                                                     modeBarButtonsToRemove = list(
                                                                       'sendDataToCloud',
                                                                       'pan2d',
                                                                       'autoScale2d',
                                                                       #'resetScale2d',
                                                                       'hoverClosestCartesian',
                                                                       'hoverCompareCartesian', 
                                                                       'select2d',
                                                                       'lasso2d',
                                                                       'drawline',
                                                                       'toggleSpikelines ',
                                                                       'zoomIn2d',
                                                                       'zoomOut2d',
                                                                       'toggleSpikelines'
                                                                     )) 
  
  return(fig)
}

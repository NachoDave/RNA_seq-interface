
library(ggplot2)
library(plotly)
library("pheatmap")
library("RColorBrewer")
library(ggforce)
library(heatmaply)

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
    geom_point(size = 0.25) + ylab("-log10(P adjusted)") + xlab("log fold change") + 
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
  
  if (is.null(df$enrichmentRatio)){
    
    p <- ggplot(df, aes(x = description, y = enrichmentScore,  text = paste("Enrichment:", sprintf("%.4g", enrichmentScore), "P Val:", sprintf("%.4g", pValue),
                                                                            "FDR:", sprintf("%.4g", FDR)))) +
      geom_col(colour = "black", aes(fill = p0p05), width = 0.7 ) + coord_flip() + 
      labs(fill = "", x = "", y = "Enrichment Score" , title = paste("Gene Table:", pthWy@geneTable), subtitle = paste("DB:", pthWy@db, "Method:", pthWy@method))
    
  } else
  {
  p <- ggplot(df, aes(x = description, y = enrichmentRatio,  text = paste("Enrichment:", sprintf("%.4g", enrichmentRatio), "P Val:", sprintf("%.4g", pValue),
                                                                          "FDR:", sprintf("%.4g", FDR)))) +
    geom_col(colour = "black", aes(fill = p0p05), width = 0.7 ) + coord_flip() + 
    labs(fill = "", x = "", y = "Enrichment Ratio" , title = paste("Gene Table:", pthWy@geneTable), subtitle = paste("DB:", pthWy@db, "Method:", pthWy@method))
             #ggtitle(paste("Gene Tab:", pthWy@geneTable, "\nDB:", pthWy@db, "Method:", pthWy@method))
}
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
 # browser()
 fig <- heatmaply(sampleDistMatrix, colors = rev(Blues(100)), grid_color = "black") %>% config(displaylogo = FALSE,
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
  # colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  # fig <- pheatmap(sampleDistMatrix,
  #          clustering_distance_rows=dst,
  #          clustering_distance_cols=dst,
  #          col=colors)
 #browser()

  
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
  fig <- heatmaply(assay(trnDt)[select,], Rowv = FALSE, Colv = FALSE, grid_color = "black", colors = heat.colors(100))  %>% config(displaylogo = FALSE,
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

pltPCA <- function(cntMat, trnsFrm, color = "None", shape = "None"){
  
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
  
  if (color == "None" & shape == "None"){
    
    color = name
    
  }
  
  if (color == "None"){color = NULL}
  if (shape == "None"){shape = NULL}
  
  pcaDt <- plotPCA(trnDt, intgroup=c(names(colData(cntMat@dds))[3: length(names(colData(cntMat@dds)))- 1]), returnData=TRUE)
  percentVar <- round(100 * attr(pcaDt, "percentVar"))
  p <- ggplot(pcaDt, aes_string(x = "PC1", y = "PC2", color = color, shape = shape)) +
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

# Venn diagram
pltVenn <- function(labelA, labelB, txt){
  df.venn <- data.frame(x = c(-.866, 0.866), y = c(0, 0), labels = c(labelA, labelB))

  
  x <- seq(-1.5, 1.5, 0.01)
  xN <- length(x)
  y <- sqrt(1.5^2 - x^2)
  
  # positions <- data.frame(
  #   id = c(rep(labelA, 4), rep(labelB, 8)),
  #   x = c(x, rev(x)c(2, 1, 1.1, 2.2, 1, 0, 1.1, 3.2)),
  #   y = c(c(-0.5, 0, 1, 0.5, -0.5, 0, 1, 0.5))
  # )
  
positions <- data.frame(
  id = c(rep(labelA, 2*xN), rep(labelB, 2*xN)),
  x = c(c(x, rev(x)) - 0.866, c(x, rev(x)) + 0.866),
  y = c(c(y, rev(-y)), c(y, rev(-y))) 
)
  
  
  # Currently we need to manually merge the two together
  #datapoly <- merge(values, positions, by = c("id"))
  
  venn <- ggplot(positions, aes(x = x, y = y, fill = factor(id))) +
    geom_polygon(aes( group = id), alpha = 0.5, colour = 'gold') + 
    coord_fixed() + 
    theme_void() + 
    labs(fill = NULL) +
    theme(legend.position = 'bottom') + 
    annotate("text", x = c(-1, 0, 1), y = c(0, 0, 0), label = txt, size = 5)
  

  # venn <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) + 
  #          geom_circle(alpha = 0.4, size = 0.5, colour = 'grey') + 
  #          coord_fixed() + 
  #          theme_void() + 
  #         theme(legend.position = 'bottom') + 
  #         scale_fill_manual(values = c('gold', 'firebrick')) +
  #         scale_colour_manual(values = c('gold', 'firebrick'), guide = FALSE) + 
  #   labs(fill = NULL) + 
  #   annotate("text", x = c(-1, 0, 1), y = c(0, 0, 0), label = txt, size = 5)
  
  a <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  
  fig = ggplotly(venn) %>% config(displaylogo = FALSE,
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
                                  )) %>% layout(xaxis = a, yaxis = a)
  return(fig)
  #return(venn)
}
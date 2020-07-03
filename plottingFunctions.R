
library(ggplot2)
library(plotly)

maplot <- function(DEseqRes, tit = "", xlims = NULL, ylims = NULL){
  browser()
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
  fig <- ggplotly(plt, tooltip = c("text", "baseMean", "log2FoldChange")) %>% config(displaylogo = FALSE,
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
    
    plt + xlim(low = xlims[1], high = xlims[2])  
    
  }
  
  if (!is.null(ylims)){
    
    plt + ylim(low = ylims[1], high = ylims[2])  
  }
  
  
  
  
  # 
  fig <- ggplotly(plt, tooltip = c("text", "baseMean", "log2FoldChange")) %>% config(displaylogo = FALSE,
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
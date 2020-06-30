
library(ggplot2)
library(plotly)

maplot <- function(DEseqRes){
  #browser()
  df <- as.data.frame(DEseqRes)
  df$`p < 0.05` <- df$padj < 0.05
  df$`p < 0.05`[is.na(df$`p < 0.05`)] = FALSE
  plt <- ggplot(df,aes(x = baseMean, y = log2FoldChange, colour = `p < 0.05`, text = paste("Gene:", gsub("\\.\\d+", "", rownames(df))))) +
    geom_point(size = 0.35) + scale_x_continuous(trans='log10') + ylim(-4, 4) + ylab("log fold change") + xlab("Mean of normalized Counts") + 
    theme(legend.position = "none") + geom_hline(yintercept = 0, linetype = "dashed") + 
    ggtitle(gsub("log2\ fold\ change\ \\(MLE\\):\ ", "", DEseqRes[[1]]@elementMetadata$description[[2]]))  + 
    scale_color_manual(values = c("red", "blue"))
  # 
  fig <- ggplotly(plt, tooltip = c("text", "baseMean", "log2FoldChange"))
  
  return(fig) 
  
}


volplot <- function(DEseqRes){
  browser()
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
    ggtitle(gsub("log2\ fold\ change\ \\(MLE\\):\ ", "", DEseqRes[[1]]@elementMetadata$description[[2]])) + 
    scale_color_manual(values = c("blue", "gray", "red"))
  # 
  fig <- ggplotly(plt, tooltip = c("text", "baseMean", "log2FoldChange"))
  
  
  return(fig)
  
}
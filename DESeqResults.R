library(DESeq2)

# class to save the results of DESeq2 analysis
setClass("DESeqResults", slots = list(dds = "DESeqDataSet", expSampNm = "character", diffAn = "logical", cutoff = "integer"))

# dds - a deseq2 dataset
# expSampNm - name of the experimental sample table used to generate the object
# vsd - vst normalized dataset
# rlog - rlog normalized dataset
# diffAn - logical to indicate if differential analysis has been performed
# cutoff - keep tracks with minimum number of reads
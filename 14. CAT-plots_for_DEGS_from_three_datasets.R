## Install package
if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("matchBox")

library("matchBox")

## Prepare the DEGs for the three datasets. In the following, we merged the three datasets by "logFC", and ordered the genes by decreasing "logFC", meaning up-regulated genes are at the top of the list.

silicosis <- as.data.frame(read.table("Silicosis_DESeq2_results.adjp0.1.FC1.5.pro-coding.xls",header=T,sep="\t"))
IPF=as.data.frame(read.table("IPFVsNorm_genes.FC1.5_adjp0.1.pro-coding.xls",header=T,sep="\t"))
COPD=as.data.frame(read.table("COPD_GSE76925_limma_output.adjp0.1.FC1.5.pro-coding.xls",header=T,sep="\t"))
silicosis$"SYMBOL"=as.character(silicosis$"SYMBOL")
IPF$"SYMBOL"=as.character(IPF$"SYMBOL")
COPD$"SYMBOL"=as.character(COPD$"SYMBOL")

merge <- list(silicosis=silicosis, IPF=IPF, COPD=COPD)

data <- mergeData(merge, idCol="SYMBOL", byCol="logFC")

##Computing CAT curves without a reference ranking

catHigh2LowNoRefByEqualRanks <- computeCat(data = data, idCol = 1, method="equalRank", decreasing=TRUE) ## up-regulated genes are at the top of the list

### Computing Probability Intervals

 PIbyRefEqualRanks03 <- calcHypPI(data=data, expectedProp=0.3) 

### Plotting the results
pdf("catHigh2LowNoRefByEqualRanks03Prop.pdf",width=10)

plotCat(catData = catHigh2LowNoRefByEqualRanks,
preComputedPI=PIbyRefEqualRanks03,
cex=1.2, lwd=1.2, cexPts=1.2, spacePts=30,
main="CAT curves for decreasing logFC",
where="center", legend=TRUE, legCex=1, ncol=1,
plotLayout = layout(matrix(1:2, ncol = 2), widths = c(2,1)))
dev.off()


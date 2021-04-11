## The DEGs of silicosis mouse lungs were downloaded from previously published paper (in its Table S1). The info of the paper as follows:
 
# Paper title: Comparative RNA-Seq transcriptome analysis on silica induced pulmonary inflammation and fibrosis in mice silicosis model
#DOI: 10.1002/jat.3587

# Before this analysis, we transfered the mouse genes to their corresponding human homologs using the information from following website  (http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt). 

library("matchBox")
human <- as.data.frame(read.table("Silicosis_DESeq2_results.adjp0.1.FC1.5.pro-coding.xls",header=T,sep="\t"))
mouse=as.data.frame(read.table("mouse_diff.renamed.txt",header=T,sep="\t"))
human$"SYMBOL"=as.character(human$"SYMBOL")
mouse$"SYMBOL"=as.character(mouse$"SYMBOL")

merge <- list(human=human, mouse=mouse)


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


############## CAT-plot for the metacore pathways of both human and mouse silicosis ##############
## Before this analysis, we performed pathway analysis using Metacore online (https://portal.genego.com).We downloaded the results and retained those pathways with FDR < 0.05 as significant ones. We ordered these pathways with increasing "FDR", meaning the most significant pathways are at the top of the list, and then used these pathways as input.

library("matchBox")
Silicosis <- as.data.frame(read.table("silicosis_diff_all_2605_metacore_FDR0.05.txt",header=T,sep="\t"))
mouse=as.data.frame(read.table("mouse_silicosis_DEGs_metacore_FDR0.05.txt",header=T,sep="\t"))

Silicosis$"Pathway"=as.character(Silicosis$"Pathway")
mouse$"Pathway"=as.character(mouse$"Pathway")


merge <- list(Silicosis=Silicosis, mouse=mouse)

data <- mergeData(merge, idCol="Pathway", byCol="FDR")

##Computing CAT curves without a reference ranking

catLow2HighNoRefByEqualRanks <- computeCat(data = data, idCol = 1, method="equalRank", decreasing=FALSE) 

### Computing Probability Intervals

 PIbyRefEqualRanks03 <- calcHypPI(data=data, expectedProp=0.3) 

### Plotting the results

pdf("catLow2HighNoRefByEqualRanks03Prop.pathway.pdf",width=10)
plotCat(catData = catLow2HighNoRefByEqualRanks,
preComputedPI=PIbyRefEqualRanks03,
cex=1.2, lwd=1.2, cexPts=1.2, spacePts=30,
main="CAT curves for increasing FDR",
where="center", legend=TRUE, legCex=1, ncol=1,
plotLayout = layout(matrix(1:2, ncol = 2), widths = c(2,1)))
dev.off()
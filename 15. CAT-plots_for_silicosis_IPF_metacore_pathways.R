## Load the significant metacore pathways for the two datasets
# Before this analysis, we performed pathway analysis using Metacore online (https://portal.genego.com).We downloaded the results and retained those pathways with FDR < 0.05 as significant ones. We ordered these pathways with increasing "FDR", meaning the most significant pathways are at the top of the list, and then used these pathways as input.

library("matchBox")
Silicosis <- as.data.frame(read.table("silicosis_diff_all_2605_metacore_FDR0.05.txt",header=T,sep="\t"))
IPF=as.data.frame(read.table("IPF_diff_all_3609_metacore_FDR0.05.txt",header=T,sep="\t"))

Silicosis$"Pathway"=as.character(Silicosis$"Pathway")
IPF$"Pathway"=as.character(IPF$"Pathway")

merge <- list(Silicosis=Silicosis, IPF=IPF)

data <- mergeData(merge, idCol="Pathway", byCol="FDR")

##Computing CAT curves without a reference ranking

catLow2HighNoRefByEqualRanks <- computeCat(data = data, idCol = 1, method="equalRank", decreasing=FALSE) ## the most significant pathways are at the top of the list


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


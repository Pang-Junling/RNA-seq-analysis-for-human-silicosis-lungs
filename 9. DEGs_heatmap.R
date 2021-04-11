## Install pheatmap package
install.packages(pheatmap)
library(pheatmap)

## Prepare the top 20 up- / down- regulated genes 
# We order the DEGs based on the adjusted p-value, and selected the top 20 up- / down- regulated genes. Prepared their FPKMs in the file "DEGs.top20.FPKM.xls".

DEG=as.matrix(read.table("DEGs.top20.FPKM.xls",header=T,row.names=1,sep="\t"))
group=as.data.frame(read.table("sample_group.txt",header=T,row.names=1))

pdf("DEGs.top20.pheatmap.pdf")

ann_colors = list(group=c(Donor="#1E90FF",Silicosis="#FF6A6A"))

pheatmap(DEG,scale="row",clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",clustering_method="ward.D",cutree_cols=2,cutree_rows=2, annotation_col = group,annotation_colors = ann_colors, border=TRUE, border_color = "grey60")

dev.off()

## sample_group.txt:
#Sample  group
#e1      Donor
#e2      Donor
#e4      Donor
#e5A     Donor
#e16     Donor
#e17     Donor
#e18     Donor
#e8      Silicosis
#e9      Silicosis
#e11     Silicosis
#e12     Silicosis
#e19     Silicosis
#e21     Silicosis
#e23     Silicosis
#e24     Silicosis
#e29     Silicosis
#e30     Silicosis
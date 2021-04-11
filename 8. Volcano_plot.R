## Install "ggplot2"
options(CRAN="https://?cloud.r-project.org/")
install.packages("ggplot2")

## head DESeq2_results.xls
#Gene    baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
#ENSRNOG00000000001      35.9339244252972        -0.035585049626701      0.0683339672503142      -0.520751992875657      0.602539549832823       NA
#ENSRNOG00000000007      1.65963758350419        0.01231228673879        0.018972167829563       0.648965729662407       0.516360527573299       NA


file="DESeq2_results.xls"

library(ggplot2)
data=read.table(file=file,header=T,row.names=1)

############### three colors ##############
data$threshold <- as.factor(ifelse(data$padj < 0.1 & abs(data$log2FoldChange) >= log2(1.5),ifelse(data$log2FoldChange > log2(1.5) ,'Up','Down'),'Not'))  ### three colors
###################

pdf("volcano_plot.pdf")

ggplot(data=data,
            aes(x=log2FoldChange, y =-log10(padj),
                colour=threshold, fill=threshold)) +
  scale_color_manual(values=c("blue","gray","red"))+
  geom_point(alpha=0.4, size=1.75) +
  xlab("log2 fold change") + ylab("-log10 adjusted p") +
  theme_bw()+
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.1), lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right")

dev.off()


## Install "limma" packge
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")

## load "limma" and perform differential gene analysis
library("limma")
data=as.matrix(read.table("GSE76925_series_matrix.txt",header=T,row.names=1))
design=as.matrix(read.table("GSE76925_design.txt",header=T,row.names=1))
fit=lmFit(data,design)
contrasts.matrix<-makeContrasts(case_control = case - control,levels=design)
huvec_fit<-contrasts.fit(fit,contrasts.matrix)
huvec_ebFit <- eBayes(huvec_fit)
result=topTable(huvec_ebFit,number=Inf,coef = 1)
diff=na.omit(result)
write.table(diff,file="GSE76925_limma_output.txt",sep="\t",quote=F)


library("DESeq2")

## Import the counts data 
header=read.table("read-counts.xls",header=F,sep="\t",nrows=1)
header=t(header[-1])
count=read.table("read-counts.xls",header=T,sep="\t",row.names=1)
colnames(count)=header
count=as.matrix(count)

## DEGs
coldata=read.table("SampleInfo.xls",header=T,sep="\t") # coldata, the first column list samples, and the following columns record the "condition" or "type" for each sample
dds=DESeqDataSetFromMatrix(countData=count,colData=coldata,design=~batch+type)
dds=DESeq(dds)   
res=results(dds)
write.table(res,"DESeq2_results.xls",quote=F,sep='\t')


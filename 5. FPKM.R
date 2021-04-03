library("DESeq2")

## Import the counts data 
header=read.table("read-counts.xls",header=F,sep="\t",nrows=1)
header=t(header[-1])
count=read.table("read-counts.xls",header=T,sep="\t",row.names=1)
colnames(count)=header
count=as.matrix(count)

## Normalize counts
coldata=read.table("SampleInfo.xls",header=T,row.names=1,sep="\t") # coldata, the first column list samples, and the following columns record the "condition" or "type" for each sample
dds=DESeqDataSetFromMatrix(countData=count,colData=coldata,design=~batch+type)  # design=~1 when the condition of all samples are the same
my.count.raw=counts(dds,normalized=F)
dds=estimateSizeFactors(dds)
my.counts.normalized=counts(dds,normalized=T)

## Import gene infomation
library("GenomicFeatures")
txdb=makeTxDbFromGFF("/data1/pangjl/ref/human/Homo_sapiens.GRCh38.84.gtf",format="gtf")
genes_list=exonsBy(txdb,by="gene")
#txdb=makeTxDbFromGFF("lncRNA.final.new.gtf",format="gtf")
#genes_list=exonsBy(txdb,by="tx",use.names=TRUE)

## Judge whether the genes_list (from gtf file) and dds (DESeqDataSet) have the same genes and the same order
if (length(rownames(dds)) == sum(rownames(dds) == names(genes_list))) {print("genes_list and dds have the same gene order")}
rowRanges(dds) <- genes_list

## Calculate the fpkm of genes
my.fpkm.from_normalized_counts <- fpkm(dds,robust = TRUE)   ### between-samples normalization
my.fpkm.from_raw_counts <- fpkm(dds,robust=FALSE)  ### within-sample normalization

## Write the results to files
write.table(my.fpkm.from_normalized_counts,file="fpkm-from-norm-counts.xls",quote=F,sep='\t',col.names=NA)
write.table(my.fpkm.from_raw_counts,file="fpkm-from-raw-counts.xls",quote=F,sep='\t',col.names=NA)


## Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
BiocManager::install("GenomicAlignments")

## Set the path of files
files.path="./bam_files/"  ### directory with bam files
file.name="read-counts.xls"   ## output file

## Read the gene annotation file 
library("GenomicFeatures")
txdb=makeTxDbFromGFF("/data/wangj/Ref/human_ensembl/grch38_snp_tran/Homo_sapiens.GRCh38.84.gtf",format="gtf")
genes_list=exonsBy(txdb,by="gene")

## count reads 
library("GenomicAlignments")
fls <- list.files(path=files.path,pattern="*.bam$")
bamlst <- BamFileList(paste(files.path,fls,sep=""),index=character(),yieldSize=100000,asMates = TRUE)
se=summarizeOverlaps(features=genes_list,reads=bamlst,mode="Union",singleEnd=FALSE,ignore.strand=TRUE,fragments=TRUE)
write.table(assays(se)$counts,file=paste(files.path,file.name,sep=""),quote=F,sep='\t',col.names = NA)



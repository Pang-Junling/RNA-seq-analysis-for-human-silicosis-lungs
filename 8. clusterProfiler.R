library(clusterProfiler)
library(org.Hs.eg.db) ## Mm, Hs, Rn

#GeneSymbol, GeneID, Ensembl,UniprotID, ncbi-proteinid...

gene_list= "DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.list" ## no header
BP_out_file="DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.BP.xls"
CC_out_file="DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.CC.xls"
MF_out_file="DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.MF.xls"
KEGG_out_file="DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.KEGG.xls"


BP_plot="DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.BP.top20.pdf"
CC_plot="DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.CC.top20.pdf"
MF_plot="DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.MF.top20.pdf"
KEGG_plot="DESeq2_results.adjp0.1.FC1.5.pro-coding.name.down.KEGG.top20.pdf"



gene=read.table(gene_list,header=F)
gene=as.vector(gene[,1])
gene_1=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")  ### change
gene=gene_1$ENTREZID

BP <- enrichGO(OrgDb="org.Hs.eg.db", gene = gene, ont = "BP", pvalueCutoff = 0.5, readable= TRUE) # pAdjustMethod = "BH", qvalueCutoff = 0.05,
write.table(as.data.frame(BP@result),file=BP_out_file,sep="\t",quote=F)


CC <- enrichGO(OrgDb="org.Hs.eg.db", gene = gene, ont = "CC", pvalueCutoff = 0.5, readable= TRUE) # pAdjustMethod = "BH", qvalueCutoff = 0.05,
write.table(as.data.frame(CC@result),file=CC_out_file,sep="\t",quote=F)


MF <- enrichGO(OrgDb="org.Hs.eg.db", gene = gene, ont = "MF", pvalueCutoff = 0.5, readable= TRUE) # pAdjustMethod = "BH", qvalueCutoff = 0.05,
write.table(as.data.frame(MF@result),file=MF_out_file,sep="\t",quote=F)


ekk <- enrichKEGG(gene= gene, organism  = 'hsa', pvalueCutoff = 0.5)  # qvalueCutoff = 0.1
write.table(as.data.frame(ekk@result),file=KEGG_out_file, sep="\t",quote=F)

pdf(BP_plot)
dotplot(BP,showCategory=20,title="Enrichment GO Top20")
barplot(BP, showCategory=20,title="EnrichmentGO") 
dev.off()

pdf(CC_plot)
dotplot(CC,showCategory=20,title="Enrichment GO Top20")
barplot(CC, showCategory=20,title="EnrichmentGO") 
dev.off()

pdf(MF_plot)
dotplot(MF,showCategory=20,title="Enrichment GO Top20")
barplot(MF, showCategory=20,title="EnrichmentGO") 
dev.off()

pdf(KEGG_plot)
dotplot(ekk,showCategory=20,title="Enrichment KEGG Top20")
barplot(ekk, showCategory=20,title="EnrichmentKEGG") 
dev.off()



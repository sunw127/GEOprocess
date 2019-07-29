rm(list = ls())
###### gpl annotation
library(GEOquery)
gpl <- getGEO('GPL96', destdir=".")
colnames(Table(gpl))
probe2gene <- Table(gpl)[,c(1,11,13)]
head(probe2gene)
save(probe2gene,gpl,file = "probe2gene_3.Rdata")
ids <- probe2gene
head(ids)
colnames(ids)=c('probe_id','symbol','transcript_id')
ids=ids[ids$symbol != '',]
dim(ids)
load("GSE19143_matrix_2.Rdata")
ids=ids[ids$probe_id %in% rownames(exprSet),]
exprSet[1:4,1:4]
exprSet <- exprSet[ids$probe_id,]
dim(exprSet)
ids$median <- apply(exprSet,1,median)  ####ids create an new col
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
exprSet <- exprSet[ids$probe_id,]
rownames(exprSet) <- ids$symbol
exprSet[1:4,1:4]
save(ids,exprSet,file = "gene_annotation_4.Rdata")

###### Genecode annotation
rm(list = ls())
load("gene_annotation_4.Rdata")
library(rtracklayer)
gene_annotation <- rtracklayer::import("gencode.v31.annotation.gtf.gz")
gene_annotation <- as.data.frame(gene_annotation)
gene_annotation <- cbind(gene_annotation$gene_type,gene_annotation$gene_name)
colnames(gene_annotation) <- c("gene_type","gene_name")
gene_annotation <- as.data.frame(gene_annotation)
gene_annotation <- gene_annotation[!duplicated(gene_annotation$gene_name),]
rownames(gene_annotation) <- gene_annotation[,2]
gene_annotation_matrix <- gene_annotation[ids$symbol,]
genetype_list <- gene_annotation_matrix[,1]
table(genetype_list)
lncRNA_matrix <- gene_annotation_matrix[grep("lncRNA",genetype_list),]
mRNA_matrix <- gene_annotation_matrix[grep("protein_coding",genetype_list),]
lncRNA_matrix <- exprSet[rownames(exprSet) %in% rownames(lncRNA_matrix),]
mRNA_matrix <- exprSet[rownames(exprSet) %in% rownames(mRNA_matrix),]
dim(lncRNA_matrix)
dim(mRNA_matrix)
save(lncRNA_matrix,mRNA_matrix,file = "lncRNA_mRNA_matrix_5.Rdata")

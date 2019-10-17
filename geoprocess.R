### step1 for download
rm(list = ls())
library(GEOquery)
gset <- getGEO('GSE42387', destdir=".",
               AnnotGPL = F,     
               getGPL = F)       
b = gset[[1]]
exprSet=exprs(b)  
head(exprSet)
#exprSet <- log2(exprSet)
pdata = pData(b)
colnames(pdata) 
length(colnames(pdata)) 
pdata[,37]   
group_list = as.character(pdata[ ,37])
design = model.matrix(~0+factor(group_list))
rownames(design) = colnames(exprSet)
colnames(design) = c("control","oxaliplatin","irinotecan")
design
dim(exprSet) 

{
  irinotecan_sensitive_expr = exprSet[,rownames(design[design[,1] == 1,])]
  dim(irinotecan_sensitive_expr)
  irinotecan_resistant_expr = exprSet[,rownames(design[design[,3] == 1,])]
  dim(irinotecan_resistant_expr)
  exprSetIrinotecan = cbind(irinotecan_sensitive_expr, irinotecan_resistant_expr)  
  dim(exprSetIrinotecan)
}
save(exprSet,exprSetIrinotecan,design,file = "download_1.Rdata")

### step2 for normalization
rm(list = ls())
load("download_1.Rdata")
library(ggplot2)
library(reshape2)
exprSet0 = exprSetIrinotecan
data_m <- melt(exprSet0)
head(data_m)
colnames(data_m) = c("Probe","Sample","Value")
p <- ggplot(data_m, aes(x=Sample, y=Value),color = Sample) + 
  geom_boxplot(aes(fill=factor(Sample))) + 
  theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
  theme(legend.position="none")
p
ggsave(filename = "bn.png")
library(preprocessCore)
exprSet1 = normalize.quantiles(exprSet0) #### 归一化处理芯片
rownames(exprSet1) = rownames(exprSet0)
colnames(exprSet1) = colnames(exprSet0)
data_m <- melt(exprSet1)
head(data_m)
colnames(data_m) = c("Probe","Sample","Value")
p <- ggplot(data_m, aes(x=Sample, y=Value),color = Sample) + 
  geom_boxplot(aes(fill=factor(Sample))) + 
  theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
  theme(legend.position="none")
p
ggsave("an.png")
exprSetNormalization = exprSet1
save(exprSetNormalization,file = "Normalization_2.Rdata")

### step3 for annotation
rm(list = ls())
library(GEOquery)
gpl <- getGEO('GPL16297', destdir=".")
colnames(Table(gpl))
probe2gene <- Table(gpl)[,c(1,3,4,5)]
head(probe2gene)
ids <- probe2gene
head(ids)
ids = na.omit(ids)
dim(ids)
colnames(ids)=c('probe_id','entrez_id','symbol','name')
ids=ids[ids$symbol != '',]
dim(ids)
load("Normalization_2.Rdata")
ids = ids[ids$probe_id %in% rownames(exprSetNormalization),]
exprSetNormalization <- exprSetNormalization[ids$probe_id,]
dim(exprSetNormalization)
dup = as.matrix(table(duplicated(ids$symbol)))#查看是否有重复
if(dup[2,] > 0){
  ids$median <- apply(exprSetNormalization,1,median)  #中位值判断
  ids = ids[order(ids$symbol,ids$median,decreasing = T),]
  ids = ids[!duplicated(ids$symbol),]
  ids = na.omit(ids)
}
exprSetNormalization = exprSetNormalization[ids$probe_id,]
dim(exprSetNormalization)
ids = ids[match(rownames(exprSetNormalization),ids$probe_id),]
rownames(exprSetNormalization) = ids[,3]
save(exprSetNormalization,ids,file = "annotation_3.Rdata")

### step4 for DE
rm(list = ls())
load("annotation_3.Rdata")
head(exprSetNormalization)
sampleTree = hclust(dist(t(exprSetNormalization)), method = "average")
pdf(file = "hclust.pdf")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
group_list = as.character(c(rep("sensitive",9),rep("resistance",9)))
design = model.matrix(~0+factor(group_list))
rownames(design) = colnames(exprSetNormalization)
colnames(design) = c("resistant","sensitive")
design
library(limma)
fit <- lmFit(exprSetNormalization,design)
contrast.matrix <- makeContrasts("resistant-sensitive",levels = design)
contrast.matrix
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
head(nrDEG)
save(nrDEG,file = "DE_4.Rdata")
rm(list = ls())
load("DE_4.Rdata")
#for volcano
if(T){
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df = nrDEG
  df$v= -log10(P.Value) 
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  df$g=ifelse(df$P.Value>0.05,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
              ifelse( df$logFC >1,'up', 
                      ifelse( df$logFC < -1,'down','stable') )
  )
  table(df$g)   #down40/stable19147/up40
  df$name = rownames(df)
  head(df)
  write.table(rownames(df[df$g == 'up',]),
              file = 'up_gene.txt',
              quote=F,row.names=F,col.names=F)
  write.table(rownames(df[df$g == 'down',]),
              file = 'down_gene.txt',
              quote=F,row.names=F,col.names=F)
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g',palette = c("#0066CC", "#616161", "#FF3300"),ylab = "log10(pValue)",
            title = "Irinotecan-resistant colon cancer Valcano Plot")
  ggsave('volcano.png')
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("#0066CC","#FF3300","#616161") )
  ggsave('MA.png')
}
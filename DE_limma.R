rm(list = ls())
load("GSE19143_eSet_1.Rdata")

##### for DE
b = gset[[1]]  ## 降级提取b
exprSet=exprs(b)  ## 获取表达矩阵，发现已是log处理后的数据。通常表达矩阵的原始数字从0但好几百万都不等，需要进行归一化处理。14488875 elements
exprSet <- log2(exprSet)
pdata = pData(b)
colnames(pdata)
group_list=as.character(pdata[,39])
table(group_list)
design <- model.matrix(~0+factor(group_list))
rownames(design) <- rownames(pdata)
colnames(design)= c("prednisolone_resistant","prednisolone_sensitive")
design
library(limma)
fit <- lmFit(exprSet,design)
contrast.matrix <- makeContrasts("prednisolone_resistant-prednisolone_sensitive",levels = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
head(nrDEG)
write.csv(nrDEG,"DEG_all.csv")
save(nrDEG,fit2,file = "GSE19143_DEG_7.Rdata")

##### for DElncRNA
rm(list = ls())
load("lncRNA_mRNA_matrix_5.Rdata")
dim(lncRNA_matrix)
load("GSE19143_eSet_1.Rdata")
b = gset[[1]] 
pdata = pData(b)
group_list=as.character(pdata[,39])
table(group_list)
design <- model.matrix(~0+factor(group_list))
rownames(design) <- rownames(pdata)
colnames(design)= c("prednisolone_resistant","prednisolone_sensitive")
design
library(limma)
fit <- lmFit(lncRNA_matrix,design)
contrast.matrix <- makeContrasts("prednisolone_resistant-prednisolone_sensitive",levels = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG_lncRNA = na.omit(tempOutput)
save(nrDEG_lncRNA,fit2,file = "nrDEG_lncRNA_8.Rdata")

##### for DEmRNA
rm(list = ls())
load("lncRNA_mRNA_matrix_5.Rdata")
dim(mRNA_matrix)
load("GSE19143_eSet_1.Rdata")
b = gset[[1]] 
pdata = pData(b)
group_list=as.character(pdata[,39])
table(group_list)
design <- model.matrix(~0+factor(group_list))
rownames(design) <- rownames(pdata)
colnames(design)= c("prednisolone_resistant","prednisolone_sensitive")
design
library(limma)
fit <- lmFit(mRNA_matrix,design)
contrast.matrix <- makeContrasts("prednisolone_resistant-prednisolone_sensitive",levels = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG_mRNA = na.omit(tempOutput)
save(nrDEG_mRNA,fit2,file = "nrDEG_mRNA_9.Rdata")

## for volcano 
if(T){
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
              ifelse( df$logFC >1,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                      ifelse( df$logFC < -1,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
  )
  table(df$g)
  df$name=rownames(df)
  head(df)
  write.table(rownames(df[df$g == 'up',]),
              file = 'up_gene.txt',
              quote=F,row.names=F,col.names=F)
  write.table(rownames(df[df$g == 'down',]),
              file = 'down_gene.txt',
              quote=F,row.names=F,col.names=F)
  
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "name", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = c('TTC9', 'AQP3', 'CXCL11','PTGS2'), #挑选一些基因在图中显示出来
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  ggsave('volcano.png')
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  ggsave('MA.png')
  
  
}


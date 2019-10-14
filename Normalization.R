library(ggplot2)
library(reshape2)
data_m <- melt(exprSet)
head(data_m)
colnames(data_m) = c("Probe","Sample","Value")
p <- ggplot(data_m, aes(x=Sample, y=Value),color = Sample) + 
  geom_boxplot(aes(fill=factor(Sample))) + 
  theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
  theme(legend.position="none")
p
ggsave(filename = "bn.png")
library(preprocessCore)
exprSet1 = normalize.quantiles(exprSet) #### 归一化处理芯片
rownames(exprSet1) = rownames(exprSet)
colnames(exprSet1) = colnames(exprSet)
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
save(exprSetNormalization,group_list,file = "Combined_1.Rdata")
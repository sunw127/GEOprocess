rm(list = ls())
library(GEOquery)
gset <- getGEO('GSE19143', destdir=".",
               AnnotGPL = T,     
               getGPL = T)       
save(gset,file='GSE19143_eSet_1.Rdata')   ## 保存到本地
load('GSE19143_eSet_1.Rdata')  ## 载入数据

b = gset[[1]]  ## 降级提取b
exprSet=exprs(b)  ## 获取表达矩阵，发现已是log处理后的数据。通常表达矩阵的原始数字从0但好几百万都不等，需要进行归一化处理。14488875 elements
exprSet <- log2(exprSet)
pdata = pData(b)  ## 使用函数?pData获取样本临床信息（如性别、年龄、肿瘤分期等等）
colnames(pdata) ## 查看列名，即查看所包含的临床信息类型，
length(colnames(pdata)) ## 查看pdata列数 
pdata[,39]   ## 查看第67列triple negative status的数据情况,发现按照"Tumor"，"not Tumor"排列
group_list = as.character(pdata[, 39])  ## 将67列改成字符
dim(exprSet) ## 查看矩阵的维度 
table(group_list) 


{
  prednisolone_sensitive_expr = exprSet[, grep("prednisolone sensitiv", group_list)]# 提出"not Tumor"所有行数据 3663225 elements
  dim(prednisolone_sensitive_expr)
  prednisolone_resistant_expr = exprSet[, grep("prednisolone resistant", group_list)]  # 提出Tumor数据 10825650
  dim(prednisolone_resistant_expr)
  exprSet = cbind(Imatinib_naive_expr, Imatinib_resistant_expr)  ## 14488875 elements
  dim(exprSet)
  ## Normal_expr：矩阵，行为探针ID，列为non_Tumor样本
  ## Tumor_expr：矩阵，行为探针ID，筛选出列为Tumor的样本
  ## exprSet变成了Normal....Tumor....排序的矩阵，那么后面的group_list也应该相应排序，才能对上信息
}

write.csv(exprSet,"GSE19143_exprSet.csv")
save(b,group_list,gset,prednisolone_sensitive_expr,prednisolone_resistant_expr,exprSet,pdata,file = "GSE19143_matrix_2.Rdata")

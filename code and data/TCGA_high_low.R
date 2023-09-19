###

###
library(tidyverse)
library(readxl)
library(DESeq2)

### 
data <- read.csv("",row.names = 1)
data <- data[rowMeans(data)>1,]              # 去除表达量过低的基因
data <- as.data.frame(t(data))

kras_expression <- data$TFAP4

# 计算KRAS基因的平均值
kras_mean <- mean(kras_expression)

# 根据平均值区分高表达组和低表达组

high_expression_group <- as.data.frame(t(data[kras_expression >= kras_mean, ]))
write.csv(high_expression_group,file = 'high_expression_group.csv')
low_expression_group <- as.data.frame(t(data[kras_expression < kras_mean, ]))
write.csv(low_expression_group,file = 'low_expression_group.csv')

countdata <- read.csv('',row.names = 1,header = T)
condition <- factor(c(rep("High",183),rep("Low",290)))
colData <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = colData, design = ~ condition)
head(dds)


dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
 #将结果用result()函数来获取
res <- results(dds1)

# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1$log2FoldChange <- -res1$log2FoldChange
write.csv(res1,file = 'CDKN1A_high_low.csv')

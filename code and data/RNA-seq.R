### RNA-seq analysis ----
### author honghuizhao

###load the packages for analysis-----
setwd("C:/Users/ZHH/Desktop/data")
library(rjson)
library(readxl)

### 获取样本信息 -----------
json <- jsonlite::fromJSON("metadata.cart.2023-10-14.json")
#id <- json$associated_entities[[1]][,1]
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  


### 获取gdc_download文件夹下的所有TSV表达文件的 路径+文件名 -----------
count_file <- list.files("",pattern = '*.tsv',recursive = TRUE)
#在count_file中分割出文件名
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})

### 创建数据表达矩阵 ---------------
matrix = data.frame(matrix(nrow=60660,ncol=0))
for (i in 1:length(count_file)){
  path = paste0('',count_file[i])
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  data <- data[3]   #取出unstranded列（第3列），即count数据，对应其它数据
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

write.csv(matrix,'COUNT_matrix.csv',row.names = TRUE)

### 设置Gene Symbol为列名的矩阵（前面得到的是Ensembl ID）------------------------------------------
path = paste0('',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
matrix0 <- cbind(gene_name,matrix)
#将gene_name列去除重复的基因，保留每个基因最大表达量结果
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)    
#将gene_name列设为行名
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-1]
write.csv(matrix0,file = 'matrix.csv')
### 分为normal和tumor矩阵--------------------------
sample <- colnames(matrix0)

normal <- c()
tumor <- c()

for (i in 1:length(sample)){
  if((substring(colnames(matrix0)[i],14,15)>10)){    #14、15位置大于10的为normal样本
    normal <- append(normal,sample[i])
  } else {
    tumor <- append(tumor,sample[i])
  }
}

tumor_matrix <- matrix0[,tumor]
normal_matrix <- matrix0[,normal]
write.csv(normal_matrix,file = 'normal_lung_matrix.csv')
write.csv(tumor_matrix,file = 'tumor_lung_matrix.csv')


### 临床数据整合 ---------------------
json <- jsonlite::fromJSON("metadata.cart.2023-08-04.json")
entity_submitter_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
sample_case <- t(rbind(entity_submitter_id,case_id))

clinical <- read.delim('clinical/clinical.tsv',header = T)
clinical <- as.data.frame(clinical[duplicated(clinical$case_id),])

clinical_matrix <- merge(sample_case,clinical,by="case_id",all.x=T)
clinical_matrix <- clinical_matrix[,-1]
write.csv(clinical_matrix,file = 'clinical_matrix.csv')

### different expression gene analysis ----------
library(DESeq2)
library(tidyverse)


lung <- read.csv('lung.csv',row.names = 1,header = T)
### 比较开始 -----------
countData <- lung[rowMeans(lung)>1,]              # 去除表达量过低的基因

condition <- factor(c(rep("Tumor",18),rep("Normal",18)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
head(dds)


dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
#将结果用result()函数来获取
res <- results(dds1)

# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

write.csv(res1,file = 'deg.csv')

count <- read.csv('deg.csv')
### RNA-seq base analysis was end ---------
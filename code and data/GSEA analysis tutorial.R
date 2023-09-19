###  

# GSEA analysis tutorial --- R language

# author Honghui Zhao 2023/5/21

###

#GS：通路描述信息
#GS DETAILS：每一个通路富集信息汇总，含上面的条形码图和表格
#SIZE：该通路中包含表达数据集文中的基因数目（经过条件筛选后的值）
#ES：富集分数
#NES：标准化后的富集分数
#NOM p-val：是对 ES 的统计学分析，用来表征富集结果的可信度
#FDR q-val：是多重假设检验校正之后的 p-value，即对NES可能存在的假阳性结果的概率估计。GSEA 对显著性的定义为 p-value<5%，FDR q-val<25%
#FWER p-val：家庭错误率；也就是说，NES代表假阳性发现的更保守估计的概率
#RANK AT MAX：当 ES 最大时，对应基因所在排序好的基因列表中所处的位置
#LEADING EDGE：tags 表示核心基因占该通路基因集的百分比；list 表示核心基因占所有基因的百分比；signal，将前 2 项统计值结合在一起计算出的富集信号强度

#人类基因集:9 个
#人类分子特征数据库 (MSigDB) 中的 33196 个基因集分为 9 个主要集合和几个子集合。

#H: hallmark gene sets
#该类别包含了由多个已知的基因集构成的超基因集，每个 H 类别的基因集都对应多个基础的其他类别的基因集。比如 HALLMARK_APOPTOSIS 对应 80 个基因集。

#C1: positional gene sets
#与人类染色体细胞遗传学条带相对应的位置基因组。

#C2: curated gene sets
#来自在线通路数据库、PubMed 出版物和领域专家知识的精选基因集。

#C3: motif gene sets
#基于对 microRNA 种子序列和预测的转录因子结合位点的基因靶点预测的调控靶基因组。

#C4: computational gene sets
#通过挖掘大量面向癌症的微阵列数据定义的计算基因集。

#C5: GO gene sets
#由相同GO术语注释的基因组成。比如 GO_LYSOSOME（溶酶体）对应 GO:0005764。

#C6: oncogenic signatures gene sets
#致癌特征基因组：直接从来自癌症基因扰动的微阵列基因表达数据中定义。

#C7: immunologic signatures gene sets
#免疫特征基因组：代表免疫系统内的细胞状态和扰动。

#C8: cell type signatures gene sets
#细胞类型特征基因组：从人体组织单细胞测序研究中确定的簇标记中收集。

#小鼠基因集:6 个
#MH，是小鼠通过同源基因与人类 H 基因集建立的具有映射关系的数据库。同理 M1、M2、M3、M5、M8 的功能对应人类数据集的 C1、C2、C3、C5、C8"

### Set the path


### load the R package for the analysis 
library(msigdbr)
library(clusterProfiler)
library(stringr)
library(tidyverse)
library(export)
library(readxl)
library(org.Hs.eg.db)
library(forcats)
library(ggstance)
library(GseaVis)
library(gggsea)
library(Hmisc)

### load the data for analysis
diff <- read.csv("",header = T)
colnames(diff)[1] <- 'symbol'
### clusterProfiler --------------------

### function of GSEA -----------------------
### Extract KEGG database information, which can be customized, and check the official website of GSEA for details
### creating the geneList for analysis
geneList <- diff$log2FoldChange
names(geneList) = diff$symbol
geneList = sort(geneList, decreasing = TRUE)

KEGG_df <- msigdbr(species = "Homo sapiens")

KEGG_df <- KEGG_df %>% 
                   dplyr::select(gs_name, gene_symbol)

set.seed(111) 
KEGG.res <- GSEA(geneList,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.15,
                 pAdjustMethod = "BH",
                 verbose = FALSE,
                 eps = 0,
                 TERM2GENE = KEGG_df)

GSEA.KEGG.result <- KEGG.res@result

### save the data of GSEA_KEGG
write.csv(GSEA.KEGG.result ,file = "GSEA.csv")
up <- read.csv("")
### GSEA Visualization of results
mygene <- up$genename
gseaNb(object = KEGG.res,
       geneSetID ='KEGG_ECM_RECEPTOR_INTERACTION',
       addPval = T,
       addGene = mygene,
       subPlot = 2)

### bacth plot
terms <- c('KEGG_VEGF_SIGNALING_PATHWAY',
           'KEGG_MAPK_SIGNALING_PATHWAY',
           'KEGG_NON_SMALL_CELL_LUNG_CANCER')

gseaNb(object = KEGG.res,
         geneSetID = terms,
         addPval = T,
         addGene = mygene)

# add segment line
dotplotGsea(data = KEGG_ges_result,topn = 10,
            order.by = 'NES',
            add.seg = T)
ggsave('GSEA-KEGG-UP-DOWN.png',width = 15,height = 12)
### GSEA_KEGG ----------------------
### creating the genelist for analysis
genesymbol <- diff$symbol
entrezID <- bitr(genesymbol,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")# turn genesymbol into enetezID

genelist <- diff$log2FoldChange
names(genelist) <- diff$symbol
genelist <- genelist[names(genelist) %in% entrezID[,1]]
names(genelist) <- entrezID[match(names(genelist),entrezID[,1]),2]
genelist <- sort(genelist,decreasing = T)

### GSEA_KEGG analysis by clusterProfiler
KEGG_gesa <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.15,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

### Extract the results and convert the IDs in the results into genesymbols
KEGG_ges_result = setReadable(KEGG_gesa,
                              OrgDb = "org.Hs.eg.db",
                              keyType = "ENTREZID")
result <- KEGG_ges_result@result
write.csv(result,file = 'GSEA-KEGG-TCGA.csv')

KEGG_gesa@result<-KEGG_gesa[order(KEGG_gesa$NES,decreasing=T)]
cnetplot(KEGG_ges_result,showCategory = c("ECM-receptor interaction","PI3K-Akt signaling pathway","Focal adhesion","TNF signaling pathway","VEGF signaling pathway"))
ggsave()

### Save the data
write.csv(GO_gesa_result,file = "GSEA_GO-TCGA.csv")

### GSEA_KEGG Visualization of results
# add segment line
dotplotGsea(data = KEGG_gesa,topn = 10,
            order.by = 'NES',
            add.seg = T)

gseaNb(object = KEGG_gesa,
       geneSetID ='hsa04010',
       addPval = T)

### GSEA_GO analysis by clusterProfiler
GO_gesa <- gseGO(
  geneList     = genelist,
  OrgDb        = org.Hs.eg.db,
  ont          = "ALL",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.15,
  verbose      = FALSE,
  seed = FALSE, by = "fgsea")

GO_gesa = setReadable(GO_gesa,
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENTREZID")


GO_gesa_result <- GO_gesa@result

GO_gesa@result<-GO_gesa[order(GO_gesa$NES,decreasing=T)]
# add segment line
dotplotGsea(data = GO_gesa,topn = 10,
            order.by = 'NES',
            add.seg = T)
ggsave('GSEA-GO-UP-DOWN.png',width = 15,height = 8)
gseaNb(object = GO_gesa,
       geneSetID ='GO:0001227',
       addPval = T)

### save the data
write.table(GO_ges_result ,file="goGSEA.txt",sep="  ",quote=F,row.names = F)
##geneList是前面整理好的排序后的ENTREZID
##OrgDb选择物种注释的数据库
##ont同enrichGO分析，可选3大分类BP、CC、MF或ALL
##nPerm置换检验的次数，默认为1000
##minGSSize富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤掉，默认为10
##maxGSSize富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤掉，默认为500
##verbose是否输出提示信息
##seed是否使结果具有可重复性
##by选择使用的统计学方法，默认为fgsea


### fgsea------------------------------------------
library(fgsea)
library(data.table)

geneset_KEGG = msigdbr(species = "Homo sapiens") %>% 
                       dplyr::select(gs_name,gene_symbol)
geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)# remove KEGG_
geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)
geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)

geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)
GSEA_geneset_KEGG <- geneset_KEGG %>% 
                     split(x = .$gene_symbol, f = .$gs_name)

df <- diff[order(diff$avg_log2FC,decreasing = T),]
ranks <- df$avg_log2FC
names(ranks) <- df$gene

## analysis of fgsea
GSEA_df <- fgsea(pathways = GSEA_geneset_KEGG, 
                 stats = ranks,
                 minSize=10,
                 maxSize=500,
                 eps=0.0)

### fgsea_KEGG Visualization of results
topPathwaysUp <- GSEA_df[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- GSEA_df[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(GSEA_geneset_KEGG[topPathways],ranks, GSEA_df,gseaParam = 0.5)
ggsave('fgesa-top10-up-down.png',width = 10,height = 10,path = path1)
graph2ppt(file = 'GSEA-table.pptx',height = 7,width = 8.5)

save.image('GSEA.RData')
fwrite(GSEA_df, file="fGSEA.csv")


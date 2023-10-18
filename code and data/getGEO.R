library(GEOquery)
library(umap)

# load series and platform data from GEO
gset <- getGEO('GSE45827', destdir = "C:/Users/ZHH/Desktop/GSE", AnnotGPL = F, getGPL = F)
gset = gset[[1]]
exprset = exprs(gset)

####----002_Expression_log2_transform----####
ex <- exprset
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  exprset <- log2(ex)
  print("log2 transform finished")}else{
    print("log2 transform not needed")}

library(limma) 
boxplot(exprset,outline=FALSE, notch=T, las=2)
exprset=normalizeBetweenArrays(exprset)
boxplot(exprset,outline=FALSE, notch=T, las=2)
exprset = as.data.frame(exprset)

####----003_ID_transform_not_in_platfomMap----####
gpl <- getGEO('GPL570',destdir ="C:/Users/ZHH/Desktop/GSE")
colnames(Table(gpl))
ids <- Table(gpl)[,c(1,11)]
exprset$ID <- rownames(exprset)

exprSet <- merge(x = ids, y = exprset, by = "ID")
library(dplyr)
exprSet=rename(exprSet,GeneSymbol ="Gene Symbol")
exprSet=exprSet[-grep('///',exprSet$GeneSymbol),]   
exprSet=avereps(exprSet[,-c(1,2)],         
                ID = exprSet$GeneSymbol)
exprSet = as.data.frame(exprSet)
write.csv(exprSet,file = "data.csv")

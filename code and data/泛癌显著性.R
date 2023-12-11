library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)
library(rstatix)
data <- read.csv("C:/Users/ZHH/Desktop/data.csv")
library(reshape2)     #  首先加载一下reshape2包
aql <- melt(data)
write.csv(data,file = 'data.csv')
data <- data %>% drop_na()
d2 <- data[-which(expr>0.3),]
library(ggpubr)

# 绘制小提琴图
p <- ggviolin(data, x = "project", y = "expr",
              color = "group", palette = "jco") +
  rotate_x_text(angle = 90)

# 添加显著性标记
p1 <- p + stat_compare_means(aes(group = group), label = "p.signif", label.y = 0.15)
p <- ggboxplot(data, x = "project", y = "expr",
               color = "group", palette = "jco")+
  rotate_x_text(angle = 90) #将x轴肿瘤名称旋转90°展示

p1 <- p + stat_compare_means(aes(group = group), label = "p.signif", label.y =0.15) 



ggsave("vin_plot.pdf",p1, height = 6, width = 10)

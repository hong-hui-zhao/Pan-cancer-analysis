### 学分绩计算 ----

### data input ----
setwd('C:/Users/ZHH/Desktop/TFAP4')
library(reshape2)     #  首先加载一下reshape2包
library(tidyverse)
b <- readxl::read_xlsx("C:/Users/ZHH/Desktop/data3.xlsx")
C <- melt(b,by = "Class")
write.csv(C,file = "data1.csv")
data <- read.csv("C:/Users/ZHH/Desktop/data1.csv",header = T)
credit <- readxl::read_xlsx("C:/Users/ZHH/Desktop/Credit.xlsx")
a <- left_join(data,credit)
a[a==0]<-NA
t <- na.omit(a)

# 然后，我们可以计算每个学生的绩点
t$GPA <- with(t, Grade * Credit / sum(Credit))

d <- t %>%
  group_by(ID) %>%
  summarise(GPA = sum(Grade * Credit) / sum(Credit))

r <- readxl::read_xlsx("C:/Users/ZHH/Desktop/data.xlsx")
dd <- left_join(r,d ,by = "ID")
write.csv(dd,file = "GPA.csv")

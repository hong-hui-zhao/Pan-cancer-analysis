###

### author honghuizhao

### We refer to the code of a public blog post (https://mp.weixin.qq.com/s/HPisZuf0efkoMHwyj_rfIA)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(treemap)
setwd("~")

### load the data for analysis
df <- read.csv("network.csv", header = T, row.names = 1)

### Create nodes and visualize
index <- c("pathway","gene")
nodes<- gather_graph_node(df,index = index, value = "Correlation",root="GO/npathway")
edges <- gather_graph_edge(df,index = index,root="GO/npathway")
graph_df <- tbl_graph(nodes,edges)

### Data visualization
ggraph(graph_df,layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(aes(color=node1.node.branch),alpha=1/3) + 
  geom_node_point(aes(size=node.size,color=node.branch),alpha=1/3) + 
  scale_size(range = c(3,20))+
  scale_size(range = c(0.5,50)) +
  geom_node_text(
    aes(
      x = 1.0175 * x,
      y = 1.0175 * y,
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf,
      color = node.branch
    ),
    size = 4, hjust = 'outward')+
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf,
        color = node.branch),
    fontface="bold",
    size=5,
    family="sans")+
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3))+
  theme_void()+
  theme(legend.position = "none")+
  scale_colour_manual(values =c('#efb306',
                                '#e8351e',
                                '#cd023d',
                                '#4e54ac',
                                '#0f8096'))#修改颜色
  


## function -------------
gather_graph_edge <- function(df,index=NULL,root=NULL){
  require(dplyr)
  if (length(index) < 2){
    stop("please specify at least two index column(s)")
  } else if (length(index)==2){
    data <- df %>% mutate(from=.data[[index[[1]]]]) %>%
      tidyr::unite(to,index,sep="/") %>%
      select(from,to) %>%
      mutate_at(c("from","to"),as.character)
  } else {
    list <- lapply(seq(2,length(index)), function(i){
      dots <- index[1:i]
      df %>% tidyr::unite(from,dots[-length(dots)],sep = "/",remove = F)  %>%
        tidyr::unite(to,dots,sep="/") %>%
        select(from,to) %>%
        mutate_at(c("from","to"),as.character)
    })
    data <- do.call("rbind",list)
  }
  data <- as_tibble(data)
  if (is.null(root)){
    return(data)
  } else {
    root_data <- df %>% group_by(.dots=index[[1]]) %>%
      summarise(count=n()) %>%
      mutate(from=root,to=as.character(.data[[index[[1]]]] )) %>%
      select(from,to)
    rbind(root_data,data)
  }
  
}


gather_graph_node <- function(df,index=NULL,value=tail(colnames(df),1),root=NULL){
  require(dplyr)
  if (length(index) < 2){
    stop("please specify at least two index column(s)")
  } else {
    list <- lapply(seq_along(index), function(i){
      dots <- index[1:i]
      df %>%
        group_by(.dots=dots) %>%
        summarise(node.size=sum(.data[[value]]),
                  node.level=index[[i]],
                  node.count=n()) %>%
        mutate(node.short_name=as.character(.data[[ dots[[length(dots)]] ]]),
               node.branch = as.character(.data[[ dots[[1]]]])) %>%
        tidyr::unite(node.name,dots,sep = "/")
    })
    data <- do.call("rbind",list) %>% as_tibble()
    data$node.level <- factor(data$node.level,levels = index)
    
    if (is.null(root)){
      return(data)
    } else {
      root_data <- data.frame(node.name=root,
                              node.size=sum(df[[value]]),
                              node.level=root,
                              node.count=1,
                              node.short_name=root,
                              node.branch=root,
                              stringsAsFactors = F)
      data <- rbind(root_data,data)
      data$node.level <- factor(data$node.level, levels = c(root,index))
      return(data)
    }
  }
}



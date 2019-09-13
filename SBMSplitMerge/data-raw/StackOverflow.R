## code to prepare `StackOverflow` dataset goes here
## kaggle login required
## library(rvest)
## read_html("https://www.kaggle.com/stackoverflow/stack-overflow-tag-network/downloads/stack_network_links.csv/1")

## retrieved this local version on 27/8/2019

library(SBMSplitMerge)
stack <- read.csv("stack_network_links.csv", stringsAsFactor=F)
nodes <- unique(c(unique(stack$source), unique(stack$target)))
N <- length(nodes)
E <- matrix(0,N,N)
for(i in 1:N){
    for(j in 1:N){
        ind <- (stack$source == nodes[i]) & (stack$target == nodes[j])
        if(any(ind)){
            E[i,j] <- stack$value[ind]
        }
    }
}
colnames(E) <- nodes
rownames(E) <- nodes
StackOverflow <- edges(E)

usethis::use_data(StackOverflow)

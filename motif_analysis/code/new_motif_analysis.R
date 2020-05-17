##########################################################
#This R script is to analysis the motif enrichment and cluster
##########################################################

library(pheatmap)
library(ggplot2)
library(plyr)
library(corrplot)
library(ggseqlogo)
library("ape")
library("ggdendro")
library("dendextend")
library("dplyr")
library("tidyselect")

##########################################################
#calculate the enrichment pvalue of top five motif by Fisher's exact test
##########################################################

human_motif<-read.table(file="motif_human_weight_top5.txt", header = F, sep = "\t")

head(human_motif)

i = 1
Pvalue <- c()
Odd_ratio <- c()
FDR <- c()
for(i in 1:dim(human_motif)[1]){
  Sum = human_motif[i,3]/human_motif[i,4]
  compare<-matrix(floor(c(human_motif[i,3],Sum*0.1,Sum*(1-human_motif[i,4]),Sum*0.9)),nr=2,dimnames=
                    list(c("sites","not sites"),c("motif","random")))
  fisher_test <- fisher.test(compare,alternative = "greater")
  Pvalue <- c(Pvalue, fisher_test$p.value)
  Odd_ratio <- c(Odd_ratio, fisher_test$estimate)
}

for(i in seq(1, dim(human_motif)[1], 5)){
  FDR <- c(FDR, p.adjust(Pvalue[i:(i+4)], method = "fdr"))
}

colnames(human_motif) <- c("motif_id", "motif", "number", "percent")
data1 <- cbind(human_motif, Odd_ratio)
data1 <- cbind(data1, Pvalue)
data1 <- cbind(data1, FDR)
write.table(data1, file = "motif_human_summary_pvalue_top5.txt", row.names = T, col.names = T, sep = "\t")

#####################################################
#cluster the enriched motifs of the RBP with function annotation
#####################################################

human_motif_prob<-read.table(file="motif_human_prob_top5.txt", header = F, sep = "\t")
human_motif_select<-read.table(file="human_motif_filter2.txt", header = T, sep = "\t")

head(human_motif_select)
dim(human_motif_prob)
human_motif_prob[1,1]

data2 = human_motif_select[,1:2]
colnames(human_motif_prob)[1] = "motif_id2"

human_motif_prob = human_motif_prob[human_motif_prob[,1]%in%human_motif_select[,2],]

data3 = merge(data2, human_motif_prob, by = "motif_id2")
data4 = data3[,2:62]
head(data4)

data4[,1] = gsub("_", "#", data4[,1])

motif_human_pro = data4

dis_matrix2 <- matrix(0, dim(motif_human_pro)[1], dim(motif_human_pro)[1])

for(i in 1:dim(motif_human_pro)[1]){
  for(j in 1:dim(motif_human_pro)[1]){
    dis_matrix2[i,j] = cal_distance(motif_human_pro[i,2:61], motif_human_pro[j,2:61])
  }
}

row.names(dis_matrix2) = data4[,1]
colnames(dis_matrix2) = data4[,1]
hdist2 = as.dist(dis_matrix2)
hc2 <- hclust(hdist2, method = "complete")

plot(as.phylo(hc2), type = "fan")
#length(hc2$order)
#write.table(data4[hc2$order,1], "cluster_order.txt")

######################################################
#plot the motif logo figure
######################################################

data3 = motif_human_pro[hc2$order,]
data3 = data4[order(data4[,1]),]

p_list <- list()
for(i in 1:161){
  ll <- as.numeric(data3[i,2:61])*1000
  matt <- matrix(ll, nrow = 6, ncol = 10, byrow = F)
  #out1 = paste(data2[i,1], "seq", sep = "_")
  out1 = data3[i,1]
  matt1 <- matt[1:4,1:10]
  row.names(matt1) <- c("A", "C", "G", "U")
  csl1 <- make_col_scheme(chars = c("A", "C", "G", "U"), cols = c("green","blue","orange","red"))
  myplot1 <- ggseqlogo(matt1, method="bits", namespace=c("A", "C", "G", "U"), col_scheme=csl1) + ylim(0,2) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + ylab(out1)
  p_list[[gaps(i)]] = myplot1
  #out1 = paste(data2[i,1], "str", sep = "_")
  matt2 <- matt[5:6,1:10]
  csl <- make_col_scheme(chars = c("P", "U"), cols = c("brown","purple"))
  row.names(matt2) <- c("P", "U")
  matt2 = entopy_matrix(matt2)
  matt2 = -matt2
  myplot2 <- ggseqlogo(matt2, method="custom", namespace=c("P", "U"), col_scheme=csl) + ylim(-2,0) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + ylab("")
  p_list[[gaps(i)+8]] = myplot2
  #ggsave(paste(out1, "seqmotif.pdf", sep = "_"), myplot)
  #dev.off()
}

p_list[[gaps(162)]] = myplot1
p_list[[gaps(163)]] = myplot1
p_list[[gaps(164)]] = myplot1
p_list[[gaps(165)]] = myplot1
p_list[[gaps(166)]] = myplot1
p_list[[gaps(167)]] = myplot1
p_list[[gaps(168)]] = myplot1
p_list[[gaps(162)+8]] = myplot2
p_list[[gaps(163)+8]] = myplot2
p_list[[gaps(164)+8]] = myplot2
p_list[[gaps(165)+8]] = myplot2
p_list[[gaps(166)+8]] = myplot2
p_list[[gaps(167)+8]] = myplot2
p_list[[gaps(168)+8]] = myplot2


do.call(gridExtra::grid.arrange,c(p_list, ncol=8))

#########################################################

cal_distance_seq <- function(l1, l2){
  l1 = as.numeric(l1)
  l2 = as.numeric(l2)
  #mat1 <- l1[13:48]
  mat1 <- l1[c(13:16,19:22,25:28,31:34,37:40,43:46)]
  for(i in 1:5){
    #mat1 = c(mat1, l2[(i*6-5):(i*6+30)])
    mat1 = c(mat1, l2[c((i*6-5):(i*6-2),(i*6+1):(i*6+4),(i*6+7):(i*6+10),
                        (i*6+13):(i*6+16),(i*6+19):(i*6+22),(i*6+25):(i*6+28))])
  }
  mat1 = matrix(mat1, ncol = 36, byrow = T)
  return(min(as.matrix(dist(as.matrix(mat1)))[2:6,1]))
}


cal_distance <- function(l1, l2){
  l1 = as.numeric(l1)
  l2 = as.numeric(l2)
  mat1 <- l1[13:48]
  for(i in 1:5){
    mat1 = c(mat1, l2[(i*6-5):(i*6+30)])
  }
  mat1 = matrix(mat1, ncol = 36, byrow = T)
  return(min(as.matrix(dist(as.matrix(mat1)))[2:6,1]))
}

entopy_matrix <- function(mat1) {
  posum <- colSums(mat1)
  for (i in 1:dim(mat1)[1]){
    mat1[i,] = mat1[i,]/posum
  }
  mat2 = mat1
  for (i in 1:dim(mat1)[2]){
    entopy_score = cal_entopy(mat1[,i])
    mat2[,i] = entopy_score * mat2[,i]
  }
  return(mat2)
}

cal_entopy <- function(ll) {
  entopy_score = log2(length(ll))
  for(i in 1:length(ll)){
    if(ll[i] != 0){
      entopy_score = entopy_score + ll[i]*log2(ll[i])
    }
  }
  return(entopy_score)
}

gaps<-function(num){
  numa = floor(num/8)
  if(num%%8 == 0){
    return(num + (numa-1)*8)
  }else{
    return(num + numa*8)
  }
}

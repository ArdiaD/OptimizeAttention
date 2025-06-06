#rm(list = ls())
library("ggplot2")
library("reshape2")
library("cowplot")
load("output/EPU/3.rda")
load("output/EPU/3_final_selection.rda")
dim_1 = dim_2 = dim_3 = list()

for(i in 1:length(res)){
  ### Label switching fix due to sorting (for duplicate detection in the algorithm) for graphical representation
  
  if(i <5) {
    dim_1[[i]] = res[[i]]$keywords_results[[2]]
    dim_2[[i]] = res[[i]]$keywords_results[[1]]
    dim_3[[i]] = res[[i]]$keywords_results[[3]]
  } else {
    dim_1[[i]] = res[[i]]$keywords_results[[1]]
    dim_2[[i]] = res[[i]]$keywords_results[[2]]
    dim_3[[i]] = res[[i]]$keywords_results[[3]]
  }
 
}

dim_1[[14]] = rownames(final_res$selection)[final_res$selection[,1] == 1]
dim_2[[14]] = rownames(final_res$selection)[final_res$selection[,2] == 1]
dim_3[[14]] = rownames(final_res$selection)[final_res$selection[,3] == 1]

dim_1_u = sort(unique(unlist(dim_1)))
dim_2_u = sort(unique(unlist(dim_2)))
dim_3_u = sort(unique(unlist(dim_3)))

dim_1_mat = matrix(0,nrow = length(dim_1_u), ncol = length(res)+1)
dim_2_mat = matrix(0,nrow = length(dim_2_u), ncol = length(res)+1)
dim_3_mat = matrix(0,nrow = length(dim_3_u), ncol = length(res)+1)

rownames(dim_1_mat) = dim_1_u
rownames(dim_2_mat) = dim_2_u
rownames(dim_3_mat) = dim_3_u
colnames(dim_1_mat) = colnames(dim_2_mat) = colnames(dim_3_mat) = c(1:14)

for(i in 1:length(res)){
  
  if(i <5) {
    dim_1_mat[match(res[[i]]$keywords_results[[2]],rownames(dim_1_mat)),i] = 1 
    dim_2_mat[match(res[[i]]$keywords_results[[1]],rownames(dim_2_mat)),i] = 1 
    dim_3_mat[match(res[[i]]$keywords_results[[3]],rownames(dim_3_mat)),i] = 1
  } else {
    dim_1_mat[match(res[[i]]$keywords_results[[1]],rownames(dim_1_mat)),i] = 1 
    dim_2_mat[match(res[[i]]$keywords_results[[2]],rownames(dim_2_mat)),i] = 1 
    dim_3_mat[match(res[[i]]$keywords_results[[3]],rownames(dim_3_mat)),i] = 1 
  }
}

dim_1_mat[match(rownames(final_res$selection)[final_res$selection[,1] == 1],rownames(dim_1_mat)),14] = 1
dim_2_mat[match(rownames(final_res$selection)[final_res$selection[,2] == 1],rownames(dim_2_mat)),14] = 1
dim_3_mat[match(rownames(final_res$selection)[final_res$selection[,3] == 1],rownames(dim_3_mat)),14] = 1

########################################

dim_1_mat_melt = melt(dim_1_mat)
dim_2_mat_melt = melt(dim_2_mat)
dim_3_mat_melt = melt(dim_3_mat)

colnames(dim_1_mat_melt) = colnames(dim_2_mat_melt) =colnames(dim_3_mat_melt) =c("token","Epoch","value")
dim_1_mat_melt$token = stringr::str_replace_all(dim_1_mat_melt$token,pattern =" ", "_")
dim_2_mat_melt$token = stringr::str_replace_all(dim_2_mat_melt$token,pattern =" ", "_")
dim_3_mat_melt$token = stringr::str_replace_all(dim_3_mat_melt$token,pattern =" ", "_")

textcol <- "grey40"

p1 <- ggplot(dim_1_mat_melt, aes(x=Epoch, y=token, fill=factor(value)))+ 
  geom_tile(colour="white", size=0.2)+
  scale_y_discrete(expand=c(0, 0))+ ylab(c("")) + xlab(c("")) +
  scale_x_continuous(expand=c(0, 0), breaks=c(1:14),labels = c(1:13,"P"))+
  scale_fill_manual(values=c("white","grey40")) +
  theme_bw(base_size=14)+ guides(fill="none")

p2 <- ggplot(dim_2_mat_melt, aes(x=Epoch, y=token, fill=factor(value)))+ 
  geom_tile(colour="white", size=0.2)+
  scale_y_discrete(expand=c(0, 0))+ ylab(c("")) + xlab(c("")) +
  scale_x_continuous(expand=c(0, 0), breaks=c(1:14),labels = c(1:13,"P"))+
  scale_fill_manual(values=c("white","grey40")) +
  theme_bw(base_size=14)+ guides(fill="none")

p3 <- ggplot(dim_3_mat_melt, aes(x=Epoch, y=token, fill=factor(value)))+ 
  geom_tile(colour="white", size=0.2)+
  scale_y_discrete(expand=c(0, 0))+ ylab(c("")) + xlab(c("Epoch")) +
  scale_x_continuous(expand=c(0, 0), breaks=c(1:14),labels = c(1:13,"P"))+
  scale_fill_manual(values=c("white","grey40")) +
  theme_bw(base_size=14)+ guides(fill="none")

min_converg_iteration = NA
for(i in 1:length(res)){
  
  min_converg_iteration[i] = min(which(res[[i]]$train_loss == min(res[[i]]$train_loss)))
}

min_converg_iteration[14] = min_converg_iteration[13]

min_converg = data.frame(steps = min_converg_iteration, K = c(1:14))

p4 <- ggplot(min_converg, aes(x=K, y=min_converg_iteration))+
  geom_line() +
  theme_bw(base_size=12)+ guides(fill="none") + xlab("Epoch") + ylab("Nb Iterations")+
  scale_x_continuous(expand=c(0.041, 0), breaks=c(1:14),labels = c(1:13,"P"))
p4             

min_loss_iteration = NA
for(i in 1:length(res)){
  
  min_loss_iteration[i] = min(res[[i]]$validation_loss)
}

min_loss_iteration[14] = 0

min_loss = data.frame(steps = min_loss_iteration, K = 1:14)

p5 <- ggplot(min_loss, aes(x=K, y=min_loss_iteration))+
  geom_line() +
  theme_bw(base_size=12)+ guides(fill="none") + xlab("Epoch") + ylab("Val loss")+
  scale_x_continuous(expand=c(0.041, 0), breaks=c(1:14),labels = c(1:13,"P"))
p5  

#### SAVE #####
pdf("figures/EPU_selection.pdf",
    height = 10,
    width = 8,
    paper = "special")

plot_grid(p5,p4,p2,p1,p3, rel_heights = c(1.5,1.5,1.5,4,1), ncol = 1, align="hv")
dev.off()



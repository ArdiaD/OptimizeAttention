rm(list = ls())

load("output/0/CPI_nowcast_PositiveLink3_lambda=0_lambda1=1_lambda3=0_h=1_final_selection_lag_FALSE.rda")
index_1 = final_res$selection

load("output/1/CPI_forcast_PositiveLink3_lambda=0_lambda1=1_lambda3=0_h=1_final_selection_lag_FALSE.rda")
index_3 = final_res$selection

load("output/3/CPI_forcast_PositiveLink3_lambda=0_lambda1=1_lambda3=0_h=3_final_selection_lag_FALSE.rda")
index_5 = final_res$selection

load("output/6/CPI_forcast_PositiveLink3_lambda=0_lambda1=1_lambda3=0_h=6_final_selection_lag_FALSE.rda")
index_7 = final_res$selection

index = list(index_1,index_3,index_5,index_7)

dim_1 = dim_2 = dim_3 = list()

for(i in 1:length(index)){
  ### Label switching fix due to sorting (for duplicate detection in the algorithm) for graphical representation
  
  
  dim_1[[i]] = rownames(index[[i]][index[[i]][,1] == 1,])
  dim_2[[i]] = rownames(index[[i]][index[[i]][,2] == 1,])
  dim_3[[i]] = rownames(index[[i]][index[[i]][,3] == 1,])
  
}

dim_1_u = sort(unique(unlist(dim_1)))
dim_2_u = sort(unique(unlist(dim_2)))
dim_3_u = sort(unique(unlist(dim_3)))

dim_1_mat = matrix(0,nrow = length(dim_1_u), ncol = length(index))
dim_2_mat = matrix(0,nrow = length(dim_2_u), ncol = length(index))
dim_3_mat = matrix(0,nrow = length(dim_3_u), ncol = length(index))

rownames(dim_1_mat) = dim_1_u
rownames(dim_2_mat) = dim_2_u
rownames(dim_3_mat) = dim_3_u
colnames(dim_1_mat) = colnames(dim_2_mat) = colnames(dim_3_mat) = c(0,1,2,3)

for(i in 1:length(index)){
  
  
  dim_1_mat[match(rownames(index[[i]][index[[i]][,1] == 1,]),rownames(dim_1_mat)),i] = 1 
  dim_2_mat[match(rownames(index[[i]][index[[i]][,2] == 1,]),rownames(dim_2_mat)),i] = 1 
  dim_3_mat[match(rownames(index[[i]][index[[i]][,3] == 1,]),rownames(dim_3_mat)),i] = 1
  
}

library("ggplot2")

###################

dim_1_mat_melt = reshape2::melt(dim_1_mat)
dim_2_mat_melt = reshape2::melt(dim_2_mat)
dim_3_mat_melt = reshape2::melt(dim_3_mat)

colnames(dim_1_mat_melt) = colnames(dim_2_mat_melt) =colnames(dim_3_mat_melt) =c("token","index","value")
textcol <- "grey40"

dim_1_mat_melt$token = stringr::str_replace_all(dim_1_mat_melt$token,pattern =" ", "_")
dim_2_mat_melt$token = stringr::str_replace_all(dim_2_mat_melt$token,pattern =" ", "_")
dim_3_mat_melt$token = stringr::str_replace_all(dim_3_mat_melt$token,pattern =" ", "_")

p1 <- ggplot(dim_1_mat_melt, aes(x=index, y=token, fill=factor(value)))+ 
  geom_tile(colour="white", size=0.2)+
  scale_y_discrete(expand=c(0, 0))+ ylab(c("")) + xlab(c("")) +  
  scale_x_continuous(expand=c(0, 0), breaks=c(0,1,2,3),labels = c(0,1,3,6))+
  scale_fill_manual(values=c("white","grey40")) +
  theme_bw(base_size=14)+ guides(fill="none")

p2 <- ggplot(dim_2_mat_melt, aes(x=index, y=token, fill=factor(value)))+ 
  geom_tile(colour="white", size=0.2)+
  scale_y_discrete(expand=c(0, 0))+ ylab(c("")) + xlab(c("")) +  
  scale_x_continuous(expand=c(0, 0), breaks=c(0,1,2,3),labels = c(0,1,3,6))+
  scale_fill_manual(values=c("white","grey40")) +
  theme_bw(base_size=14)+ guides(fill="none")

p3 <- ggplot(dim_3_mat_melt, aes(x=index, y=token, fill=factor(value)))+ 
  geom_tile(colour="white", size=0.2)+
  scale_y_discrete(expand=c(0, 0))+ ylab(c("")) + xlab(c("Horizon")) + 
  scale_x_continuous(expand=c(0, 0), breaks=c(0,1,2,3),labels = c(0,1,3,6))+
  scale_fill_manual(values=c("white","grey40")) +
  theme_bw(base_size=14)+ guides(fill="none")

library("cowplot")

#### SAVE #####
pdf("figures/forecast_selection.pdf",
    height = 12,
    width = 8,
    paper = "special")
plot_grid(p1, p2, p3, rel_heights = c(2,2,2), ncol = 1,align="v")
dev.off()

source(file = "functions/compute_attention_measure.R")

load(file = "data/dfm_filtered_resolved_unigram_ManualvocFilt_sentiment_accronym.rda") # for IJF Editor only
#load(file = "data/dfm_filtered_resolved_unigram_ManualvocFilt_sentiment_accronym_pseudo.rda") # for other users

#remove duplicate form old and new data
dat = dat[-which(unlist(lapply(stringr::str_split(dat@docvars$docname_,"_"), function(x) x[2])) == "WSJ" & dat@docvars$date >= as.Date("2021-08-01")),]
dat = dat[,order(colnames(dat))]

dat = dat[dat@docvars$sentiment <= quantile(dat@docvars$sentiment,0.3333),]
dat@docvars$sentiment = (-dat@docvars$sentiment)

freq_y = "month"
dfm_sub = dat
dfm_sub$frequency_date = floor_date(x = dfm_sub@docvars$date, unit = freq_y)

mat_doc_freq = Matrix.utils::dMcast(dfm_sub@docvars,
                                    formula = frequency_date  ~ docid_)

files = c("output/0/CPI_nowcast_PositiveLink3_lambda=0_lambda1=1_lambda3=0_h=1_final_selection_lag_FALSE.rda",
          "output/1/CPI_forcast_PositiveLink3_lambda=0_lambda1=1_lambda3=0_h=1_final_selection_lag_FALSE.rda",
          "output/3/CPI_forcast_PositiveLink3_lambda=0_lambda1=1_lambda3=0_h=3_final_selection_lag_FALSE.rda",
          "output/6/CPI_forcast_PositiveLink3_lambda=0_lambda1=1_lambda3=0_h=6_final_selection_lag_FALSE.rda")

get_topic = function(file, dat,mat_doc_freq, findk = TRUE, k){
  
  load(file)
  dfm_sub = dat[, rownames(final_res$selection)]
  freq_y = "month"
  dfm_sub$frequency_date = floor_date(x = dfm_sub@docvars$date, unit = freq_y)
  y = compute_attention_measure(dfm = dfm_sub,
                                omega = final_res$selection,
                                gamma = t(matrix(rep(1,nrow(final_res$selection)))),
                                mat_doc_freq = mat_doc_freq,freq = freq_y,
                                sentiment = FALSE
  )
  
  
  sub_dat = dat[which(y$sel),]
  
  sub_dat_stm = stm::asSTMCorpus(sub_dat)
  
  set.seed(1234)
  if(findk){
    topics = stm::searchK(sub_dat_stm$documents,vocab = sub_dat_stm$vocab,K = k)
  } else {
    topics = stm::stm(sub_dat_stm$documents,vocab = sub_dat_stm$vocab,K = k)
  }
  return(list(topics = topics, sel = y$sel,attention = y$unormalized_attention))
}

top_1_k = get_topic(files[1], dat = dat, mat_doc_freq = mat_doc_freq, findk = TRUE, k = 4:20)
top_2_k = get_topic(files[2], dat = dat, mat_doc_freq = mat_doc_freq, findk = TRUE, k = 4:20)
top_3_k = get_topic(files[3], dat = dat, mat_doc_freq = mat_doc_freq, findk = TRUE, k = 4:20)
top_4_k = get_topic(files[4], dat = dat, mat_doc_freq = mat_doc_freq, findk = TRUE, k = 4:20)

top_dat = data.frame(c(unlist(top_1_k$topics$results$exclus),unlist(top_2_k$topics$results$exclus),unlist(top_3_k$topics$results$exclus),unlist(top_4_k$topics$results$exclus)),
           c(unlist(top_1_k$topics$results$semcoh),unlist(top_2_k$topics$results$semcoh),unlist(top_3_k$topics$results$semcoh),unlist(top_4_k$topics$results$semcoh)),rep(4:20,4),sort(rep(c(0,1,3,6),17)))

colnames(top_dat) = c("exclus","semcoh","topics","horizon")
top_dat$horizon = as.factor(top_dat$horizon)

ggplot2::ggplot(top_dat,aes(x = exclus,y = semcoh,group = horizon,col = horizon,label = topics))+
  geom_line() +
  geom_label()

top_1 = get_topic(files[1],dat = dat, mat_doc_freq = mat_doc_freq, findk = FALSE,k = 8)
top_2 = get_topic(files[2],dat = dat, mat_doc_freq = mat_doc_freq, findk = FALSE,k = 10)
top_3 = get_topic(files[3],dat = dat, mat_doc_freq = mat_doc_freq, findk = FALSE,k = 11)
top_4 = get_topic(files[4],dat = dat, mat_doc_freq = mat_doc_freq, findk = FALSE,k = 8)

##################################################################################
topics = list(top_h_0 = top_1, top_h_1 = top_2, top_h_3 = top_3,top_h_6 = top_4)
save(topics,file = "output/topics.rda")

rio::export(rbind(cbind("h=0",round(colMeans(topics$top_h_0$topics$theta)*100,2),apply(stm::sageLabels(topics$top_h_0$topics,n = 20)$marginal$frex,MARGIN = 1,FUN = function(x) paste0(x,collapse = ", "))),
      cbind("h=1",round(colMeans(topics$top_h_1$topics$theta)*100,2),apply(stm::sageLabels(topics$top_h_1$topics,n = 20)$marginal$frex,MARGIN = 1,FUN = function(x) paste0(x,collapse = ", "))),
      cbind("h=3",round(colMeans(topics$top_h_3$topics$theta)*100,2),apply(stm::sageLabels(topics$top_h_3$topics,n = 20)$marginal$frex,MARGIN = 1,FUN = function(x) paste0(x,collapse = ", "))),
      cbind("h=6",round(colMeans(topics$top_h_6$topics$theta)*100,2),apply(stm::sageLabels(topics$top_h_6$topics,n = 20)$marginal$frex,MARGIN = 1,FUN = function(x) paste0(x,collapse = ", ")))),file = "tables/topics.csv")




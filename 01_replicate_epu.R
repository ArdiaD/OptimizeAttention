#rm(list = ls())
library("quanteda")
library("doRNG")
library("doParallel")
library("caret")
library("GA")
library("crayon")
library("doFuture")
library("future")
library("Matrix")
library("data.table")
library("lubridate")
library("grr")
library("Matrix.utils")
load("data/wv_keywords_fintext.rda")
source(file = "functions/custom_CX_Mut.R")
source(file = "functions/custom_ga.R")
source(file = "functions/lexicon_expansion.R")
source(file = "functions/compute_attention_measure.R")
source(file = "functions/ga_attention_all.R")
source(file = "functions/loss_functions.R")

#keep those parameters for replication

monitor = TRUE
seed = 1234
lambda_0 = 0
lambda_1 = 0
lambda_2_base = 0
lag = NULL
max_word = 15
maxiter = 200
ncore = 6
batch_size = 100
popSize = ncore*batch_size
bootstrap = FALSE
CV = TRUE
start_date = "2000-01-01"
load(file = "data/dfm_filtered_resolved_unigram_ManualvocFilt_sentiment_accronym.rda")
dat = dat[1:839484,]
sign_coef = 1
conditioning_variable = NULL
f_transformation = function(x) {x}

#remove duplicate form old and new data
dat = dat[-which(unlist(lapply(stringr::str_split(dat@docvars$docname_,"_"), function(x) x[2])) == "WSJ" & dat@docvars$date >= as.Date("2021-08-01")),]
dat = dat[,order(colnames(dat))]

freq_y = "month"
dfm_sub = dat
dfm_sub$frequency_date = floor_date(x = dfm_sub@docvars$date, unit = freq_y)

mat_doc_freq = Matrix.utils::dMcast(dfm_sub@docvars,
                                    formula = frequency_date  ~ docid_)

keywords <- colnames(dat)

# creating EPU index from data

keywords_EPU_baker <- list(
  E = c("economy", "economic"),
  P = c("congress", "legislation","white house","regulation","federal reserve","deficit"),
  U = c("uncertain","uncertainty")
)

dfm_sub = dat[, unlist(keywords_EPU_baker)]
freq_y = "month"
dfm_sub$frequency_date = floor_date(x = dfm_sub@docvars$date, unit = freq_y)

mat_doc_freq = Matrix.utils::dMcast(dfm_sub@docvars,
                                    formula = frequency_date  ~ docid_)

omega_EPU = matrix(0,nrow = length(unlist(keywords_EPU_baker)), ncol = length(keywords_EPU_baker))
v_loc <- rep(NA, length(unlist(keywords_EPU_baker)))
names(v_loc) <- unlist(keywords_EPU_baker)
v_loc_list = list()

for (i in 1:length(keywords_EPU_baker)) {
  v_loc[unlist(keywords_EPU_baker) %in% keywords_EPU_baker[[i]]] <- i
  v_loc_list[[i]] = v_loc
}

for (j in 1:length(v_loc_list)) {
  for (i in 1:length(v_loc_list[[j]])) {
    omega_EPU[i,v_loc_list[[j]][i]] = 1
  }
}

rownames(omega_EPU) = unlist(keywords_EPU_baker)

y = compute_attention_measure(dfm = dfm_sub,
                              omega = omega_EPU,
                              gamma = t(matrix(rep(1,nrow(omega_EPU)))),
                              mat_doc_freq = mat_doc_freq,freq = freq_y,
                              sentiment = FALSE
)$unormalized_attention


y = data.frame(date = rownames(y), value = y)

#########################################################

y_fit = y[y$date >= start_date,]
set.seed(1234)

y_valid = y_fit
y_train = y_fit

y_test = y[y$date >= start_date,]

for (ndim in 3:3) {
  set.seed(seed)
  keywords = colnames(dat)
  res = list()
  l = 1
  suggestion_gamma = NULL
  for(i in ndim:max_word){
    
    if(i >= (ndim*2)){
      lambda_2 = lambda_2_base
    } else {
      lambda_2 = 0
    }
    
    lambda_ga = ga_attention_all(dfm = dat,
                                 y_train = y_train,
                                 y_valid = y_valid,
                                 y_test = y_test,
                                 x = NULL,
                                 freq_y = freq_y,
                                 keywords = keywords,
                                 bootstrap = bootstrap,
                                 CV = CV,
                                 maxiter = maxiter,
                                 popSize = popSize,
                                 monitor = monitor,
                                 crossover = custom_crossover,
                                 mutation = custom_mutation,
                                 ndim = ndim,
                                 ncore = ncore,
                                 k = i,
                                 batch_size = batch_size,
                                 lambda_0 = lambda_0,
                                 lambda_1 = lambda_1,
                                 lambda_2 = lambda_2,
                                 suggestion_gamma = suggestion_gamma,
                                 sentiment = FALSE,
                                 sign = sign_coef,
                                 f_transformation = f_transformation)
   
    #setting up suggestion for old results
    
    
    omega = matrix(0,nrow = length(keywords)*ndim, ncol = ndim)
    start = 1
    for(f in 1:ncol(omega)){
      omega[start:(start+length(keywords)-1),f] = 1
      start = start + length(keywords)
    }
    
    rownames(omega) = rep(keywords, ndim)
    
    #fix cases where two solution give the same solution
    lambda_ga$sol = lapply(lambda_ga$sol,FUN = function(x) x[1,,drop = FALSE])
    keywords_results_omega = omega[which(lambda_ga$sol[[which.min(lambda_ga$validation_loss)]] == 1),, drop =FALSE]
    key = list()
    for(k in 1:ndim){
      key[[k]] = rownames(keywords_results_omega)[keywords_results_omega[,k] ==1]
    }
    keywords_results = key
    suggestion_gamma = lambda_ga$sol[[which.min(lambda_ga$validation_loss)]]
    
    sol = do.call(rbind, lambda_ga$sol)
    sol = Matrix::Matrix(sol,sparse = TRUE)
    lambda_ga$sol = sol
    
    res[[l]] = list(train_loss = lambda_ga$train_loss,
                    validation_loss = lambda_ga$validation_loss,
                    sol = sol,
                    keywords = keywords,
                    keywords_results = keywords_results, 
                    keywords_results_omega = keywords_results_omega, 
                    ndim = ndim)
    l = l+1
    
    mk = match(colnames(dat),colnames(wv_keywords))
    mk = mk[!is.na(mk)]
    keywords_expanded = near_terms(keywords = keywords_results,
                                   keep = 0.2,
                                   wordvectors = wv_keywords[,mk],
                                   keep.score = FALSE,
                                   keep.original = TRUE)
    keywords_expanded = unique(unlist(keywords_expanded))
    keywords = keywords_expanded
    keywords = sort(keywords)
    suggestion_gamma = matrix(0,nrow = 1, ncol = length(keywords)*ndim)
    
    for(k in 1:ndim){
      for(w in 1:length(key[[k]])){
        suggestion_gamma[which(!is.na(match(rep(keywords,ndim),key[[k]][w])))[k]] = 1
      }
    }
  }

  
  save(res,file = paste0("output/EPU/",ndim,".rda"))
  ##################################################################
  
  ########## Pruning ############
  
  intermediary_omega = res[[length(res)]]$keywords_results_omega
  n <- nrow(intermediary_omega)
  l <- rep(list(0:1), n)
  l = expand.grid(l)
  l = l[-1,]
  dfm_sub = dat[,rownames(intermediary_omega)]
  xx = NA
  batches = split(1:nrow(l), ceiling(seq_along(1:nrow(l))/300))
  cl <- parallel::makeCluster(ncore)
  registerDoParallel(cl)
  for(i in 1:length(batches)){
    cat(i)
    news_index = compute_attention_measure(dfm = dfm_sub,
                                           omega = res[[length(res)]]$keywords_results_omega,
                                           gamma = l[batches[[i]],, drop = FALSE],
                                           mat_doc_freq = mat_doc_freq,
                                           freq = freq_y,
                                           sentiment = FALSE
    )$unormalized_attention
    zeros = colSums(news_index<= 0.001)/nrow(news_index)
    for(j in 1:ncol(news_index)){
      xx[batches[[i]][j]] = penalized_MSEK_card(y = y_fit,
                                                news_variable = f_transformation(news_index[,j]),
                                                x = lag,
                                                zeros = zeros[j],
                                                bootstrap = FALSE,
                                                CV = TRUE,
                                                lambda_0 = lambda_0,
                                                lambda_1 = lambda_1,
                                                lambda_2 = lambda_2,
                                                sign = sign_coef,
                                                omega_select = res[[length(res)]]$keywords_results_omega[unlist(l[batches[[i]],, drop = FALSE][j,])==1,, drop = FALSE])
    }
  }
  stopCluster(cl)
  final_omega = res[[length(res)]]$keywords_results_omega[l[which.min(xx),]==1,]
  
  if(is.vector(final_omega)){
    
    final_omega = as.matrix(final_omega)
  }
  dfm_sub = dat[,rownames(final_omega)]
  final_res = list(selection = final_omega,
                   index =   compute_attention_measure(dfm = dfm_sub,
                                                       omega = final_omega,
                                                       gamma = t(matrix(rep(1,nrow(final_omega)))),
                                                       mat_doc_freq = mat_doc_freq,
                                                       freq = freq_y,
                                                       sentiment = FALSE),
                   training_data = y_fit,
                   test_data = y_fit,
                   error_select = min(xx),
                   original_data = y_fit)
  
  final_res$error_training = MSE(y = final_res$training_data, news_variable = f_transformation(final_res$index$unormalized_attention[,1]), x = NULL)$error
  final_res$error_test = MSE(y = final_res$test_data, news_variable = f_transformation(final_res$index$unormalized_attention[,1]), x = NULL)$error
  
  save(final_res,file = paste0("output/EPU/",ndim,"_final_selection.rda"))
}


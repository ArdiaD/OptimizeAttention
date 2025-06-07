#rm(list = ls())
library("quanteda")
library("doRNG")
library("doParallel")
library("caret")
library("GA")
library("crayon")
library("doFuture")
library("future")
load("data/wv_keywords_fintext.rda")
source(file = "functions/custom_CX_Mut.R")
source(file = "functions/custom_ga.R")
source(file = "functions/lexicon_expansion.R")
source(file = "functions/compute_attention_measure.R")
source(file = "functions/ga_attention_all.R")
source(file = "functions/loss_functions.R")
library("Matrix")
library("data.table")
library("lubridate")

# Set parameters and run full script for replication

if (!exists("h")) {
	h = 3 # h = 1, h = 3, h = 6 for replication	
}
if (!exists("lambda_1")) {
	lambda_1 = 1 # 0 or 1 depending on the setup for replication
}
if (!exists("use_lag")) {
	use_lag = TRUE # TRUE or FALSE depending on replication setup
}
monitor = TRUE # keep that to TRUE for replication
seed = 1234 # keep that to 1234 for replication
lambda_0 = 0 # keep that to 0 for replication
lambda_2_base = 0 # keep that to 0 for replication
max_word = 15 # keep that to 15 for replication
maxiter = 200 # keep that to 200 for replication
ncore = 6 # keep that to 6 for replication
batch_size = 100 # keep that to 100 for replication
popSize = ncore*batch_size 
bootstrap = FALSE
CV = TRUE
start_date = "2000-01-01" # keep that to 2000-01-01 for replication
load(file = "data/dfm_filtered_resolved_unigram_ManualvocFilt_sentiment_accronym.rda") # for IJF Editor only
#load(file = "data/dfm_filtered_resolved_unigram_ManualvocFilt_sentiment_accronym_pseudo.rda") # for other users
dat = dat[1:839484,]
sign_coef = 1 # keep that to 1 for replication
test_thresh_oos = "2019-12-31" # keep that to 2019-12-31 for replication

#################################################################################

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

keywords <- colnames(dat)

# y series
fred = rio::import("data/CPILFESL.csv")
CPI = fred
CPI$DATE = as.Date(parsedate::parse_date(CPI$DATE))
CPI = CPI[,c("DATE","CPI")]
colnames(CPI) = c("date","value")
CPI$value[(h+1):nrow(CPI)] = (1200/h)*log(CPI$value[(1+h):nrow(CPI)]/CPI$value[1:(nrow(CPI)-h)])
CPI = CPI[-1:-h,]
freq_y = "month"

#lag data
lag = data.frame(date = CPI$date[1:(nrow(CPI)-h)] , CPI$value[1:(nrow(CPI)-h)])
lag = lag[lag$date >= "2000-01-01" & lag$date <= "2024-08-01",]
CPI = CPI[CPI$date >= "2000-01-01" & CPI$date <= "2024-08-01",]
CPI$value[1:(nrow(CPI)-h)] = CPI$value[(h+1):nrow(CPI)]
CPI = CPI[1:nrow(lag),]

y_fit = CPI[CPI$date <= test_thresh_oos & CPI$date >= "2000-01-01",]
y_train = y_fit

y_test = CPI[CPI$date > test_thresh_oos,]
BEIR = read.csv("data/T5YIEM.csv")
colnames(lag) = c("date","value")
colnames(BEIR) = c("date","value")

conditioning_variable = rbind(data.frame(date = lag$date[1:36],value = lowess(lag[1:36,"value"])$y),BEIR)
f_transformation = function(x){ 
  sqrt(x)*sign(conditioning_variable[match(as.Date(names(x)),conditioning_variable$date),"value"] - 2)
}
################### ITERRATIVE ########################


##################################################################

if(isTRUE(use_lag)){
  lag = lag
} else {
  lag = NULL
} 

for(ndim in 3:3){
  set.seed(1234)
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
    
    set.seed(seed)
    lambda_ga = ga_attention_all(dfm = dat,
                                 y_train = y_train,
                                 y_valid = y_valid,
                                 y_test = y_test,
                                 x = lag,
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
  save(res,file = paste0("output/",h,"/CPI_forcast_PositiveLink",ndim,"_lambda=",lambda_0,"_lambda1=",lambda_1,"_lambda3=",lambda_2_base,"_h=",h,"_lag_",use_lag,".rda"))
  
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
    zeros = colSums(news_index <= 0.001)/nrow(news_index)
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
                   test_data = y_test,
                   error_select = min(xx),
                   original_data = CPI)
  
  final_res$error_training = MSE(y = final_res$training_data, news_variable = f_transformation(final_res$index$unormalized_attention[,1]), x = lag)$error
  final_res$error_test = MSE(y = final_res$test_data, news_variable = f_transformation(final_res$index$unormalized_attention[,1]), x = lag)$error
  
  save(final_res,file = paste0("output/",h,"/CPI_forcast_PositiveLink",ndim,"_lambda=",lambda_0,"_lambda1=",lambda_1,"_lambda3=",lambda_2_base,"_h=",h,"_final_selection_lag_",use_lag,".rda"))
}


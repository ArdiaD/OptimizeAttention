rm(list = ls())

# Bootstrap procedure for significance test
DM_boot = function(loss,B){
  d_loss = loss[,2] - loss[,1]
  d_loss_demeaned = d_loss - mean(d_loss)
  
  
  
  n = length(d_loss)
  k = ar(d_loss)$order
  k = max(k,2)
  v = ceiling(length(d_loss)/k)
  
  indexes_b = lapply(1:B, MCS:::boot.block, v = v, n = n, 
                     k = k)
  
  d_mean_resampled = lapply(indexes_b, 
                            
                            function(x, d) {
                              mean(d[x])
                            },
                            d = d_loss_demeaned)
  
  d_HAC_resampled = lapply(indexes_b, 
                           function(x, d) {
                             nse::nse.andrews(x = d[x], lag.prewhite = NULL)
                           },
                           d = d_loss_demeaned)
  
  t_b_DM = unlist(d_mean_resampled)/unlist(d_HAC_resampled)
  d_mean = mean(d_loss)
  d_sd = nse::nse.andrews(x = d_loss, lag.prewhite = NULL)
  t_DM = d_mean / d_sd
  
  return(sum(t_b_DM <= t_DM)/B)
}

############ For h = 1, 3, 6 ##############
library("quanteda")
library("GA")
library("crayon")
library("lubridate")
library("nse")
source(file = "functions/compute_attention_measure.R")
h_list = c(1,3,6)
res = matrix(NA, nrow = 10,ncol = 4)
confidence =  matrix(0, nrow = 10,ncol = 4)
DM_p_list = vector(mode = "list",length = 4)

for(k in 1:length(h_list)){
  h = h_list[k]
  # y series
  fred = rio::import("data/CPILFESL.csv")
  CPI = fred
  CPI$DATE = as.Date(parsedate::parse_date(CPI$DATE))
  CPI = CPI[,c("DATE","CPI")]
  colnames(CPI) = c("date","value")
  CPI$value[(h+1):nrow(CPI)] = (1200/h)*log(CPI$value[(1+h):nrow(CPI)]/CPI$value[1:(nrow(CPI)-h)])
  CPI = CPI[-1:-h,]
  freq_y = "month"
  
  # threshold for oos analysis
  test_thresh_oos = "2019-12-01"
  
  lag = data.frame(date = CPI$date[1:(nrow(CPI)-h)] , CPI$value[1:(nrow(CPI)-h)])
  lag = lag[lag$date >= "2000-01-01" & lag$date <= "2024-08-01",]
  CPI = CPI[CPI$date >= "2000-01-01" & CPI$date <= "2024-08-01",]
  CPI$value[1:(nrow(CPI)-h)] = CPI$value[(h+1):nrow(CPI)]
  CPI = CPI[1:nrow(lag),]
  
  
  y_fit = CPI[CPI$date >= "2000-01-01",]
  
  set.seed(1234)
  
  y_test = CPI[CPI$date > test_thresh_oos,]
  BEIR = read.csv("data/T5YIEM.csv")
  colnames(lag) = c("date","value")
  colnames(BEIR) = c("date","value")
  
  conditioning_variable = rbind(data.frame(date = lag$date[1:36],value = lowess(lag[1:36,"value"])$y),BEIR)
  f_transformation = function(x){ 
  
    sqrt(x)*sign(conditioning_variable[match(as.Date(names(x)),conditioning_variable$date),"value"] - 2)
  }
  
  start_OOS = which(y_fit$date == as.Date(test_thresh_oos)) +1
  err = matrix(NA, nrow = length(start_OOS:nrow(y_fit)),ncol = 10)
  dat = cbind(y = y_fit$value, lag = lag$value)
  dat = as.data.frame(dat)
  
  for(i in start_OOS:nrow(y_fit)){
    
    reg = lm(formula = y ~ 1, data = dat[1:(i-1),])
    err[i-start_OOS+1,1] = (predict(reg,dat[i,]) - dat$y[i])^2
    
  }
  
  
  files = list.files(paste0("output/",h))
  files = files[stringr::str_detect(files,pattern = "final")]
  
  for(j in 1:length(files)){
    load(paste0("output/",h,"/",files[j]))
    dat = cbind(y = y_fit$value, 
                lag = lag$value,
                attention = f_transformation(final_res$index$unormalized_attention[match(as.character(y_fit$date), 
                                                                                         rownames(final_res$index$unormalized_attention)),]))
    dat = as.data.frame(dat)
    for(i in start_OOS:nrow(y_fit)){
      
      reg = lm(formula = y ~ 1 + attention, data = dat[1:(i-1),])
      err[i-start_OOS+1,j+2] = (predict(reg,dat[i,]) - dat$y[i])^2
    }
  }
  
  for(i in start_OOS:nrow(y_fit)){
    reg = lm(formula = y ~ 1 + lag, data = dat[1:(i-1),])
    err[i-start_OOS+1,2] = (predict(reg,dat[i,]) - dat$y[i])^2
  }
  
  for(j in 1:length(files)){
    load(paste0("output/",h,"/",files[j]))
    dat = cbind(y = y_fit$value, 
                lag = lag$value,
                attention = f_transformation(final_res$index$unormalized_attention[match(as.character(y_fit$date),
                                                                                         rownames(final_res$index$unormalized_attention)),]))
    dat = as.data.frame(dat)
    for(i in start_OOS:nrow(y_fit)){
      
      reg = lm(formula = y ~ 1 + lag + attention, data = dat[1:(i-1),])
      err[i-start_OOS+1,j+length(files)+2] = (predict(reg, dat[i,]) - dat$y[i])^2
    }
  }
  
  
   res[,k+1] = sqrt(colMeans(err))

  set.seed(1234)
  DM_p= matrix(NA,10,10)
  for(h in 1:10){
    for(l in 1:10) {
      if( l != h ){
        DM_p[h,l] = DM_boot(err[,c(h,l)],B = 2000)
      }
    }
  }
  colnames(DM_p) = rownames(DM_p) = c("mean","AR(1)","l1=0","l1=0 lag","l1=1","l1=1 lag","AR(1) l1=0"," AR(1) l1=0 lag","AR(1) l1=1","AR(1) l1=1 lag")
  DM_p_list[[k+1]] = DM_p
}

############ For h = 0 ##############

h = 1 # for lag of nowcast

fred = rio::import("data/CPILFESL.csv")
CPI = fred
CPI$DATE = as.Date(parsedate::parse_date(CPI$DATE))
CPI = CPI[,c("DATE","CPI")]
colnames(CPI) = c("date","value")
CPI$value[(h+1):nrow(CPI)] = (1200/h)*log(CPI$value[(1+h):nrow(CPI)]/CPI$value[1:(nrow(CPI)-h)])
CPI = CPI[-1:-h,]
freq_y = "month"

# threshold for oos analysis

lag = CPI
lag$value[2:(nrow(CPI))] = lag$value[1:(nrow(CPI)-1)]
CPI = CPI[1:nrow(lag),]
CPI = CPI[CPI$date >= "2000-01-01" & CPI$date <= "2024-08-01",]
lag = lag[lag$date >= "2000-01-01" & lag$date <= "2024-08-01",]

y_fit = CPI[CPI$date >= "2000-01-01",]

set.seed(1234)

y_test = CPI[CPI$date > test_thresh_oos,]
BEIR = read.csv("data/T5YIEM.csv")
colnames(lag) = c("date","value")
colnames(BEIR) = c("date","value")

conditioning_variable = rbind(data.frame(date = lag$date[1:36],value = lowess(lag[1:36,"value"])$y),BEIR)

f_transformation = function(x){ 
  sqrt(x)*sign(conditioning_variable[match(as.Date(names(x)),conditioning_variable$date),"value"] - 2)
}

start_OOS = which(y_fit$date == as.Date(test_thresh_oos)) +1
err = matrix(NA, nrow = length(start_OOS:nrow(y_fit)),ncol = 10)
dat = cbind(y = y_fit$value, lag = lag$value)
dat = as.data.frame(dat)

for(i in start_OOS:nrow(y_fit)){
  
  reg = lm(formula = y ~ 1 , data = dat[1:(i-1),])
  err[i-start_OOS+1,1] = (predict(reg,dat[i,]) - dat$y[i])^2
  
}

for(i in start_OOS:nrow(y_fit)){
  
  reg = lm(formula = y ~ 1 + lag, data = dat[1:(i-1),])
  err[i-start_OOS+1,2] = (predict(reg,dat[i,]) - dat$y[i])^2
  
}

files = list.files(paste0("output/",0))
files = files[stringr::str_detect(files,pattern = "final")]

for(j in 1:length(files)){
  load(paste0("output/",0,"/",files[j]))
  dat = cbind(y = y_fit$value, 
              lag = lag$value,
              attention = f_transformation(final_res$index$unormalized_attention[match(as.character(y_fit$date), 
                                                                                       rownames(final_res$index$unormalized_attention)),]))
  dat = as.data.frame(dat)
  for(i in start_OOS:nrow(y_fit)){
    
    reg = lm(formula = y ~ 1 + attention, data = dat[1:(i-1),])
    err[i-start_OOS+1,j+2] = (predict(reg,dat[i,]) - dat$y[i])^2
  }
}

for(j in 1:length(files)){
  load(paste0("output/",0,"/",files[j]))
  dat = cbind(y = y_fit$value, 
              lag = lag$value,
              attention = f_transformation(final_res$index$unormalized_attention[match(as.character(y_fit$date),
                                                                                       rownames(final_res$index$unormalized_attention)),]))
  dat = as.data.frame(dat)
  for(i in start_OOS:nrow(y_fit)){
    
    reg = lm(formula = y ~ 1 + lag + attention, data = dat[1:(i-1),])
    err[i-start_OOS+1,j+length(files)+2] = (predict(reg, dat[i,]) - dat$y[i])^2
  }
}

res[,1] = sqrt(colMeans(err))

set.seed(1234)

DM_p= matrix(NA,10,10)
for(h in 1:10){
  for(l in 1:10) {
    if( l != h ){
      DM_p[h,l] = DM_boot(err[,c(h,l)],B = 2000)
    }
  }
}

colnames(DM_p) = rownames(DM_p) = c("mean","AR(1)","l1=0","l1=0 lag","l1=1","l1=1 lag","AR(1) l1=0"," AR(1) l1=0 lag","AR(1) l1=1","AR(1) l1=1 lag")
DM_p_list[[1]] = DM_p

##################

res = round(res,2)


final_res = res
final_res[2:5,] = res[3:6,]
final_res[6,] = res[2,]

rownames(final_res) = rownames(DM_p_list[[1]])[c(1,3,4,5,6,2,7,8,9,10)]

# When reporting, only report the model that are equal within the same base model
# for instance if the best model is under the AR(1) benchmark, only report the equal model with AR(1) benchmark

#AR(1) l1=1 lag
which(DM_p_list[[1]][,which.min(res[,1])] > 0.1)
#AR(1) l1=1 /=/ AR(1) l1=1 lag
which(DM_p_list[[2]][,which.min(res[,2])] > 0.1)
#l1=1 /=/ l1=1 lag 
which(DM_p_list[[3]][,which.min(res[,3])] > 0.1)
# l1=0 lag /=/ l1=1 lag (not including from equivalence from AR(1) type model, still beat baseline mean and AR(1) only model)
which(DM_p_list[[4]][,which.min(res[,4])] > 0.1)

print(xtable::xtable(final_res),file = "tables/results.txt")

MSE = function(news_variable = NULL, y, x = NULL,bootstrap = FALSE, CV = FALSE){

  news_variable = as.matrix(news_variable)
  id = match(rownames(news_variable),as.character(y$date))
  news_variable = news_variable[!is.na(id),]
  news_variable = as.matrix(news_variable)
  y = y[id[!is.na(id)],2]
  y = scale(y)
 
  if(bootstrap){
 
    theta.fit <- function(x,y){lsfit(x,y)}
    theta.predict <- function(fit,x){
      cbind(1,x)%*%fit$coef
    }
    sq.err <- function(y,yhat) { (y-yhat)^2}
    
    if(!is.null(x)){
      id = match(rownames(news_variable),as.character(x$date))
      x = x[id,-1]
      x = cbind(x,news_variable)
      x = as.matrix(x)
  
      if(CV){
        old <- .Random.seed
        set.seed(2)
        train.control <- caret::trainControl(method = "repeatedcv", 
                                      number = 10, repeats = 3)
        data.cv = as.data.frame(cbind(y,x))
        
        model <- caret::train(V1 ~., data = data.cv, method = "lm",
                              trControl = train.control)
        .Random.seed <<- old
        res = summary(lm(y ~ news_variable))
        return(list(error = model$results$RMSE,
                    sign =  res$coefficients[nrow(res$coefficients),1],
                    res = res))
      } else {
        old <- .Random.seed
        set.seed(2)
        error =  bootstrap::bootpred(x,y,200,theta.fit,theta.predict,
                 err.meas=sq.err)[[3]]
        .Random.seed <<- old
        res = summary(lm(y ~ x))
        
        return(list(error = error,
                    sign =  res$coefficients[nrow(res$coefficients),1],
                    res = res))
      }
    } else {
      if(CV){
        old <- .Random.seed
        set.seed(2)
        train.control <- caret::trainControl(method = "repeatedcv", 
                                             number = 10, repeats = 3)
        data.cv = as.data.frame(cbind(y,news_variable))
        
        model <- caret::train(V1 ~., data = data.cv, method = "lm",
                              trControl = train.control)
        .Random.seed <<- old
        res = summary(lm(y ~ news_variable))
        if(dim(res$coefficients)[1] == 1){
          return(list(error = model$results$RMSE,
                      sign =  0))
        } else {
          return(list(error = model$results$RMSE,
                      sign =  res$coefficients[nrow(res$coefficients),1],
                      res = res))
        }
      } else {
        old <- .Random.seed
        set.seed(2)
        error =  bootstrap::bootpred(news_variable,y,200,theta.fit,theta.predict,
                 err.meas=sq.err)[[3]]
        .Random.seed <<- old
        res = summary(lm(y ~ news_variable))
        if(dim(res$coefficients)[1] == 1){
          return(list(error = error,
                      sign =  0))
        } else {
          return(list(error = error,
                      sign =  res$coefficients[nrow(res$coefficients),1],
                      res = res))
        }
      }
    }
  } else {
    if(!is.null(x)){
 
      id = match(rownames(news_variable),as.character(x$date))
      x = x[id,-1]
      x = cbind(x,news_variable)
      x = as.matrix(x)
      
      if(CV){
        old <- .Random.seed
        set.seed(2)
        train.control <- caret::trainControl(method = "repeatedcv", 
                                             number = 10, repeats = 3)
        data.cv = as.data.frame(cbind(y,x))
        
        model <- caret::train(V1 ~., data = data.cv, method = "lm",
                              trControl = train.control)
        .Random.seed <<- old
        res = summary(lm(y ~ news_variable))
        return(list(error = model$results$RMSE,
                    sign =  res$coefficients[nrow(res$coefficients),1],
                    res = res))
      } else {
      
        res = summary(lm(y ~ x))
    
        return(list(error = mean(res$residuals^2),
                    sign =  res$coefficients[nrow(res$coefficients),1],
                    res = res))
      }
    } else {
      
      if(CV){
        old <- .Random.seed
        set.seed(2)
        train.control <- caret::trainControl(method = "repeatedcv", 
                                             number = 10, repeats = 3)
        data.cv = as.data.frame(cbind(y,news_variable))
        
        model <- caret::train(V1 ~., data = data.cv, method = "lm",
                              trControl = train.control)
        .Random.seed <<- old
        res = summary(lm(y ~ news_variable))
        if(dim(res$coefficients)[1] == 1){
          return(list(error = model$results$RMSE,
                      sign =  0,
                      res = res))
        } else {
          return(list(error = model$results$RMSE,
                      sign =  res$coefficients[nrow(res$coefficients),1],
                      res = res))
        }
      } else {
        
        res = summary(lm(y ~ news_variable))
        if(dim(res$coefficients)[1] == 1){
          return(list(error = mean(res$residuals^2),
                      sign =  0,
                      res = res))
        } else {
          return(list(error = mean(res$residuals^2),
                      sign =  res$coefficients[nrow(res$coefficients),1],
                      res = res))
        }
      }
    }
  }
}



penalized_MSEK_card = function(news_variable,
                               y, 
                               x = NULL,
                               lambda_0 = 0.25, 
                               lambda_1 = 0,
                               lambda_2 = 0,
                               zeros = 0,  
                               sign = 1, 
                               omega_select, 
                               bootstrap = FALSE,
                               CV = FALSE){

    res = MSE(news_variable = news_variable, y = y, x =  x, bootstrap, CV)

   test2 = stringr::str_split(rownames(omega_select),pattern = " ")
  #   
   id_which_ngrams = which(sapply(test2, length) > 1)
  #   
  id.penalize = FALSE
  if(length(id_which_ngrams) != 0){
      for(i in 1:length(id_which_ngrams)){
       id.match = match(test2[[id_which_ngrams[i]]],rownames(omega_select))
       
       if(all(is.na(id.match))){
         cond = FALSE
       } else {
         id.match = id.match[!is.na(id.match)]
         cond = sum(colSums(omega_select[c(id.match,id_which_ngrams[i]),]) > 0)
         if(cond >1){
           cond = TRUE
         } else{
           cond = FALSE
         }
       }
       
        id.penalize = id.penalize | cond
     }
  }
  id.penalize = as.numeric(id.penalize)
  
    return(res$error  +
             1 * ifelse(sign > 0, yes =  -1*min(0,sign(res$sign)), no =   1*max(0,sign(res$sign)))+
             lambda_0 * sum(rowSums(aggregate(omega_select, by = list(rownames(omega_select)), FUN = sum)[,-1, drop = FALSE]) > 1) +
             lambda_0 * id.penalize +
             lambda_1*zeros + 
             lambda_2 * as.numeric(max(colSums(omega_select)) - min(colSums(omega_select)) > 2)
    )
}


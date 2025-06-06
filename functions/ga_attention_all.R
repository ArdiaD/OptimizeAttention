ga_attention_all = function(dfm, 
                        y_train,
                        y_valid,
                        y_test,
                        ndim = 3,
                        x = NULL,
                        freq_y = "month",
                        keywords,
                        bootstrap = TRUE,
                        CV = FALSE,
                        maxiter = 100,
                        popSize = 100,
                        monitor,
                        ncore,
                        k,
                        batch_size,
                        lambda_0 = 0.25,
                        lambda_1 = 0,
                        lambda_2 = 0,
                        sentiment = FALSE,
                        crossover = GA::gabin_uCrossover,
                        mutation = GA::gabin_raMutation,
                        suggestion_gamma = NULL,
                        sign = 1,
                        f_transformation = function(x) x){
  
  library(doFuture)
  cl <- parallel::makeCluster(ncore)
  registerDoParallel(cl)
  
  dfm_sub <- dfm_select(dfm, pattern = unlist(keywords))
  
  dfm_sub$frequency_date = floor_date(x = dfm_sub@docvars$date, unit = freq_y)
  
  mat_doc_freq = Matrix.utils::dMcast(dfm_sub@docvars,
                                      formula = frequency_date  ~ docid_)
  

  
  omega = matrix(0,nrow = length(keywords)*ndim, ncol = ndim)
  start = 1
  for(i in 1:ncol(omega)){
    omega[start:(start+length(keywords)-1),i] = 1
    start = start + length(keywords)
  }
 
  rownames(omega) = rep(keywords, ndim)
  n_feat = nrow(omega)
  validation_loss = rep(NA, maxiter)
  train_loss = rep(NA, maxiter)
  test_loss = rep(NA, maxiter)
  old_xx = NA

  
  parallel::clusterExport(varlist = c("conditioning_variable","omega","dfm_sub","mat_doc_freq","batch_size","compute_attention_measure","sentiment","f_transformation"),
                          cl = cl,envir = environment())
  
  loss_function = function(selection, best = FALSE, k) {

    
    batches = split(1:nrow(selection), ceiling(seq_along(1:nrow(selection))/batch_size))
    
    dfm_sub_sub = dfm_sub[,  rownames(omega)[colSums(selection) != 0]]
    omega_sub = omega[colSums(selection) != 0,, drop = FALSE]
    selection_sub = selection[,colSums(selection) != 0]
    keywords_selected = list()
    
    cat(summary(rowSums(selection_sub)),"\n")

    res = foreach(i = 1:length(batches)) %dorng% {
      library(quanteda)
      dfm_sub_sub_sub = dfm_sub[,  rownames(omega)[colSums(selection[batches[[i]],]) != 0]]
      omega_sub_sub = omega[colSums(selection[batches[[i]],]) != 0,, drop = FALSE]
      selection_sub_sub = selection[batches[[i]],colSums(selection[batches[[i]],]) != 0]
      
      feature = compute_attention_measure(dfm = dfm_sub_sub_sub, 
                                          omega = omega_sub_sub,
                                          gamma = selection_sub_sub,
                                          mat_doc_freq = mat_doc_freq,
                                          attention_rule = NULL,
                                          sentiment = sentiment)
      zeros = colSums(feature$unormalized_attention <= 0.001)/nrow(feature$unormalized_attention)
      feature$unormalized_attention = apply(feature$unormalized_attention,MARGIN = 2,f_transformation)
      if(sentiment){
        feature$sentiment = apply(feature$sentiment,MARGIN = 2,f_transformation)
      }
      gc()
      return(list(unormalized_attention = feature$unormalized_attention, sentiment = feature$sentiment, zeros = zeros))
    }
    
    feature = list()
    feature$unormalized_attention = do.call(cbind,lapply(res, FUN = function(x) x$unormalized_attention))
    feature$sentiment = do.call(cbind,lapply(res, FUN = function(x) x$sentiment))
    zeros = do.call(c,lapply(res, FUN = function(x) x$zeros))
    
    res = foreach(i = 1:nrow(selection_sub),.combine = c,.export = c("sentiment","sign","penalized_MSEK_card","bootstrap","MSE","x","y_train","lambda_0","lambda_1","lambda_2")) %dorng% {
   
      xx = diag(selection_sub[i,]) %*% omega_sub
      rownames(xx) = rownames(omega_sub)
      xx = xx[rowSums(xx) != 0,, drop = FALSE]
      if(nrow(xx) == 0){
        res[i] = 100
      } else {
      Cat_0 = sum(colSums(diag(selection_sub[i,]) %*% omega_sub) == 0)
      
      if(!sentiment){
        
        
          res = penalized_MSEK_card(news_variable = feature$unormalized_attention[,i],
                                       x = x,
                                       zeros  =zeros[i],
                                       sign = sign,
                                       y = y_train,
                                       omega_select = xx,
                                       lambda_0 = lambda_0,
                                       lambda_1 = lambda_1,
                                       lambda_2 = lambda_2,
                                       bootstrap = bootstrap)
        } else {
          res = penalized_MSEK_card(news_variable = feature$sentiment[,i],
                                       x = x,
                                       zeros  =zeros[i],
                                       k = k,
                                       sign = sign,
                                       y = y_train,
                                       omega_select = xx,
                                       lambda_0 = lambda_0,
                                       lambda_1 = lambda_1,
                                       lambda_2 = lambda_2,
                                       bootstrap = bootstrap)
        }
      }
      res
    }
    if(best){
  
      xx = diag(selection_sub[which.min(res),]) %*% omega_sub
      rownames(xx) = rownames(omega_sub)
      xx = xx[rowSums(xx) != 0,, drop = FALSE]
      
      diff = base::setdiff(rownames(old_xx), rownames(xx))
      cat(crayon::red(diff),"\n")
      diff2 = base::setdiff(rownames(xx),rownames(old_xx))
      cat(crayon::green(diff2),"\n")
      old_xx <<- xx
      
      if(!sentiment){
       
        train_loss[which(is.na(train_loss))[1]] <<- res[which.min(res)]
        
        validation_loss[which(is.na(validation_loss))[1]] <<-  penalized_MSEK_card(news_variable = feature$unormalized_attention[,which.min(res)],
                                                                                     x = x,
                                                                                     zeros  =zeros[which.min(res)],
                                                                                     sign = sign,
                                                                                     y = y_train,
                                                                                     omega_select = xx,
                                                                                     lambda_0 = lambda_0,
                                                                                     lambda_1 = lambda_1,
                                                                                     lambda_2 = lambda_2,
                                                                                     bootstrap = FALSE,
                                                                                     CV = CV)
          
        test_loss[which(is.na(test_loss))[1]] <<-  penalized_MSEK_card(news_variable = feature$unormalized_attention[,which.min(res)],
                                                                         x = x,
                                                                         zeros  =zeros[which.min(res)],
                                                                         sign = sign,
                                                                         y = y_test,
                                                                         omega_select = xx,
                                                                         lambda_0 = lambda_0,
                                                                         lambda_1 = lambda_1,
                                                                         lambda_2 = lambda_2,
                                                                         bootstrap = FALSE)
        
        print(MSE(news_variable = feature$unormalized_attention[,which.min(res)], y = y_test, x =  x,bootstrap =  FALSE)$res)
        
      } else {
      
        train_loss[which(is.na(train_loss))[1]] <<- res[which.min(res)]
        
        validation_loss[which(is.na(validation_loss))[1]] <<-  penalized_MSEK_card(news_variable = feature$sentiment[,which.min(res)],
                                                                                     x = x,
                                                                                     zeros  =zeros[which.min(res)],
                                                                                     sign = sign,
                                                                                     y = y_train,
                                                                                     omega_select = xx,
                                                                                     lambda_0 = lambda_0,
                                                                                     lambda_1 = lambda_1,
                                                                                     lambda_2 = lambda_2,
                                                                                     bootstrap = FALSE,
                                                                                     CV = CV)
        
        test_loss[which(is.na(train_loss))[1]] <<-  penalized_MSEK_card(news_variable = feature$sentiment[,which.min(res)],
                                                                         x = x,
                                                                         zeros  =zeros[which.min(res)],
                                                                         sign = sign,
                                                                         y = y_test,
                                                                         omega_select = xx,
                                                                         lambda_0 = lambda_0,
                                                                         lambda_1 = lambda_1,
                                                                         lambda_2 = lambda_2,
                                                                         bootstrap = FALSE)
        
        print(MSE(news_variable = feature$sentiment[,which.min(res)], y = y_test, x =  x,bootstrap =  FALSE)$res)
      }
      
     zoo::plot.zoo(zoo::zoo(cbind(train_loss,validation_loss,test_loss),order.by = 1:length(train_loss)),plot.type = "single",col = c("blue","green","red"))


     print(xx)

     cat(sum(selection_sub[which.min(res),]),"\n")
    }
    return(res)
  }
  

    ids = sapply(1:popSize,FUN = function(x) {sample(1:(n_feat/ndim),
                                                     size = k/ndim,
                                                     prob = rep(1/(n_feat/ndim), times = (n_feat/ndim)),
                                                     replace = FALSE)})
    if(!is.matrix(ids)){
      ids = t(as.matrix(ids))
    }
                 
    suggestions =     t(apply(ids,MARGIN = 2,FUN = function(x) {
                                                                xx = rep(0,times = (n_feat/ndim))
                                                                xx[x] = 1
                                                                return(xx)}))
    if(ndim != 1){
      for(i in 2:ndim){
        
        ids = sapply(1:popSize,FUN = function(x) {sample(1:(n_feat/ndim),
                                                         size = k/ndim,
                                                         prob = rep(1/(n_feat/ndim), times = (n_feat/ndim)),
                                                         replace = FALSE)})
        
        if(!is.matrix(ids)){
          ids = t(as.matrix(ids))
        }
        
        suggestions2 =     t(apply(ids,MARGIN = 2,FUN = function(x) {
                                                                      xx = rep(0,times = (n_feat/ndim))
                                                                      xx[x] = 1
                                                                      return(xx)}))
        
        suggestions  = cbind(suggestions,suggestions2)
      }
    }                                                                   

    colnames(suggestions) = rownames(omega)
    if(!is.null(suggestion_gamma)){
      suggestions[1:(popSize),] = t(replicate((popSize), suggestion_gamma[1,])) 
    }
    

    GA <- ga_custom(type = "binary",
                    fitness = function(x, best = FALSE, k) {-loss_function(x, best, k)}, 
                    nBits = n_feat,  
                    maxiter = maxiter, 
                    popSize = popSize,
                    parallel = FALSE,
                    crossover = crossover,
                    selection = GA::gabin_tourSelection,
                    mutation = mutation,
                    suggestions = suggestions,
                    k = k,
                    ndim = ndim,
                    keepBest = TRUE,
                    monitor = monitor,
                    batch = TRUE,
                    pop_names =  rownames(omega))
    
 
  
  if(is.matrix(GA@solution)){
    sol = GA@solution[1,]
  } else {
    sol = GA@solution
  }
  stopCluster(cl)
  optimal = list(sol = GA@bestSol,
                 train_loss = train_loss,
                 validation_loss = validation_loss)
  
  return(optimal)
}


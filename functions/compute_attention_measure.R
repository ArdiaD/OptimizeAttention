
compute_attention_measure = function(dfm,
                                     omega,
                                     gamma,
                                     freq = "month",
                                     mat_doc_freq = NULL,
                                     attention_rule = NULL,
                                     sentiment = FALSE){

  if(is.null(mat_doc_freq)){
    dfm$frequency_date = lubridate::floor_date(x = dfm@docvars$date, unit = freq)
    
    mat_doc_freq = Matrix.utils::dMcast(dfm@docvars,
                                        formula = frequency_date  ~ docid_)
  }
  library(Matrix)
  n = nrow(gamma)
  gamma_diag = apply(gamma,MARGIN = 1,simplify = FALSE, FUN = function(x) array(Matrix::Diagonal(x=x),dim = c(ncol(gamma),ncol(gamma),1)))
  gamma_diag = array(unlist(gamma_diag, use.names = FALSE),dim = c(ncol(gamma),ncol(gamma),nrow(gamma)))
  gamma_diag <- Matrix::Matrix(gamma_diag, nrow = ncol(gamma),sparse = TRUE)

  gamma_diag = as(gamma_diag, "dgCMatrix")
  btranform = Matrix::bdiag(rep(list(rep(1,ncol(omega))),n))
  
  test = Matrix::bdiag(rep(list(omega),n))
  
  test2 = gamma_diag %*% test
  
  test3 = dfm %*% test2
  test3@x = rep(1,length(test3@x))
  
  if(nrow(gamma) == 1) {
    ndim = matrix(1/as.vector((colSums(test2) >0)%*%btranform),1,1)
  } else {
    ndim = Matrix::Matrix(diag(1/as.vector((colSums(test2) >0)%*%btranform)),sparse = TRUE)
  }
  test5 =  test3 %*% btranform %*% ndim
  test5 = test5 == 1
  test6 =  mat_doc_freq %*% test5

  
  test7 = Matrix::Matrix(diag(1/rowSums(mat_doc_freq)),sparse = TRUE) %*% test6
  
  if(sentiment){

    sentiment = Matrix::Diagonal(x = as.vector(dfm@docvars$sentiment)) %*% test5
    
    sentiment =  (mat_doc_freq %*% sentiment)
    sentiment[is.na(sentiment)] = 0
    unormalized_attention = as.matrix(test7)
   
    rownames(unormalized_attention) =   rownames(sentiment) =  rownames(mat_doc_freq)
    return(list(unormalized_attention = unormalized_attention,
                sentiment = sentiment,
                sel = test5))
  } else {
    unormalized_attention = as.matrix(test7)
    rownames(unormalized_attention) =   rownames(mat_doc_freq)
    return(list(unormalized_attention = unormalized_attention,
                sentiment = 0,
                sel = test5))
  }

}

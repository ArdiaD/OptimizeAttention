glove_to_text2vec<- function(file_name) {
  glove <- scan(file = file_name, what="", sep="\n")
  
  vals <- vector(mode = "list", length(glove))
  names <- character(length(glove))
  
  for(i in 1:length(glove)) {
    if(i %% 1000 == 0) {print(i)}
    this_vec <- glove[i]
    this_vec_unlisted <- unlist(strsplit(this_vec, " "))
    this_vec_values <- as.numeric(this_vec_unlisted[-1])
    this_vec_name <- this_vec_unlisted[1]
    
    vals[[i]] <- this_vec_values
    names[[i]] <- this_vec_name
  }
  
  WordVectorMatrix <- data.frame(vals)
  names(WordVectorMatrix) <- names
  
  return (WordVectorMatrix)
}


near_terms <- function(keywords,
                       wordvectors = wordvectors, 
                       keep=ncol(wordvectors),
                       keep.score = FALSE,
                       keep.original = TRUE) {
  lapply(keywords, 
         FUN = f_near_terms, 
         keep = keep, 
         wordvectors = wordvectors,
         keep.score = keep.score,
         keep.original)
}

f_loc_wv <- function(term, wordvectors) {
  if(term %in% colnames(wordvectors)){
    wordvectors[,term]
  } else {
    rep(NA,times = nrow(wordvectors))
  }
}


f_near_terms <- function(terms,
                         wordvectors,
                         keep,
                         keep.score,
                         keep.original){

  res = find_sim_wvs(this_wv = f_avg_wv(wordvectors = wordvectors,
                                        terms = terms),
                     all_wvs = wordvectors)
  res = data.frame(keyword = names(res), similarity = res)
  rownames(res) = NULL

  if(keep < 1 ){
    res = res[res$similarity > keep,]
  } else {
    res = head(res, keep)
  }
  if(keep.original){
    res = rbind(data.frame(keyword = terms,similarity = 1),res)
    res = res[!duplicated(res[,1]),]
  }
  if(keep.score){
    return(res)
  } else {
    return(res[,1])
  }
  return(res)
}

f_avg_wv <- function(wordvectors, terms, aggregate_n_grams = FALSE) {
  
  if(aggregate_n_grams){
    terms = unlist(stringr::str_split(terms, pattern =  " "))
  }
  current_wv = sapply(terms, FUN = f_loc_wv, wordvectors =  wordvectors)
  if(dim(current_wv)[2] == 1){
    current_wv
  } else {
    rowMeans(current_wv,na.rm = TRUE)
  }
}

find_sim_wvs <- function(this_wv, all_wvs) {
  this_wv_mat <- matrix(this_wv, ncol=length(this_wv), nrow=1)
  all_wvs_mat <- as.matrix(all_wvs)
  
  if(dim(this_wv_mat)[[2]] != dim(all_wvs_mat)[[2]]) {

    all_wvs_mat <- t(all_wvs_mat)
  }
  
  cos_sim = f_compute_cosine_similarities(x = all_wvs_mat,
                                          y = this_wv_mat, 
                                          norm = TRUE)
  
  sorted_cos_sim <-  cos_sim[order(cos_sim[,1], decreasing = T),1]

  return(sorted_cos_sim)
  
}


f_compute_cosine_similarities = function (x, y, normalize = TRUE) {
  if(normalize){
    x = text2vec::normalize(x, norm = "l2")
    y = text2vec::normalize(y, norm = "l2")
    res = tcrossprod(x, y)
  } else {
    res = tcrossprod(x, y)
  }
  res
}

solve_duplicated = function(keywords){
  keywords_expanded_add = list()
  for(i in 1:length(keywords)){
    
    keywords_expanded_add[[i]] = cbind(keywords[[i]],label = names(keywords)[i])
  }
  keywords_expanded_add = rbindlist(keywords_expanded_add)
  keywords_expanded_add = keywords_expanded_add[order(keywords_expanded_add$similarity,decreasing = TRUE),]
  keywords_expanded_add = keywords_expanded_add[!duplicated(keywords_expanded_add$keyword),]
  keywords_expanded = split(keywords_expanded_add$keyword,factor(keywords_expanded_add$label,levels = names(keywords)))
  return(keywords_expanded)
}

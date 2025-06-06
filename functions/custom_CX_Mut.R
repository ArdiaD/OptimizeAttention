custom_crossover = function(object, parent, k = NULL, ndim) {
  
  parents = object@population[parent, , drop = FALSE]

  if(all(parents[1,] == parents[2,])){
    return(list(children = parents, fitness = rep(NA, 2)))
  }
   
  sector1 = sample(1:ndim, size = 1)
  sector2 = sample(1:ndim, size = 1)
  
  end1 = sector1*ncol(parents)/ndim
  start1 = sector1*ncol(parents)/ndim - ncol(parents)/ndim +1
  
  end2 = sector2*ncol(parents)/ndim
  start2 = sector2*ncol(parents)/ndim - ncol(parents)/ndim +1
  
  child1 = parents[1,]
  child2 = parents[2,]
  id_1 = which(parents[1,start1:end1] != 0)
  id_2 = which(parents[2,start2:end2] != 0)
  id_1_filt = id_1
  id_2_filt = id_2
  
  id_1_filt = id_1[!(id_1 %in% id_2)]
  id_2_filt = id_2[!(id_2 %in% id_1)]
  
  if(length(id_1_filt) == 0){
    return(list(children = parents, fitness = rep(NA, 2)))
  }
  if(length(id_2_filt) == 0){
    return(list(children = parents, fitness = rep(NA, 2)))
  }
  

  if(length(id_1_filt) > 1){
     id_1_filt = sample(id_1_filt,size = 1)
  } 
  if(length(id_2_filt) > 1){
     id_2_filt = sample(id_2_filt,size = 1)
  }
  
  child1[start1-1+id_1_filt] = 0
  child1[start1-1+id_2_filt] = 1

  child2[start2-1+id_2_filt] = 0
  child2[start2-1+id_1_filt] = 1
  
  children = rbind(child1, child2)
 

  return(list(children = children, fitness = rep(NA, 2)))
  
}

custom_mutation = function (current_gene, k, pmutation, n.mutation, ndim) {

   mutate <- parent <- current_gene
   
   for(z in 1:n.mutation){
     ex = NA
     for(i in 1:ndim){
       section = i
       end = section*length(parent)/ndim
       start =  section*length(parent)/ndim - length(parent)/ndim + 1
       ex[i] = sum(mutate[start:end])
     }
     
     section = sample(1:ndim, size = 1)
     
     end = section*length(parent)/ndim
     start =  section*length(parent)/ndim - length(parent)/ndim + 1
     
     id = which(mutate[start:end] == 1)
     mutate[start-1+id] = 0
     j <- sample(start:end, size = 1)
     mutate[j] = 1
   }
 
  return(mutate)
}



switch_mutation = function(current_gene, k, pmutation, n.mutation = k, ndim){

  if(ndim == 1){
    return(current_gene)
  }
  current_gene = current_gene
  n <- length(current_gene)
  mutate = current_gene

  section = sample(1:ndim,size = 1)
  
  end = (1:ndim)*length(current_gene)/ndim

  if(length(end) != 1){
    start = c(1,1:(ndim-1)*length(current_gene)/ndim+1)
  } else {
    start = 1
  }

  ids = list()
  for(i in 1:ndim){
    ids[[i]] = which(mutate[start[i]:end[i]] == 1)
  }
  l.section = sapply(ids,FUN = length)
  if(all(l.section <= 1)){
    return(mutate)
  } else {
    p.section = which(l.section>1)
    if(length(p.section) == 1){
      section = p.section
    } else {
      section = sample(p.section,  size = 1) 
    }
  }

  id = sample(ids[[section]], size = 1)
  
  mutate[id+start[section]-1] = 0
  
  if(ndim == 2){
    section2 = (1:ndim)[-section]
  } else {
    section2 =  sample((1:ndim)[-section],size = 1)
  }
  
  if(mutate[id+start[section2]-1] == 1){
    return(current_gene)
  }
  mutate[id+start[section2]-1] = 1

  return(mutate)
}


evolve_mutation = function(current_gene, k, pmutation, n.mutation = k, pop_names, ndim){

  current_gene = current_gene
  n <- length(current_gene)
  mutate = current_gene

  id = which(current_gene == 1)
  ngrams = vector()
  k = 0
  while(length(ngrams) == 0){
    if(length(id) > 1){
      id_sel = sample(id,size = 1)
    } else {
      id_sel = id
    }
    sector = sort(rep(1:ndim,times = n/ndim))[id_sel]
    selected = pop_names[id_sel]
    
    ngrams = unique(pop_names[which(stringr::str_detect(pop_names,pattern = paste0(" ",selected,"$")) |  stringr::str_detect(pop_names,pattern = paste0("^",selected," ")))])
    k = k+1
    if(k > 10){
      break
    }
  }
  if(length(ngrams) == 0){
    return(current_gene)
  }
  
  mutate[id_sel] = 0
  if(length(ngrams) > 1){
    ngrams =  sample(ngrams, size = 1)
  }
  id2 = which(pop_names == ngrams)
  id3 = id2[sector]
  mutate[id3] = 1
  

  return(mutate)  
}



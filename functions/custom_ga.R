##############################################################################
#                                                                            #
#                        GENETIC ALGORITHMS in R                             #
#                                                                            #
##############################################################################

ga_custom <- function(type = c("binary", "real-valued", "permutation"), 
               fitness, ...,
               lower, upper, nBits,
               population = gaControl(type)$population,
               selection = gaControl(type)$selection,
               crossover = gaControl(type)$crossover, 
               mutation = gaControl(type)$mutation,
               popSize = 50, 
               pcrossover = 0.80, 
               pmutation = 0.1, 
               elitism = base::max(1, round(popSize*0.05)), 
               updatePop = FALSE,
               postFitness = NULL,
               maxiter = 100,
               ndim,
               run = maxiter,
               maxFitness = Inf,
               names = NULL,
               suggestions = NULL,
               optim = FALSE,
               optimArgs = list(method = "L-BFGS-B", 
                                poptim = 0.05,
                                pressel = 0.5,
                                control = list(fnscale = -1, maxit = 100)),
               keepBest = FALSE,
               parallel = FALSE,
               batch = FALSE,
               monitor = if(interactive()) gaMonitor else FALSE,
               seed = NULL,
               k = NULL,
               pop_names = NULL) {
  
  call <- match.call()
  
  type <- match.arg(type, choices = eval(formals(ga)$type))
  
  if(!is.function(population)) population <- get(population)
  if(!is.function(selection))  selection  <- get(selection)
  if(!is.function(crossover))  crossover  <- get(crossover)
  if(!is.function(mutation))   mutation   <- get(mutation)
  
  if(missing(fitness)){ 
    stop("A fitness function must be provided") 
    }
  if(!is.function(fitness)) { 
    stop("A fitness function must be provided")
    }
  if(popSize < 10) { 
    warning("The population size is less than 10.") 
    }
  if(maxiter < 1) { 
    stop("The maximum number of iterations must be at least 1.")
    }
  if(elitism > popSize) { 
    stop("The elitism cannot be larger that population size.")
  }
  elitism <- as.integer(elitism)
  if(pcrossover < 0 | pcrossover > 1){
    stop("Probability of crossover must be between 0 and 1.") 
  }
  if(is.numeric(pmutation)){ 
    if(pmutation < 0 | pmutation > 1) {
      stop("If numeric probability of mutation must be between 0 and 1.")
    } else if(!is.function(population)) {
        stop("pmutation must be a numeric value in (0,1) or a function.")
    }
  }
  
  # check for min and max arguments instead of lower and upper
  callArgs <- list(...)
  if(any("min" %in% names(callArgs))) {
    lower <- callArgs$min
    callArgs$min <- NULL
    warning("'min' arg is deprecated. Use 'lower' instead.")
  }
  if(any("max" %in% names(callArgs))) {
    upper <- callArgs$max
    callArgs$max <- NULL
    warning("'max' arg is deprecated. Use 'upper' instead.")
  }
  
  if(missing(lower) & missing(upper) & missing(nBits)){ 
    stop("A lower and upper range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!") 
  }
  
  # check GA search type 
  switch(type, 
         "binary"  = { nBits <- as.vector(nBits)[1]
         lower <- upper <- NA
         nvars <- nBits 
         if(is.null(names))
           names <- paste0("x", 1:nvars)
         },
         "real-valued" = { lnames <- names(lower)
         unames <- names(upper)
         lower <- as.vector(lower)
         upper <- as.vector(upper)
         nBits <- NA
         if(length(lower) != length(upper))
           stop("lower and upper must be vector of the same length!")
         nvars <- length(upper)
         if(is.null(names) & !is.null(lnames))
           names <- lnames
         if(is.null(names) & !is.null(unames))
           names <- unames
         if(is.null(names))
           names <- paste0("x", 1:nvars)
         },
         "permutation" = { lower <- as.vector(lower)[1]
         upper <- as.vector(upper)[1]
         nBits <- NA
         nvars <- length(seq.int(lower,upper))
         if(is.null(names))
           names <- paste0("x", 1:nvars)
         }
  )
 
  # check suggestions
  if(is.null(suggestions)){ 
    suggestions <- matrix(nrow = 0, ncol = nvars)
  }else { 
    if(is.vector(suggestions)) {
      if(nvars > 1) {
        suggestions <- matrix(suggestions, nrow = 1)
      } else  {       
        suggestions <- matrix(suggestions, ncol = 1) 
      }
    } else{ 
      suggestions <- as.matrix(suggestions) 
    }
    if(nvars != ncol(suggestions)){
      stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
    }
  }
  
  # check monitor arg
  if(is.logical(monitor)){ 
    if(monitor){ monitor <- gaMonitor
    }
  }
  if(is.null(monitor)){
    monitor <- FALSE
    } 
  
  # if optim merge provided and default args for optim()
  if(optim) { # merge default and provided parameters
    optimArgs.default <- eval(formals(ga)$optimArgs)
    optimArgs.default$control[names(optimArgs$control)] <- optimArgs$control
    optimArgs$control <- NULL
    optimArgs.default[names(optimArgs)] <- optimArgs
    optimArgs <- optimArgs.default; rm(optimArgs.default)
    if(any(optimArgs$method == c("L-BFGS-B", "Brent"))){ 
      optimArgs$lower <- lower
       optimArgs$upper <- upper 
    } else { 
      optimArgs$lower <- -Inf
      optimArgs$upper <- Inf }
      optimArgs$poptim <- min(max(0, optimArgs$poptim), 1)
      optimArgs$pressel <- min(max(0, optimArgs$pressel), 1)
      optimArgs$control$maxit <- as.integer(optimArgs$control$maxit)
    # ensure that optim maximise the fitness
      if(is.null(optimArgs$control$fnscale)){
         optimArgs$control$fnscale <- -1
      }
    if(optimArgs$control$fnscale > 0){
      optimArgs$control$fnscale <- -1*optimArgs$control$fnscale
    }
  }
  
  # set seed for reproducibility  
  if(!is.null(seed)) set.seed(seed)
  i. <- NULL # dummy to trick R CMD check 

  fitnessSummary <- matrix(as.double(NA), nrow = maxiter, ncol = 6)
  colnames(fitnessSummary) <- names(gaSummary(rnorm(10)))
  
  bestSol <- if(keepBest) {
    vector(mode = "list", length = maxiter)
  } else {
    list()
  }       
  
  Fitness <- rep(NA, popSize)
  
  object <- new("ga", 
                call = call, 
                type = type,
                lower = lower, 
                upper = upper, 
                nBits = nBits, 
                names = if(is.null(names)) character() else names,
                popSize = popSize,
                iter = 0, 
                run = 1, 
                maxiter = maxiter,
                suggestions = suggestions,
                population = matrix(), 
                elitism = elitism, 
                pcrossover = pcrossover, 
                pmutation = if(is.numeric(pmutation)) pmutation else NA,
                optim = optim,
                fitness = Fitness, 
                summary = fitnessSummary,
                bestSol = bestSol)
  
  if(maxiter == 0){
    return(object)
  }
  
  # generate beginning population
  Pop <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  ng <- min(nrow(suggestions), popSize)
  if(ng > 0) { 
    Pop[1:ng,] <- suggestions
  }
  # fill the rest with a random population
  if(popSize > ng){ 
    Pop[(ng+1):popSize,] <- population(object)[1:(popSize-ng),] 
  }
  object@population <- Pop

  
  for(i in 2:popSize) {  
    while(sum(Pop[i,]) < k){
      
      j <- sample(1:ncol(Pop), size = 1)
      Pop[i,j] = 1
    }
  }
  
  # start iterations
  for(iter in seq_len(maxiter)) {
    
    object@iter <- iter

    end = 1:ndim*ncol(Pop)/ndim
    if(length(end) != 1){
      start = c(1,1:(ndim-1)*ncol(Pop)/ndim+1)
    } else {
      start = 1
    }
    
    for(i in 1:nrow(Pop)){
      min_word = c()
      for(j in 1:length(start)){
        min_word[j]=  min(which(Pop[i,start[j]:end[j]]==1))
      }
      ord_dim = order(min_word)
      Pop_temp = Pop[i,]
      for(j in 1:length(start)){
        Pop[i,start[j]:end[j]] = Pop_temp[start[ord_dim[j]]:end[ord_dim[j]]]
      }
    }
  
    fit = rep(NA,times = nrow(Pop))
    if(batch) {
      fit <- DescTools::StripAttr(do.call(fitness, c(list(Pop, TRUE, k))))
      Fitness <- fit
    } 


    
    # update object
    object@population <- Pop
    object@fitness <- Fitness
    
    
    if(keepBest) { 
      object@bestSol[[iter]] <- unique(Pop[Fitness == max(Fitness, na.rm = TRUE),, drop=FALSE]) 
    }
    
    # update iterations summary
    fitnessSummary[iter,] <- gaSummary(object@fitness)
    object@summary <- fitnessSummary
    
    if(is.function(monitor)) { 
      monitor(object) 
    }
    
    # check stopping criteria
    if(iter > 1){
      object@run <- garun(fitnessSummary[seq(iter),1])
    }
    if(object@run >= run) {
      break
    } 
    if(max(Fitness, na.rm = TRUE) >= maxFitness) {
      break
    }
    if(object@iter == maxiter){
      break
    } 
    
    Pop =   object@population 
    Fitness = object@fitness
    ord <- order(Fitness, decreasing = TRUE)
    

    object@population <- Pop
    object@fitness <- Fitness
   
    PopSorted <- Pop[ord,,drop=FALSE]
    FitnessSorted <- Fitness[ord]
    
    pos = matrix(NA, ncol = 2, nrow = ndim)
    
    # selection
    if(is.function(selection)){ 
      sel <- selection(object)
      Pop <- sel$population
      Fitness <- sel$fitness
    } else { 
      sel <- sample(1:popSize, size = popSize, replace = TRUE)
      Pop <- object@population[sel,]
      Fitness <- object@fitness[sel]
    }
    
    object@population <- Pop
    object@fitness <- Fitness
    orig_Fitness = Fitness
  
    
    for(i in seq_len(popSize)) { 
      if(0.05 > runif(1)){
        Mutation <- switch_mutation(Pop[i,], k, pmutation, n.mutation = 1,  ndim = ndim)
        Pop[i,] <- Mutation
        Fitness[i] <- NA
        n = length(Mutation)
      }
    }
    
    pcrossover = 1-((sum(duplicated(Pop))/nrow(Pop)))
 
    
    # crossover

    ##################################################

    cat("crossover prob = ", pcrossover, "\n")
    object@population <- Pop
    object@fitness <- Fitness
    
    nmating <- floor(popSize/2)
    mating <- matrix(sample(1:(2*nmating), size = (2*nmating)), ncol = 2)
    
    for(i in seq_len(nmating)){ 
      if(pcrossover > runif(1)){ 
        parents <- mating[i,]
        Crossover <- crossover(object, parents, k, ndim  = ndim)
        Pop[parents,] <- Crossover$children
        Fitness[parents] <- Crossover$fitness
      }
    }     
    
    
    object@population <- Pop
    object@fitness <- Fitness

   
 
    for(i in seq_len(popSize)) { 
      if(0.05 > runif(1)){
         Mutation <- evolve_mutation(Pop[i,], k, pmutation, n.mutation = 1, pop_names = pop_names, ndim = ndim)
         Pop[i,] <- Mutation
          Fitness[i] <- NA
          n = length(Mutation)
        }
      }
    
    pmutation =  (sum(duplicated(Pop))/nrow(Pop))
    cat("Mutation prob = ", pmutation, "\n")
    
    for(i in seq_len(popSize)) { 
      if(pmutation > runif(1)){

        Mutation <- custom_mutation(Pop[i,], k, pmutation, n.mutation = sample(1:ndim,size = 1,prob = c(ndim:1/sum(1:ndim))), ndim = ndim)

        n = length(Mutation)
        
        Pop[i,] <- Mutation
        Fitness[i] <- NA
      }
    } 
    
    
    id_dup = which(duplicated(Pop))
    for(i in 1:length(id_dup)) { 
        
        Mutation <- custom_mutation(Pop[id_dup[i],], k, pmutation, n.mutation = sample(1:ndim,size = 1,prob = c(ndim:1/sum(1:ndim))), ndim = ndim)
        
        n = length(Mutation)
        
        Pop[id_dup[i],] <- Mutation
        Fitness[id_dup[i]] <- NA
      
    } 

    old_pop =  object@population
    old_fit =  object@fitness 
    
    ##########################################################################
    
    
    for(i in 1:popSize) {  
      while(sum(Pop[i,]) < k){
        
        j <- sample(1:n, size = 1)
        Pop[i,j] = 1
      }
    }

    # elitism
    if(elitism > 0) {
      ord <- order(object@fitness, na.last = TRUE)
      u <- which(!duplicated(PopSorted, margin = 1))
      Pop[ord[1:elitism],] <- PopSorted[u[1:elitism],]
      Fitness[ord[1:elitism]] <- FitnessSorted[u[1:elitism]]
      object@population <- Pop
      object@fitness <- Fitness
    } 

    aa = colSums(Pop)
    names(aa) = pop_names
    aa = aa[aa != 0]
    print(sort(aa,decreasing = TRUE)[1:20])
  }

    # if(is.function(monitor)) 
  #   { cat("\n"); flush.console() }
  
  # in case of premature convergence remove NAs from summary 
  # fitness evalutations
  object@summary <- na.exclude(object@summary)
  attr(object@summary, "na.action") <- NULL
  
  # get solution(s)
  object@fitnessValue <- max(object@fitness, na.rm = TRUE)
  valueAt <- which(object@fitness == object@fitnessValue)
  solution <- object@population[valueAt,,drop=FALSE]
  if(nrow(solution) > 1)
  { # find unique solutions to precision given by default tolerance
    eps <- gaControl("eps")
    solution <- unique(round(solution/eps)*eps, margin = 1)
  }
  colnames(solution) <- parNames(object)
  object@solution <- solution
  if(keepBest)
    object@bestSol <- object@bestSol[!sapply(object@bestSol, is.null)] 

  # return an object of class 'ga'
  return(object)
}

setClassUnion("numericOrNA", members = c("numeric", "logical"))

setClass(Class = "ga", 
         representation(call = "language",
                        type = "character",
                        lower = "numericOrNA", 
                        upper = "numericOrNA", 
                        nBits = "numericOrNA", 
                        names = "character",
                        popSize = "numeric",
                        iter = "numeric", 
                        run = "numeric", 
                        maxiter = "numeric",
                        suggestions = "matrix",
                        population = "matrix",
                        elitism = "numeric", 
                        pcrossover = "numeric", 
                        pmutation = "numericOrNA",
                        optim = "logical",
                        fitness = "numericOrNA",
                        summary = "matrix",
                        bestSol = "list",
                        fitnessValue = "numeric",
                        solution = "matrix"
         ),
         package = "GA" 
) 

setMethod("print", "ga", function(x, ...) str(x))

setMethod("show", "ga",
          function(object)
          { cat("An object of class \"ga\"\n")
            cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
            cat("Available slots:\n")
            print(slotNames(object))
          }) 

summary.ga <- function(object, ...)
{
  nvars <- ncol(object@population)
  varnames <- parNames(object)
  domain <- NULL
  if(object@type == "real-valued")
  { domain <- rbind(object@lower, object@upper)
  rownames(domain) <- c("lower", "upper")
  if(ncol(domain) == nvars) 
    colnames(domain) <- varnames
  }
  suggestions <- NULL
  if(nrow(object@suggestions) > 0) 
  { suggestions <- object@suggestions
  dimnames(suggestions) <- list(1:nrow(suggestions), varnames) 
  }
  
  out <- list(type = object@type,
              popSize = object@popSize,
              maxiter = object@maxiter,
              elitism = object@elitism,
              pcrossover = object@pcrossover,
              pmutation = object@pmutation,
              domain = domain,
              suggestions = suggestions,
              iter = object@iter,
              fitness = object@fitnessValue,
              solution = object@solution)
  class(out) <- "summary.ga"
  return(out)
}

setMethod("summary", "ga", summary.ga)

print.summary.ga <- function(x, digits = getOption("digits"), ...)
{
  dotargs <- list(...)
  if(is.null(dotargs$head)) dotargs$head <- 10
  if(is.null(dotargs$tail)) dotargs$tail <- 2
  if(is.null(dotargs$chead)) dotargs$chead <- 10
  if(is.null(dotargs$ctail)) dotargs$ctail <- 2
  
  cat(cli::rule(left = crayon::bold("Genetic Algorithm"), 
                width = min(getOption("width"),40)), "\n\n")
  # cat("+-----------------------------------+\n")
  # cat("|         Genetic Algorithm         |\n")
  # cat("+-----------------------------------+\n\n")
  
  cat("GA settings: \n")
  cat(paste("Type                  = ", x$type, "\n"))
  cat(paste("Population size       = ", x$popSize, "\n"))
  cat(paste("Number of generations = ", x$maxiter, "\n"))
  cat(paste("Elitism               = ", x$elitism, "\n"))
  cat(paste("Crossover probability = ", format(x$pcrossover, digits = digits), "\n"))
  cat(paste("Mutation probability  = ", format(x$pmutation, digits = digits), "\n"))
  #
  if(x$type == "real-valued")
  { cat(paste("Search domain = \n"))
    do.call(".printShortMatrix", 
            c(list(x$domain, digits = digits), 
              dotargs[c("head", "tail", "chead", "ctail")]))
  }
  #
  if(!is.null(x$suggestions))
  { cat(paste("Suggestions =", "\n"))
    do.call(".printShortMatrix", 
            c(list(x$suggestions, digits = digits), 
              dotargs[c("head", "tail", "chead", "ctail")]))
  }
  #
  cat("\nGA results: \n")
  cat(paste("Iterations             =", format(x$iter, digits = digits), "\n"))
  cat(paste("Fitness function value =", format(x$fitness, digits = digits), "\n"))
  if(nrow(x$solution) > 1) 
  { cat(paste("Solutions = \n")) }
  else
  { cat(paste("Solution = \n")) }
  do.call(".printShortMatrix", 
          c(list(x$solution, digits = digits), 
            dotargs[c("head", "tail", "chead", "ctail")]))
  #
  invisible()
}


plot.ga <- function(x, y, ylim, cex.points = 0.7,
                    col = c("green3", "dodgerblue3", adjustcolor("green3", alpha.f = 0.1)),
                    pch = c(16, 1), lty = c(1,2), legend = TRUE,
                    grid = graphics:::grid, ...)
{
  object <- x  # Argh.  Really want to use 'object' anyway
  is.final <- !(any(is.na(object@summary[,1])))
  iters <- if(is.final) 1:object@iter else 1:object@maxiter
  summary <- object@summary
  if(missing(ylim)) 
  { ylim <- c(max(apply(summary[,c(2,4)], 2, 
                        function(x) min(range(x, na.rm = TRUE, finite = TRUE)))),
              max(range(summary[,1], na.rm = TRUE, finite = TRUE))) 
  }
  
  plot(iters, summary[,1], type = "n", ylim = ylim, 
       xlab = "Generation", ylab = "Fitness value", ...)
  if(is.final & is.function(grid)) 
  { grid(equilogs=FALSE) }
  points(iters, summary[,1], type = ifelse(is.final, "o", "p"),
         pch = pch[1], lty = lty[1], col = col[1], cex = cex.points)
  points(iters, summary[,2], type = ifelse(is.final, "o", "p"),
         pch = pch[2], lty = lty[2], col = col[2], cex = cex.points)
  if(is.final)
  { polygon(c(iters, rev(iters)), 
            c(summary[,4], rev(summary[,1])), 
            border = FALSE, col = col[3]) }
  else
  { title(paste("Iteration", object@iter), font.main = 1) }
  if(is.final & legend)
  { inc <- !is.na(col)
  legend("bottomright", 
         legend = c("Best", "Mean", "Median")[inc], 
         col = col[inc], pch = c(pch,NA)[inc], 
         lty = c(lty,1)[inc], lwd = c(1,1,10)[inc], 
         pt.cex = c(rep(cex.points,2), 2)[inc], 
         inset = 0.02) }
  
  out <- data.frame(iter = iters, summary)
  invisible(out)
}

setMethod("plot", "ga", plot.ga)

setGeneric(name = "parNames", 
           def = function(object, ...) { standardGeneric("parNames") }
)

setMethod("parNames", "ga",
          function(object, ...)
          { 
            names <- object@names
            nvars <- ncol(object@population)
            if(length(names) == 0)
            { names <- paste("x", 1:nvars, sep = "") }
            return(names)
          })

gaSummary <- function(x, ...)
{
  # compute summary for each step
  x <- na.exclude(as.vector(x))
  q <- fivenum(x)
  c(max = q[5], mean = mean(x), q3 = q[4], median = q[3], q1 = q[2], min = q[1])
}

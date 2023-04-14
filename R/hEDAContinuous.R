#'hybrid Estimation of Distribution Algorithm (hEDA) for continuous strata
#'@export
hEDAContinuous<-function (frame, err, suggestions =NULL,
                          Temp=0.01,initialStr, decrement_constant=0.95, end_time =140,
                          jsize=10,length_of_markov_chain =5,
                          SAArun=TRUE,SAAiters=1000,
                          popSize = 200, iters = 100, mutationChance = NA, elitism = NA,
                          addStrataFactor=0.1, EDAfreq=EDAfreq,
                          verbose = FALSE, dominio=dominio,minnumstrat=2,kmax_percent=0.025,ProbNewStratum=0.0001,
                          strcens=FALSE,writeFiles=FALSE, showPlot=TRUE, minTemp = 0.000005, realAllocation=TRUE){

  # stringMin <- rep(1, nrow(stra))
  # stringMax <- rep(initialStr, nrow(stra))
  require("SamplingStrata")

  stra <- buildStrataDF(frame,progress=FALSE, verbose=FALSE)
  nvar<-length(grep("CV",names(err)))
  stringMin <- rep(1, nrow(stra))
  stringMax <- rep(initialStr, nrow(stra))

  vars = nrow(stra)
  if (is.na(mutationChance)) {
    mutationChance = 1/(vars + 1)
  }
  if (is.na(elitism)) {
    elitismR = floor(popSize/5)
  }else {
    elitismR = floor(popSize*elitism)
  }
  cv_values<-c(as.numeric(err[,2:(nvar+1)]))
  if (verbose)
  {cat("Testing the sanity of parameters...\n")}
  if (length(stringMin) != length(stringMax)) {
    stop("The vectors stringMin and stringMax must be of equal length.")
  }
  if (popSize < 5) {
    stop("The population size must be at least 5.")
  }
  if (iters < 1) {
    stop("The number of iterations must be at least 1.")
  }
  if (!(elitismR < popSize)) {
    stop("The population size must be greater than the elitism.")
  }

  model=NULL


 # frame<-dataset
 # dataset<-frame

  nvar<-length(grep("CV",names(err)))
  evaluateRcpp<-function(sugg){

    newstr<-aggrStrata_RcppOpen(stra, nvar, sugg, censiti,
                       dominio)
    newstr<-as.data.frame(newstr)
    res <- sum(unlist(bethel_alfa(newstr, err[,2:(nvar+1)],
                                  minnumstrat=minnumstrat,maxiter = 200, maxiter1 = 25,
                                  realAllocation=realAllocation)[1]))
    return(res)
  }
  evaluateRcppMem <- memoise::memoise(evaluateRcpp)
  censiti<-0

  if (vars > 0) {
    if (!is.null(suggestions)) {
      if (verbose)
        cat("Adding suggestions to first population...\n")
      #population = matrix(nrow = popSize, ncol = vars)
      population<-matrix(nrow = popSize, ncol = nrow(stra))
      suggestionCount = ceiling(popSize*elitism)
      for (id in 1:suggestionCount) {
        population[id, ] = suggestions$suggestions
      }
      if (verbose)
        cat("Filling others with random values in the given domains...\n")
      for (var in 1:vars) {


        population[(suggestionCount + 1):popSize, var] = stringMin[var] +
          runif(popSize - suggestionCount) * (stringMax[var] -
                                                stringMin[var])
      }
    }    else {
      if (verbose)
      { cat("Generating stratifications...\n")}

      population<-matrix(nrow = popSize, ncol = nrow(stra))


      #population[1,]<-cluster::pam(stra[,3:(2+nvar)],initialStr)$cluster
      population[1,]<-  e1071::cmeans(stra[,3:(2+nvar)],
                    initialStr, iter.max = 1000,
                    method = "cmeans")$cluster
      pop_i<-2
      groups <- as.factor(population[1,])
       levels(groups) <- c(1:length(levels(groups)))
      for (pop_i in 2:popSize) { # don't mutate the best
        population[pop_i,]<-population[1,]
        population[pop_i,sample(vars,1)]<- as.numeric(sample(levels(groups),1))
        }

}

    #frame<-dataset
    # bestEvals = rep(NA, iters)
    bestEvals<-NULL
    meanEvals = NULL
    evalVals = rep(NA, popSize)
    Strata = rep(NA, popSize)
    popgrp<-list()

    AllIters<-0
    # elitism= floor(popSize/5)
    minnumstr=minnumstrat


    for (object in 1:popSize) {
      if (is.na(evalVals[object])) {
        #---- Modification:
        res <- evaluateRcppMem(population[object,])
        evalVals[object] = res
        #---- End modification:
        if (verbose) cat(".");
      }
    }
    # cat("Min evals ", min(evalVals),"\n")

    bestEvals<- c(bestEvals,min(evalVals))

    #plot(bestEvals,type="l")
    SolutionPopulation<-population

    reorderedPop<-SolutionPopulation[order(evalVals),]



    for (it in 1:iters) {
      Alliters<-it
      if (verbose==TRUE){
        cat(paste("Starting iteration", it, "\n"))}



      meanEvals[it] = mean(evalVals)

      if(SAArun==TRUE && (it %% SAAiters)==0){

        reorderedPop<- reorderedPop
        solution<- reorderedPop[1,]
        sugg1<-suggestions

          sugg1$suggestions<- solution
          res<-SAA(stra, err,
                   sugg1,
                   Temp,initialStrata=initialStr, decrement_constant, end_time,
                   showSettings, jsize,length_of_markov_chain,
                   verbose, dominio,minnumstrat,kmax_percent,ProbNewStratum,
                   strcens,writeFiles, showPlot=FALSE, minTemp, realAllocation)
          # cat("SAA sample size", res$best,"\n")
          bestEvals<- c(bestEvals,res$best)
          solution<-res$solution
          AllIters<-AllIters+res$solutions_generated
          reorderedPop[1,]<-solution


      }else if((it %% EDAfreq)==0){

        probsTable<-list()


        #evalVals = apply(population,1,newevaluateRcppMem)

        # population<-reorderedPop
        # SolutionPopulation<-population
        #
        # reorderedPop<-SolutionPopulation[order(evalVals),]

        elitePop<-as.matrix(reorderedPop[1:elitismR,])
        if (class(elitePop)[1]=="numeric"){ nColumns<-length(elitePop)
        for(ja in 1:nColumns){
          probsTable[[ja]]<-table(elitePop)/sum(table(elitePop))
        }
        }else{
        nColumns<-ncol(elitePop)
        for(jb in 1:nColumns){
          probsTable[[jb]]<-table(elitePop[,jb])/sum(table(elitePop[,jb]))
        }}



        for(ic in (elitismR+1):nrow(reorderedPop)){
          #  for(ic in 2:nrow(population)){
          for(jc in 1:ncol(reorderedPop)){
            reorderedPop[ic,jc]<-as.numeric(names(probsTable[[jc]]))[sample(length(probsTable[[jc]])
                                                                            ,1,prob=probsTable[[jc]])]
          }
        }



        #}
        if (mutationChance > 0) {
          #                    if (verbose) cat("  applying mutations... ");
          # cat("  applying mutations... ");
          mutationCount = 0;
          for (object in (elitismR+1):popSize) { # don't mutate the best
            for (var in 1:nrow(stra)) {
              if (runif(1) < mutationChance) { # ok, do mutation
                genoma <- as.factor(reorderedPop[object,])
                levels(genoma) <- c(1:length(levels(genoma)))
                reorderedPop[object,] <- genoma
                if (runif(1) <= (1-addStrataFactor)) {
                  mutation <- as.numeric(sample(levels(genoma),1))
                }
                else  {
                  mutation <- max(as.numeric(levels(genoma)))+1
                }

                # apply mutation, and delete known evalutation value
                #&n bsp;
                reorderedPop[object,var] = mutation;
                Strata[object] = NA;
                mutationCount = mutationCount + 1;
              }
            }
          }


        }
      }
      if (verbose==TRUE){
        cat("Calucating evaluation  values... ")}
        population<-reorderedPop
        for (object in 1:popSize) {
          #---- Modification:
          res <- evaluateRcppMem(population[object,])
          evalVals[object] = res
          #---- End modification:
          if (verbose) cat(".");
        }
        # cat("Min evals ", min(evalVals),"\n")

        bestEvals<- c(bestEvals,min(evalVals))
        meanEvals[it] = mean(evalVals)
        SolutionPopulation<-population

        reorderedPop<-SolutionPopulation[order(evalVals),]





     # population<-reorderedPop

      #}

        #plot(bestEvals,type="l")
    }

  }

  # cat("Min evals ", min(evalVals),"\n")

  result = list(stringMin = stringMin,
                stringMax = stringMax, popSize = popSize, iters = iters,
                suggestions = suggestions, population = population, elitism = elitismR,
                mutationChance = mutationChance, evaluations = evalVals, solution=solution,
                best = min(evalVals),#min(apply(population,1,evaluateRcppMem)),
                samplesizes=bestEvals, mean = meanEvals, TotalIterations=AllIters)

  return(result)

}

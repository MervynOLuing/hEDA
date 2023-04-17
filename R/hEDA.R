#'hybrid Estimation of Distribution Algorithm (hEDA)
#'@export
hEDA<-function (stra, err, suggestions =NULL,
                Temp=0.01,initialStr, decrement_constant=0.95, end_time =140,
                jsize=10,length_of_markov_chain =5,
                SAArun=TRUE,SAAiters=1000,
                popSize = 200, iters = 100, mutationChance = NA, elitism = NA,
                addStrataFactor=0.1, EDAfreq=EDAfreq,
                verbose = FALSE, dominio=dominio,minnumstrat=2,kmax_percent=0.025,ProbNewStratum=0.0001,
                strcens=FALSE,writeFiles=FALSE, showPlot=TRUE, minTemp = 0.000005, realAllocation=TRUE){

  stringMin <- rep(1, nrow(stra))
  stringMax <- rep(initialStr, nrow(stra))
  nvar<-length(grep("CV",names(err)))
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

  evaluateRcpp<-function(sugg){


    newstr<-aggrStrata(stra, nvar, floor(sugg), censiti,
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
      population = matrix(nrow = popSize, ncol = vars)
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
      { cat("Gerating stratifications using fuzzy clustering...\n")}

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
   # bestEvals = rep(NA, iters)
    bestEvals<-NULL
    meanEvals = NULL
    evalVals = rep(NA, popSize)
    Strata = rep(NA, popSize)
    popgrp<-list()

    AllIters<-0
    # elitism= floor(popSize/5)
    minnumstr=minnumstrat

    evalVals = apply(population,1,evaluateRcppMem)

    SolutionPopulation<-population

    reorderedPop<-SolutionPopulation[order(evalVals),]

    for (iter in 1:iters) {
      AllIters<-iter
      if (verbose==TRUE){
        cat(paste("Starting iteration", iter, "\n"))}


      if(SAArun==TRUE && (iter %% SAAiters)==0){
        reorderedPop<- reorderedPop
        #tot<-min(evalVals)


        #SolutionPopulation<-population

        #reorderedPop<-SolutionPopulation[order(evalVals),]
        # solution<- reorderedPop[which(evalVals==min(evalVals))[1],]
        solution<- reorderedPop[1,]
        sugg1<-suggestions
        ia<-1
        resSample<-NULL
        for(ia in 1:elitismR){
         # cat("ia ", ia, "\n")
          sugg1$suggestions<-reorderedPop[ia,]
          # outstrcor <- aggrStrata(stra, nvar,sugg1$suggestions, censiti,
          #                         dominio=dominio)
          # 
          # dimsamp <- nrow(outstrcor )
          # if (strcens == TRUE)
          #   outstrcor  <- rbind(outstrcor , cens)
          # dimens <- nrow(outstrcor )
          # outstrcor <-as.data.frame(outstrcor )
          # res1<-bethel_alfa(outstrcor , errors,realAllocation = realAllocation)
          # soluz1 <- res1[[1]]
          # cat("poptot test", sum(soluz1),"\n")

          res<-SAA(stra, err,
                   sugg1,
                   Temp,initialStrata=initialStr, decrement_constant, end_time,
                   showSettings, jsize,length_of_markov_chain,
                   verbose, dominio,minnumstrat,kmax_percent,ProbNewStratum,
                   strcens,writeFiles, showPlot=FALSE, minTemp, realAllocation)
           #cat("SAA sample size", res$best,"\n")
           resSample<-c(resSample,res$best)
           bestEvals<- c(bestEvals,min(resSample))
           # outstrcor <- aggrStrata(stra, nvar,res$solution, censiti,
           #                         dominio=dominio)
           #
           # dimsamp <- nrow(outstrcor )
           # if (strcens == TRUE)
           #   outstrcor  <- rbind(outstrcor , cens)
           # dimens <- nrow(outstrcor )
           # outstrcor <-as.data.frame(outstrcor )
           # res2<-bethel_alfa(outstrcor , errors,realAllocation = realAllocation)
           # soluz2 <- res2[[1]]
           # cat("Res test", sum(soluz2),"\n")
          solution<-res$solution
          AllIters<-AllIters+res$solutions_generated
          reorderedPop[ia,]<-solution
        }


      }else if((iter %% EDAfreq)==0){

        probsTable<-list()
        elitePop<-reorderedPop[1:elitismR,]
        if (class(elitePop)[1]=="numeric"){ nColumns<-length(elitePop)
        for(ja in 1:nColumns){
          probsTable[[ja]]<-table(elitePop)/sum(table(elitePop))
        }
        }else{nColumns<-ncol(elitePop)
        for(jb in 1:nColumns){
          probsTable[[jb]]<-table(elitePop[,jb])/sum(table(elitePop[,jb]))
        }
        }


        for(ic in (elitismR+1):nrow(reorderedPop)){
          #  for(ic in 2:nrow(population)){
          for(jc in 1:ncol(reorderedPop)){
            reorderedPop[ic,jc]<-as.numeric(names(probsTable[[jc]]))[sample(length(probsTable[[jc]])
                                                                            ,1,prob=probsTable[[jc]])]
          }
        }



        #}

      }

      if (mutationChance > 0) {
        #                    if (verbose) cat("  applying mutations... ");
        # cat("  applying mutations... ");
        mutationCount = 0;
        for (object in (elitismR+1):popSize) { # don't mutate the best
          for (var in 1:vars) {
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

      if (verbose==TRUE){
        cat("Calucating evaluation  values... ")}
      
      population<-reorderedPop
      evalVals = apply(population,1,evaluateRcppMem)
      #cat("Min evals ", min(evalVals),"\n")
      
      bestEvals<- c(bestEvals,min(evalVals))
      meanEvals[iter] = mean(evalVals)
      #plot(bestEvals,type="l")
      SolutionPopulation<-population
      
      reorderedPop<-SolutionPopulation[order(evalVals),]
      population<-reorderedPop
      #}


    }

  }
  evalVals = apply(population,1,evaluateRcppMem)
  result = list(stringMin = stringMin,
                stringMax = stringMax, popSize = popSize, iters = iters,
                suggestions = suggestions, population = population, elitism = elitismR,
                mutationChance = mutationChance, evaluations = evalVals,
                best = min(evalVals),#min(apply(population,1,evaluateRcppMem)),
                samplesizes=bestEvals, mean = meanEvals, TotalIterations=AllIters)

  return(result)

}

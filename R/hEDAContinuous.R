#'hybrid Estimation of Distribution Algorithm (hEDA) for continuous strata
#'@export
newhEDAContinuous<-function (frame, err, suggestions =NULL,
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
  ncuts = (initialStr - 1)
  stringMin = rep(0,ncuts*sum(grepl("X",colnames(frame))))
  stringMax = rep(1,ncuts*sum(grepl("X",colnames(frame))))
  stra <- buildStrataDF(frame,progress=FALSE, verbose=FALSE)
  nvar<-length(grep("CV",names(err)))
  vars = length(stringMin)
  nvars<-nrow(stra)
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
  nvar<-length(grep("CV",names(err)))
  model=NULL

  
  frame<-dataset
  
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
      { cat("Starting with random values in the given domains...\n")}
      pop = matrix(nrow = popSize, ncol = vars)
      # fill values
      for (var in 1:vars) {
        pop[,var] = stringMin[var] +
          runif(popSize)*(stringMax[var]-stringMin[var]);
      }
      i<-0
      population<-matrix(nrow = popSize, ncol = nrow(frame))
      
      dataset<-frame
      fr<-frame
      for (object in 1:popSize) {
        nX <- sum(grepl("X",colnames(fr)))
        for(i in 1:nX){
          eval(parse(text=paste("fr$ZZ",i,"<- fr$X",i,sep="")))
        }
        v<-pop[object,]
        for(j in 1:nX){
          ini=(j-1)*(NROW(v)/nX)+1
          fin=j*(NROW(v)/nX)
          eval(parse(text=paste("v",j,"<-string[ini:fin]*max(fr$ZZ",j,")",sep="")))
          eval(parse(text=paste("x",j,"_cuts<-as.data.frame(v",j,"[order(v",j,")],stringsAsFactors = TRUE)",sep="")))
          eval(parse(text=paste("x",j,"_cuts<-as.data.frame(rbind(min(fr$ZZ",j,")",",x",j,"_cuts,max(fr$ZZ",j,")),stringsAsFactors = TRUE)",sep="")))
          eval(parse(text=paste("x",j,"_cuts$lim<-x",j,"_cuts$`v",j,"[order(v",j,")]`",sep="")))
          eval(parse(text=paste("x",j,"_cuts$`v",j,"[order(v",j,")]`<-NULL",sep="")))
          eval(parse(text=paste("fr$X",j," <- NULL",sep="")))
        }
        
        for(i in 1:(ncuts+1)) {
          eval(parse(text=paste("fr$c",i,"<-0",sep="")))
          for(j in 1:nX) {
            eval(parse(text=paste("fr$c",i,"<-ifelse((fr$ZZ",j,">=x",j,"_cuts$lim[",i,"] & fr$ZZ",j,"<= x",j,"_cuts$lim[",i+1,"]),",i,",fr$c",i,")",sep="")))
          }  
        }
        fr$X1=apply(fr[,c((ncol(fr)-ncuts):ncol(fr))],1,max)
        fr$X1 <- as.factor(fr$X1)
        if (max(levels(fr$X1)) > length(levels(fr$X1))) {
          levels(fr$X1) <- c(1:length(levels(fr$X1)))
          fr$X1 <- droplevels(fr$X1)
          fr$X1 <- as.numeric(fr$X1)
        }
        fr$X1 <- as.numeric(fr$X1)
        population[object,]<-fr$X1
      }
      
    }
    frame<-dataset
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
      plot(bestEvals,type="l")
      SolutionPopulation<-population
      
      reorderedPop<-SolutionPopulation[order(evalVals),]
      
      meanEvals[it] = mean(evalVals)

      if(SAArun==TRUE && (it %% SAAiters)==0){

        reorderedPop<- reorderedPop
        #tot<-min(evalVals)


        #SolutionPopulation<-population

        #reorderedPop<-SolutionPopulation[order(evalVals),]
        # solution<- reorderedPop[which(evalVals==min(evalVals))[1],]
        solution<- reorderedPop[1,]
        sugg1<-suggestions
       

          # cat("ia ", ia, "\n")
          sugg1$suggestions<- solution
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
          # cat("SAA sample size", res$best,"\n")
          resSample<-c(resSample,res$best)
          bestEvals<- c(bestEvals,res$best)
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
          reorderedPop[1,]<-solution


      }else if((it %% EDAfreq)==0){

        probsTable<-list()


        #evalVals = apply(population,1,newevaluateRcppMem)

        # population<-reorderedPop
        # SolutionPopulation<-population
        # 
        # reorderedPop<-SolutionPopulation[order(evalVals),]

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
      }

  

      if (verbose==TRUE){
        cat("Calucating evaluation  values... ")}

     # population<-reorderedPop

      #}


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

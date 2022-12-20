#'Simulated Annealing Algorithm (SAA)
#'@export
SAA<-function (strata, errors, suggestions = NULL,
               Temp=NA,initialStrata, decrement_constant=0.95, end_time =140,
               showSettings = FALSE, jsize=100,length_of_markov_chain =50,
               verbose = FALSE, dominio=NULL,minnumstrat=NULL,kmax_percent=0.025,ProbNewStratum=NA,
               strcens=FALSE,writeFiles=FALSE, showPlot=TRUE, minTemp = 0.0005, realAllocation=TRUE)
{



  stringMin <- rep(1, nrow(strata))
  if(is.na(initialStrata)){
    stringMax <- (rep(ceiling(sqrt(nrow(strata)) *1), nrow(strata)))}else{  stringMax <- (rep(ceiling(initialStrata), nrow(strata)))}


  vars = length(stringMin)

  if (verbose)
  {cat("Testing the sanity of parameters...\n")}
  if (length(stringMin) != length(stringMax)) {
    stop("The vectors stringMin and stringMax must be of equal length.")
  }

  if (jsize < 1) {
    stop("The number of Markov Chains must be at least 1.")
  }

  start_time <- proc.time()

  if (vars > 0) {
    nvar <- length(grep("CV", names(errors)))
    strata_means<-as.matrix(strata[,3:(nvar+2)])

    ndom <- nrow(errors)
    vars = length(stringMin)


    if (!is.null(suggestions)){
      solution<-suggestions$suggestion
    }else {
      solution=sample.int(stringMax,vars,replace=TRUE,prob=NULL)
      #    solution<-cluster::pam(strata[,3:(nvar+2)],ceiling(sqrt(nrow(strata))),pamonce = 3)$cluster
    }

    kmax=ceiling(nrow(strata)*kmax_percent)




    if(is.na(length_of_markov_chain)){L_size<-2*(nrow(strata))}else{L_size<-length_of_markov_chain}

    soluz <- NULL
    v <- NULL
    dimens <- NULL
    censiti <- 0
    solution<-floor(solution)
    if( "matrix"%in%class(solution)){solution<-solution[1,]}
    strcor <- aggrStrata_RcppOpen(strata, nvar, solution, censiti,dominio=dominio)
    #strcor <- aggrStrata(strata, nvar, solution, censiti,
     #                    dominio=dominio)

    dimsamp <- nrow(strcor)
    if (strcens == TRUE)
      strcor <- rbind(strcor, cens)
    dimens <- nrow(strcor)
    # alfa<-c(rep(1/nvar, nvar))
    strcor<-as.data.frame(strcor)
    res<-bethel_alfa(strcor, errors,realAllocation = realAllocation)
    soluz <- res[[1]]
    alfa<- res[[2]]
    best_alfa<-alfa
    tot <- sum(soluz)
    # cat("Original sample size", tot, "\n")
    best_tot<-tot
    best_sol<-solution
    iters<-0
    bestIters<-NA
    x_axis <- NULL     # x axis
    y_axis <- NULL# y axis
    deltaList<-0
    j<-1
    k<-0
    newtotVector<-NULL
    L_chain<-tot
    Prev_tot<-tot
    newTemp<-Temp

    current_time<-proc.time() - start_time
    Nrow=(jsize*length_of_markov_chain)
    solutionPop<-matrix(NA, Nrow, ncol = nrow(strata))

    while(j < (jsize+1) & current_time[3] < end_time & newTemp > minTemp){


      Prev_tot<-min(L_chain)


      if (j ==1){(k<-kmax)}else{k=1}
      if (runif(1)<=1/(length_of_markov_chain)){
        groups<-unique(solution)
        newgroup<-max(groups)+1
        for (i in 1:length(solution)){
          if (runif(1)<=ProbNewStratum){
            solution[i]<-newgroup
          }}
      }
      alfa=c(rep(1/nvar, nvar))
      for(i in 1:L_size){


        iters<-iters+1
        #cat("iter   ",iters,"\n")
        newsolution<-solution

        groups<-unique(newsolution)
        grps<-groups[sample(length(groups),2,replace=FALSE)]
        orig_group<-grps[1]

        orig_group_records<- strata_means[which(newsolution %in% orig_group),]

        if (is.null(nrow(orig_group_records))){
          num_rec_orig<-length(orig_group_records)
        } else{
          num_rec_orig<-nrow(orig_group_records)
        }

        orig_loc <- sample(num_rec_orig,k,replace=TRUE)


        replace_group<-grps[2]
        newsolution[which(newsolution %in% orig_group)[orig_loc]]<-replace_group

        str<-NULL
        str <- cbind(strata, newsolution)
        strataReplace <- str[which(str$newsolution==replace_group),]
        strataOrig <-str[which(str$newsolution==orig_group),]
        strataDelta<-rbind(strataReplace,strataOrig)


       # strcorDelta <-aggrStrata(strataDelta, nvar=nvar,strataDelta$newsolution, censiti=censiti,
         #                        dominio=dominio)
        strcorDelta <-aggrStrata_RcppOpen(strataDelta, nvar=nvar,strataDelta$newsolution, censiti=censiti,
                                 dominio=dominio)
        strcorDelta <-as.data.frame(strcorDelta)
        newstrcor<-NULL
        newstrcor<-strcor[-which(strcor$STRATO %in% c(orig_group,replace_group)),]
        newstrcor<-rbind(newstrcor,strcorDelta)
        newstrcor<-newstrcor[order(newstrcor$STRATO),]
        dimsamp <- nrow(newstrcor)
        if (strcens == TRUE)
          strcor <- rbind(newstrcor, cens)
        dimens <- nrow(newstrcor)
        newstrcor<-as.data.frame(newstrcor)
        if(round(res[[2]],2)[1]>0.9){alfa=c(rep(1/nvar, nvar))}
        if(round(res[[2]],2)[2]>0.9){alfa=c(rep(1/nvar, nvar))}
        res<-second_bethel_alfa(newstrcor, errors,minnumstrat = 2, maxiter = 200,
                                maxiter1 = 25, alfa=alfa,realAllocation =realAllocation)
        soluz <- res[[1]]
        alfa<-res[[2]]

        newtot <- sum(soluz)



        delta<- newtot - tot
        newtotVector<-c(newtotVector,newtot)



        if (delta <=0) {
          solution<-newsolution
          tot<-newtot
          strcor<-newstrcor
          #solutionPop[iters,]<-solution
          # cat("Current best tot (updated) ", tot, "\n")
          best_sol<-newsolution
          best_tot<-tot
          best_alfa<-alfa
          #   outstrcor <- aggrStrata(strata, nvar, best_sol, censiti,
          #                           dominio=dominio)
          #
          #   dimsamp <- nrow(outstrcor )
          #   if (strcens == TRUE)
          #     outstrcor  <- rbind(outstrcor , cens)
          #   dimens <- nrow(outstrcor )
          #   outstrcor <-as.data.frame(outstrcor )
          #   res<-bethel_alfa(outstrcor , errors,minnumstrat = 2, maxiter = 20000,
          #                    maxiter1 = 25, realAllocation = realAllocation)
          #   soluz <- res[[1]]
          #   cat("### check current value ", sum(soluz), "\n")
          #   checktot<-sum(soluz)
          # if (!all.equal(checktot,tot)) { stop() }

        }else if ((exp((-delta/(newTemp))) > runif(1))){
          deltaList <- c(deltaList, delta)

          solution<-newsolution
          tot<-newtot
          strcor<-newstrcor
          #solutionPop[iters,]<-solution

        }else{
          solution<-solution
          tot<-tot
          strcor<-strcor
          # solutionPop[iters,]<-solution
          # cat("Else Tot ", tot, "\n")
        }


        y_axis <- c(y_axis, tot)

        x_axis <- c(x_axis, iters)

        L_chain<-c(L_chain,tot)






        current_time<-proc.time() - start_time

        i<-i+1

        if (j ==1){(k=ceiling(k*0.99))}else{k=1}
      }



      j<-j+1
      if(showPlot==TRUE){
        plot(y_axis,type="l",xlab="Number of Solutions",ylab="Sample Size")
      }

      newTemp<-decrement_constant*newTemp
    }

  }


  #best_sol<-solutionPop[bestIters,]

  # outstrcor <- aggrStrata(strata, nvar, best_sol, censiti,
  #                         dominio=dominio)
  outstrcor <-  aggrStrata_RcppOpen(strata, nvar, best_sol, censiti,
                      dominio=dominio)
  dimsamp <- nrow(outstrcor )
  if (strcens == TRUE)
    outstrcor  <- rbind(outstrcor , cens)
  dimens <- nrow(outstrcor )
  outstrcor <-as.data.frame(outstrcor )
  res<-bethel_alfa(outstrcor , errors,minnumstrat = 2, maxiter = 200,
                   maxiter1 = 25, realAllocation = realAllocation)
  soluz <- res[[1]]
  #cat("final sol value is ", sum(soluz), "\n")

  if (verbose)
  { cat(" done.\n")}

  result = list(maxiterations= jsize, j_reached=j-1, solutions_generated =iters,
                solution = best_sol,
                best = sum(soluz),samplesize=y_axis, deltaList=deltaList,
                Final_temperature=newTemp, best_alfa=best_alfa)


  return(result)
}

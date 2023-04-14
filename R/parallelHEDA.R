#'parallelhEDAContinuous
#'@export
parallelhEDAContinuous<-function(frame, cv,
                                 sugg,
                                 Temp=0.0001,initialStrata=15, decrement_constant=0.95, end_time =140,
                                 jsize=5,length_of_markov_chain =50,
                                 SAArun=TRUE,SAAiters=50,
                                 popSize = 20, iters = 100, mutationChance = 0.01, elitism = 0.1,
                                 addStrataFactor=0.000001, EDAfreq=1,
                                 verbose = FALSE, dominio=dominio,minnumstrat=2,kmax_percent=0.025,ProbNewStratum=0.0001,
                                 strcens=FALSE,writeFiles=FALSE, showPlot=TRUE, minTemp = 0.000005, realAllocation=TRUE){
 
if (writeFiles == TRUE) {
    dire <- getwd()
    direnew <- paste(dire, "/output", sep = "")
    if (dir.exists(direnew))
      unlink(direnew, recursive = TRUE)
    if (!dir.exists(direnew))
      dir.create(direnew)
  }
  require("foreach")
  require("parallel")
  require("doParallel")
  require("SamplingStrata")
  #thanks to: https://towardsdatascience.com/getting-started-with-parallel-programming-in-r-d5f801d43745
  dom<-unique(frame$domainvalue)
  ndom<-length(dom)
  cores<-(detectCores())-1
  if (ndom < cores) {cores <-ndom}
  cl <- parallel::makeCluster(cores)
  # Activate cluster for foreach library

  doParallel::registerDoParallel(cl)

  #ptm <- proc.time()
  r <- foreach::foreach(i = 1:ndom,
                        .combine = rbind,
                        #.packages = c("hEDA")
                        .packages = c("hEDA","Rcpp2doParallel","SamplingStrata")#,.verbose = TRUE
  ) %dopar% {


    if (!is.null(sugg)){
      suggestions =sugg[which(sugg$domainvalue==dom[i]),];
    }else {suggestions<-NULL}

    nvar=length(grep("CV",names(cv[i,])))
    fr=frame[which(frame$domainvalue==dom[i]),];err=cv[i,];
    Temp=Temp;nStrat=initialStrata[i]; decrement_constant=decrement_constant; end_time =end_time;
    jsize=jsize;length_of_markov_chain =length_of_markov_chain;
    SAArun=SAArun;SAAiters=SAAiters;
    popSize = popSize ; iters = iters; mutationChance = mutationChance; elitism =elitism ;
    addStrataFactor=addStrataFactor; EDAfreq=EDAfreq;
    verbose = verbose; dominio=dom[i];minnumstrat=minnumstrat;kmax_percent=kmax_percent;ProbNewStratum=ProbNewStratum;
    strcens=strcens;writeFiles=writeFiles; showPlot=showPlot; minTemp = minTemp; realAllocation=realAllocation
    censiti <-0
   hEDAContinuous(fr, err, suggestions ,
                   Temp,initialStr=nStrat, decrement_constant, end_time,
                   jsize,length_of_markov_chain,
                   SAArun,SAAiters,
                   popSize, iters, mutationChance, elitism,
                   addStrataFactor, EDAfreq,
                   verbose,dominio=dom[i],minnumstrat,kmax_percent,ProbNewStratum,
                   strcens,writeFiles, showPlot, minTemp, realAllocation)



  }
  #  time_foreach[3]
  # Stop cluster to free up resources
  parallel::stopCluster(cl)
  #  proc.time() - ptm
  #thanks to: https://towardsdatascience.com/getting-started-with-parallel-programming-in-r-d5f801d43745

  #  proc.time() - ptm

  if(showPlot==TRUE){


    for (i in 1:ndom){
      plot(unlist(r[i,]$samplesizes),type="l",xlab="Steps",ylab="Sample Size")

      title(paste("Domain #", dom[i], " - Sample cost",

                  round(min(unlist(r[i,]$best)), 2)),

            col.main = "red")
    }

  }

  if (writeFiles == TRUE) {

    if(ndom==1){

      stmt <- paste("png(filename = file.path(direnew, 'plotdom", 1, ".png'),height=5, width=7, units='in', res=144)", sep = "")

      eval(parse(text = stmt))



      plot(unlist(r[1,]$samplesizes),type="l",xlab="Number of Solutions",ylab="Sample Size")

      title(paste("Domain #", i, " - Sample cost",

                  round(min(unlist(r[1,]$best)), 2)),

            col.main = "red")

      if (writeFiles == TRUE)  dev.off()



    }else{

      for (i in 1:ndom){


        stmt <- paste("png(filename = file.path(direnew, 'plotdom", i, ".png'),height=5, width=7, units='in', res=144)", sep = "")

        eval(parse(text = stmt))

        plot(unlist(r[i,]$samplesizes),type="l",xlab="Steps",ylab="Sample Size")

        title(paste("Domain #", i, " - Sample cost",

                    round(min(unlist(r[i,]$best)), 2)),
              col.main = "red")

        if (writeFiles == TRUE)  dev.off()

      }
    }


  }

  # colnames(r)
  # [1] "stringMin"       "stringMax"       "popSize"         "iters"           "suggestions"     "population"
  # [7] "elitism"         "mutationChance"  "evaluations"     "best"            "samplesizes"     "mean"
  # [13] "TotalIterations"
  result = list( popSize=r[,3], population=r[,6],samplesizes=r[,11],steps=r[,4],suggestions=r[,5], SampleSize=r[,10])

  return(result)
}


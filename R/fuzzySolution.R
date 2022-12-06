#'fuzzySolution
#'@export
fuzzySolution<-function (strata, cv,m=2,
                         minClusters=2,
                      maxclusters=20, scale=TRUE,censiti=0)
{
require(e1071)
  nvar<-length(grep("CV",names(cv)))
DomSolutionList<-list()
ndom<-length(unique(strata$DOM1))
dom<-unique(strata$DOM1)
Domcost_vec<-rep(0,ndom)
for(DOM in 1:ndom){
  dom1<-strata[which(strata$DOM1==dom[DOM]),]
  data_train <- dom1[,3:(2+nvar)]
  if(scale==TRUE){
    #data_train_matrix <- scale(as.matrix(data_train))
    data_train_matrix <- scale(as.matrix(data_train))
  }else{data_train_matrix <- as.matrix(data_train)
  }

  solutionList<-list()
  costList<-list()

  cost_vec<-rep(0,(maxclusters))
  if (maxclusters >= (nrow(data_train_matrix)-2)){maxclusters<-(nrow(data_train_matrix)-2)}
  for(i in 1:maxclusters){
    #needs more than one cluster, hence k+1
    set.seed(1234)
    cluster.pm.fuzzycmeans <- e1071::cmeans(data_train_matrix, 1+i,iter.max = 100, m=m,  method="cmeans")$cluster
    cluster_assignment <- cluster.pm.fuzzycmeans
    strcor <- aggrStrata(strata[which(strata$DOM1==dom[DOM]),], nvar,
                         cluster_assignment, censiti,
                         dom[DOM])
    cost <- sum(SamplingStrata::bethel(strcor, cv[which(cv$domainvalue==dom[DOM]),],realAllocation = TRUE))
    cost_vec[i]<-cost
    solutionList[[i]]<-cluster_assignment
    costList[[i]]<-cost
    # cat("Clusters ",length(unique(cluster_assignment)), " cost ", cost, "\n")
  }

  DomSolutionList[[DOM]]<-unlist(solutionList[which(unlist(costList)==min(unlist(costList)))][1])
  Domcost_vec[DOM]<-min(unlist(costList))
}



OutputList<-list()
OutputList[[1]]<-DomSolutionList
OutputList[[2]]<-Domcost_vec
#sum(unlist(OutputList[[2]]))
#sum(unlist(Domcost_vec))
return(OutputList)
}

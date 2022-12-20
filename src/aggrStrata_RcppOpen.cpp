#include <RcppArmadillo.h>
using namespace Rcpp;
//using namespace sugar;
using namespace std;
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <algorithm>
//#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include <cmath>
#include <signal.h>
#include <unistd.h>
#include <pthread.h>
#include <array>
#include <float.h>
#include <sstream>
#include <RcppArmadilloExtensions/sample.h>

#include <map>
#include <string>
#include <string_view>

int n, i, j, k, p;
double num[100], sumOf=0.0, average;


// [[Rcpp::export]]
IntegerVector Rcpp_sortDelta(IntegerVector x) {
  IntegerVector y=x;
  // Order the elements of x by sorting y
  // First create a vector of indices
  IntegerVector idx = seq_along(x) - 1;
  // Then sort that vector by the values of y
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
  // And return x in that order
  return x[idx];
}

double sumof2=0.0;
// [[Rcpp::export]]
NumericMatrix   aggrStrata_RcppOpen(DataFrame strata, int nvar, IntegerVector vett,
                                    int censiti, int dominio) {

  strata["vett"] = vett;
  IntegerVector  varloop= seq(1, nvar);
  NumericVector N = strata["N"];
  NumericVector COST = strata["COST"];
  NumericVector CN = N *COST;
  double TMtot;

  IntegerVector STRATO = Rcpp_sortDelta(unique(vett));
  NumericMatrix x(vett.size(), nvar);
  NumericMatrix y(vett.size(),nvar);
  for (int i=0; i<nvar;i++) {
    y(_,i)=NumericVector(strata[(2+i)]);
  }

  int numCols = y.ncol();
  int numRows = y.nrow();
  int numGroups = STRATO.size();
  NumericVector sumVec;
  NumericVector AvVec(numGroups);
  NumericMatrix AvMat(nvar, numGroups);
  NumericMatrix TM(vett.size(),nvar);

  // double sumof2;
  NumericMatrix Centroids(numGroups,numCols);
  for(i = 0; i < (numCols); ++i)
  {

    for(k = 0; k < numGroups; ++k)
    {
      n=0;
      sumOf=0.0;
      sumof2=0.0;
      for(j = 0; j < numRows; ++j)
      {

        if (vett[j]==STRATO[k])
        {
          TMtot= y(j,i) *N[j];
          // TM(j,i)= TMtot;
          ;          sumOf += TMtot;
          sumof2+=N[j];
        }
      }
      average = sumOf/sumof2;
      AvMat(i, k)=average;
      Centroids(k,i)=average;
    }
  }


  for(i = 0; i < (numCols); ++i)
  {

    for(k = 0; k < numGroups; ++k)
    {
      for(j = 0; j < numRows; ++j){

        if (vett[j]==STRATO[k])
        {

          TM(j,i)= y(j,i) *N[j];
        }
      }
      average = sumOf/sumof2;
      AvMat(i, k)=average;
    }
  }



  // NumericMatrix Centroids(numGroups,numCols);
  //
  //
  //   for(i = 0; i < (numCols); ++i)
  //   {
  //     Centroids.row(i)=AvMat.column(i);
  //   }
  //
  NumericVector SumSquares(numGroups);
  NumericVector Nh(numGroups);

  NumericMatrix std_dev(vett.size(),nvar);
  for (int i=0; i<nvar;i++) {
    std_dev(_,i)=NumericVector(strata[((2+nvar)+i)]);
  }

  // NumericMatrix SquareMat =transpose(AvMat);


  double  TVARdouble;


  NumericMatrix TMt(nvar,numGroups);
  NumericMatrix TVAR(vett.size(),nvar);
  //double sumof2=0.0;
  double  TMtdouble;
  for(i = 0; i < (nvar); ++i)
  {

    for(k = 0; k < numGroups; ++k)
    {
      sumOf=0.0;
      TVARdouble=0.0;
      TMtdouble=0.0;
      //"TVAR1 <- strata$S1**2 * (strata$N - 1)"
      for(j = 0; j < numRows; ++j)
      {

        if (vett[j]==STRATO[k])
        {

          TVAR(j,i)= (pow(std_dev(j,i),2))*(N[j]-1);
          TVARdouble+= pow(std_dev(j,i),2);
          TMtdouble+= (y(j,i) *N[j]);
        }
      }

      TMt(i,k) = TMtdouble;
    }





  }

  for(k = 0; k < numGroups; ++k)
  {
    sumof2=0.0;
    for(j = 0; j < numRows; ++j)
    {

      if (vett[j]==STRATO[k])
      {
        sumof2+=N[j];
      }
      //  Nhg=sumof2;
    }
    Nh[k]=sumof2;
  }

  NumericMatrix diff(vett.size(),nvar);
  for(i = 0; i < (nvar); ++i)
  {

    for(k = 0; k < numGroups; ++k)
    {
      //"TVAR1 <- strata$S1**2 * (strata$N - 1)"
      for(j = 0; j < numRows; ++j)
      {

        if (vett[j]==STRATO[k])
        {
          //  "strwrk$diff2 <- strwrk$N * ((1/strwrk$N)*strwrk$TM2 - (1/strwrk$Nt)*strwrk$TM2t)**2"
          //2*((((1/2)*25)-(1/1028)*(53191)))^2
          diff(j,i)=N[j]*(pow((((1/N[j])*TM(j,i))-((1/Nh[k])*TMt(i,k))),2));



        }
      }
    }
  }



  //}


  NumericMatrix aggrDiff(nvar,numGroups);
  double diffT;
  for(i = 0; i < (nvar); ++i)
  {

    for(k = 0; k < numGroups; ++k)
    {

      diffT=0.0;
      //"TVAR1 <- strata$S1**2 * (strata$N - 1)"
      for(j = 0; j < numRows; ++j)
      {

        if (vett[j]==STRATO[k])
        {
          diffT+= diff(j,i);
        }

      }
      aggrDiff(i,k) =  diffT;

    }
  }

  NumericMatrix TVARt(nvar,numGroups);
  double TVARtdouble;
  for(i = 0; i < (nvar); ++i)
  {
    for(k = 0; k < numGroups; ++k)
    {
      TVARtdouble=0.0;
      //"TVAR1 <- strata$S1**2 * (strata$N - 1)"
      for(j = 0; j < numRows; ++j)
      {
        if (vett[j]==STRATO[k])
        {
          TVARtdouble+= TVAR(j,i);
        }
      }
      TVARt(i,k) =  TVARtdouble;
    }
  }

  NumericMatrix copyTVARt= TVARt;
  NumericMatrix S_Mat(nvar,numGroups);
  double Nhk;
  double TVARt_ik;
  double diff_ik;
  double Var_ik;
  double Tvar_plus_diff=0.0;
  Tvar_plus_diff=0.0+0.0;
  NumericMatrix Var_Mat(nvar,numGroups);
  for(i = 0; i < (nvar); ++i)
  {
    for(k = 0; k < numGroups; ++k)
    {
      //"S2 <- round(sqrt((1/strwrkagg$N)*(strwrkagg$TVAR2 + strwrkagg$diff2)),digits=4)"
      Nhk=Nh[k];
      TVARt_ik=copyTVARt(i,k);
      diff_ik=aggrDiff(i,k);
      Tvar_plus_diff=(TVARt_ik+diff_ik);
      ////Rcout << Tvar_plus_diff << "\n";
      Var_ik=((1/Nhk)*(TVARt_ik+diff_ik));

      //Var_Mat(i,k)=Var_ik;
      S_Mat(i,k)=sqrt(Var_ik);
    }
  }


  NumericVector CostAv(k);
  double CostSum;
  double count;

  for(k = 0; k < numGroups; ++k)
  {

    CostSum=0.0;
    count=0.0;
    //"TVAR1 <- strata$S1**2 * (strata$N - 1)"
    for(j = 0; j < numRows; ++j)
    {

      if (vett[j]==STRATO[k])
      {
        CostSum+= COST[j];
        count+=1;
      }

    }
    CostAv[k] =  CostSum/count;

  }
  IntegerVector DOM1(numGroups);
  for(k = 0; k < numGroups; ++k)
  {
    DOM1[k]=dominio;
  }

  IntegerVector CENS(numGroups);
  for(k = 0; k < numGroups; ++k)
  {
    CENS[k]=censiti;
  }
  // SquareMat(k,i) = sumOf/sumof2;

  //names(aggrStrata(strata,nvar,a,0,1))

  CharacterVector names((5+(2*nvar)));
  names[0]="STRATO";

  names[(1+(2*nvar))]="N";
  names[(2+(2*nvar))]="DOM1";
  names[(3+(2*nvar))]="COST";
  names[(4+(2*nvar))]="CENS";
  //
  Rcpp::List tmp(names.size());
  tmp[0] = STRATO;
  tmp[(1+(2*nvar))]=Nh;
  tmp[(2+(2*nvar))]=DOM1;
  tmp[(3+(2*nvar))]=CostAv;
  tmp[(4+(2*nvar))]=CENS;
  Rcpp::NumericMatrix strcor(numGroups,names.size());
  strcor.column(0) = STRATO;
  strcor.column(1+(2*nvar))=Nh;
  strcor.column(2+(2*nvar))=DOM1;
  strcor.column(3+(2*nvar))=CostAv;
  strcor.column(4+(2*nvar))=CENS;

  for(i = 0; i < (nvar); ++i)
  {

    names[1+i]="M" + std::to_string(i+1);
    names[((nvar+1)+i)]="S" + std::to_string(i+1);
    strcor.column(1+i)=round(Centroids.column(i),4);
    strcor.column((nvar+1)+i)=round(S_Mat.row(i),4);
  }



  //NumericMatrix T_AvMat=transpose(AvMat);

  // DataFrame strcor = DataFrame::create( Named("STRATO") = STRATO,
  //                    // AvMat,S_Mat,
  //                    Named("N") =Nh ,
  //                    Named("DOM1")=DOM1,
  //                    Named("COST")=CostAv,
  //                    Named("CENS") =CENS);

  //Rcpp::DataFrame result(tmp);
  //result.attr("names") = names;
  colnames(strcor) = names;
  return strcor;

}

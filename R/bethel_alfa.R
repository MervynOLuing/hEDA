# Adaptation of code from: G. Barcaroli, M. Ballin, H. Odendaal, D. Pagliuca, E. Willighagen, D. Zardetto (2020). SamplingStrata:
#   Optimal Stratification of Sampling Frames for Multipurpose Sampling Surveys. R package version 1.5-2 URL
# https://CRAN.R-project.org/package=SamplingStrata
#
# Giulio Barcaroli (2014). SamplingStrata: An R Package for the Optimization of Stratified Sampling.
# Journal of Statistical Software, 61(4), 1-24. URL http://www.jstatsoft.org/v61/i04/. doi
# 10.18637/jss.v061.i04
bethel_alfa <- function (stratif, errors, minnumstrat = 2, maxiter = 200, maxiter1 = 25,
          printa = FALSE, realAllocation = FALSE, epsilon = 0.00000000001)
{
  colnames(stratif) <- toupper(colnames(stratif))
  colnames(errors) <- toupper(colnames(errors))
  #checkData(strata = stratif, errors = errors)
  ordina_variabili <- function(dati, prefisso, n_var) {
    if (!is.data.frame(dati))
      stop()
    as.matrix(dati[, paste(prefisso, 1:n_var, sep = ""),
                   drop = FALSE])
  }
  iter1 <- 0
  val <- NULL
  m <- NULL
  s <- NULL
  cv <- NULL
  nstrat <- nrow(stratif)
  nvar <- length(grep("CV", names(errors)))
  ndom <- nrow(errors)
  varloop <- c(1:nvar)
  strloop <- c(1:nstrat)
  med <- ordina_variabili(stratif, "M", nvar)
  esse <- ordina_variabili(stratif, "S", nvar)
  if (ncol(med) != ncol(esse))
    stop(print("Error: Number of M variables don't match the number of S variables"))
  if (ncol(med) != nvar)
    stop(print("Error: Number of variables don't match the number of planned CV"))
  N <- as.vector(stratif$N)
  cens <- as.vector(stratif$CENS)
  cens[N < minnumstrat] <- 1
  cost <- as.vector(stratif$COST)
  if (is.null(cost))
    cost <- rep(1, nstrat)
  if (is.null(cens))
    cens <- rep(0, nstrat)
  nocens <- 1 - cens
  if (sum(cens) == length(cens)) {
    warning(print("Warning: Variable CENS always equal 1"))
  }
  nom_dom <- sapply(1:ndom, function(i) paste("DOM", i, sep = ""))
  dom <- ordina_variabili(stratif, "DOM", ndom)
  nvalues <- sapply(nom_dom, function(vari) {
    val <- c(val, nlevels(as.factor(dom[, vari])))
  })
  crea_disj <- function(data, vars) {
    out <- NULL
    sapply(vars, function(vari) {
      col <- as.factor(data[, vari])
      out <<- cbind(out, outer(col, levels(col), function(y,
                                                          x) ifelse(y == x, 1, 0)))
    })
    out
  }
  disj <- crea_disj(stratif, nom_dom)
  nc <- ncol(disj)
  for (i in 1:nc) {
    m <- cbind(m, disj[, i] * med)
    s <- cbind(s, disj[, i] * esse)
  }
  cvDom <- NULL
  cvDom2 <- NULL
  for (k in 1:ndom) {
    cvx <- ordina_variabili(errors[k, ], "CV", nvar)
    ndomvalues <- c(1:nvalues[k])
    for (k1 in ndomvalues) {
      cv <- cbind(cv, cvx)
      cvDom <- c(cvDom, rep(nom_dom[k], length(cvx)))
      cvDom2 <- c(cvDom2, rep(levels(as.factor(dom[, k]))[k1],
                              length(cvx)))
    }
  }
  nvar <- ncol(cv)
  varloop <- c(1:nvar)
  varfin <- c(rep(0, nvar))
  totm <- c(rep(0, nvar))
  alfa2 <- c(rep(0, nvar))
  crea_a <- function() {
    numA <- (N^2) * (s^2) * nocens
    denA1 <- colSums(t(t(N * m) * c(cv)))^2
    denA2 <- colSums(N * (s^2) * nocens)
    denA <- denA1 + denA2 + epsilon
    a <- t(t(numA)/denA)
    return(a)
  }
  chromy <- function(alfatot, diff, iter, alfa, alfanext,
                     x) {
    while (diff > epsilon && iter < maxiter) {
      iter <- iter + 1
      den1 <- sqrt(rowSums(t(t(a) * c(alfa))))
      den2 <- sum(sqrt(rowSums(t(t(a * cost) * c(alfa)))))
      x <- sqrt(cost)/(den1 * den2 + epsilon)
      alfatot <- sum(c(alfa) * (t(a) %*% x)^2)
      alfatot[alfatot == 0] <- epsilon
      alfanext <- c(alfa) * (t(a) %*% x)^2/alfatot
      diff <- max(abs(alfanext - alfa))
      alfa <- alfanext
      alfa2 <<- alfanext
    }
    if (realAllocation == FALSE)
      n <- ceiling(1/x)
    if (realAllocation == TRUE)
      n <- 1/x
    return(list(n,alfa))
  }
  check_n <- function() {
    n[n < minnumstrat] <- pmin(minnumstrat, N)[n < minnumstrat]
    n
  }
  a <- crea_a()
 res <- chromy(0, 999, 0, c(rep(1/nvar, nvar)), c(rep(0, nvar)),
              array(0.1, dim = c(nstrat, 1)))
 n<-unlist(res[1])
 alfa<-unlist(res[2])
  contx <- sum(n > N)
  cens[n > N] <- 1
  nocens <- 1 - cens
  n <- check_n()
  while (contx > 0 && iter1 < maxiter1) {
    iter1 <- iter1 + 1
    a <- crea_a()
    res <- chromy(0, 999, 0, c(rep(1/nvar, nvar)), c(rep(0, nvar)),
                  array(0.1, dim = c(nstrat, 1)))
    n<-unlist(res[1])
    alfa<-unlist(res[2])
    contx <- sum(n > N)
    cens[n > N] <- 1
    nocens <- 1 - cens
    n <- check_n()
  }
  n <- (nocens * n) + (cens * N)
  if (printa == TRUE) {
    stampa_confronto <- function(n, N, strato) {
      nomi <- c("STRATUM", "POPULATION", "BETHEL", "PROPORTIONAL",
                "EQUAL")
      df <- NULL
      df <- cbind(df, as.character(strato), N, n, ceiling(sum(n) *
                                                            N/sum(N)), ceiling(sum(n)/length(n)))
      tot <- apply(matrix(as.numeric(df[, 2:5]), ncol = 4),
                   2, sum)
      df <- rbind(df, c("TOTAL", tot))
      colnames(df) <- nomi
      df
    }
    calcola_cv <- function() {
      NTOT <- c(rep(0, nvar))
      CVfin <- c(rep(0, nvar))
      NTOT <- colSums((m > 0) * N)
      varfin <- rowSums(t((s * N)^2 * (1 - round(n)/N)/round(n))/NTOT^2)
      totm <- rowSums(t(m * N))
      CVfin <- round(sqrt(varfin/(totm/NTOT)^2), digits = 4)
      return(CVfin)
    }
    calcola_sensibilita <- function() {
      t <- g <- 0
      for (i in 1:nstrat) {
        t <- sum(alfa2 * a[i, ])
        g <- g + sqrt(cost[i] * t)
      }
      g <- g^2
      sens <- 2 * 0.1 * alfa2 * g
      return(sens)
    }
    stampa_cv <- function() {
      CVfin <- calcola_cv()
      sens <- calcola_sensibilita()
      domcard <- c(rep(0, ndom))
      dom <- as.data.frame(dom, stringsAsFactors = TRUE)
      for (i in 1:ndom) domcard <- c(domcard, nlevels(dom[,
                                                          i]))
      tit_cv <- c("TYPE", "DOMAIN/VAR.", "PLANNED CV ",
                  "ACTUAL CV", "SENSITIVITY 10%")
      outcv <- cbind(as.vector(cvDom), paste(as.vector(cvDom2),
                                             "/V", c(1:ncol(med)), sep = ""), as.vector(cv),
                     as.vector(CVfin), as.vector(ceiling(sens)))
      colnames(outcv) <- tit_cv
      return(outcv)
    }
    strato <- stratif$STRATO
    if (is.null(strato))
      strato <- paste("STR", 1:nstrat, sep = "")
    confr <- stampa_confronto(n, N, strato)
    outcv <- stampa_cv()
    attr(n, "confr") <- confr
    attr(n, "outcv") <- outcv
  }
  return(list(n,alfa))
}

# Thanks to: G. Barcaroli, M. Ballin, H. Odendaal, D. Pagliuca, E. Willighagen, D. Zardetto (2020). SamplingStrata:
#   Optimal Stratification of Sampling Frames for Multipurpose Sampling Surveys. R package version 1.5-2 URL
# https://CRAN.R-project.org/package=SamplingStrata
#
# Giulio Barcaroli (2014). SamplingStrata: An R Package for the Optimization of Stratified Sampling.
# Journal of Statistical Software, 61(4), 1-24. URL http://www.jstatsoft.org/v61/i04/. doi
# 10.18637/jss.v061.i04
#'aggrStrata
#'@export 
aggrStrata <- function(strata, nvar, vett, censiti, dominio) {
  strata <- cbind(strata, vett)
  varloop <- c(1:nvar)
  N <- strata$N
  CN <- strata$COST * strata$N 
  string2 <- ""
  string3 <- ""
  string5 <- ""
  string6 <- ""
  string8 <- ""
  string9 <- ""
  string10 <- ""
  string11 <- ""
  varloop <- c(1:nvar)
  for (k1 in varloop) {
    statement <- paste("TM", k1, " <- strata$M", k1, " * strata$N",
                       sep = "")
    eval(parse(text = statement))
    statement <- paste("TVAR", k1, " <- strata$S", k1, "**2 * (strata$N - 1)",
                       sep = "")
    eval(parse(text = statement))
    string2 <- paste(string2, "TM", k1, ",", sep = "")
    string3 <- paste(string3, "TVAR", k1, ",", sep = "")
    string5 <- paste(string5, "'TM", k1, "',", sep = "")
    string6 <- paste(string6, "'TM", k1, "t',", sep = "")
    string8 <- paste(string8, "'diff", k1, "',", sep = "")
    string9 <- paste(string9, "'TVAR", k1, "',", sep = "")
    string10 <- paste(string10, "M", k1, ",", sep = "")
    string11 <- paste(string11, "S", k1, ",", sep = "")
  }
  strwrk <- NULL
  strwrk2 <- NULL
  strwrkagg <- NULL
  strcor <- NULL
  statement <- paste("strwrk <- data.frame(gruppo=vett,", string2,
                     string3, "N,CN)", sep = "")
  eval(parse(text = statement))
  statement <- paste("strwrk2 <- aggregate(strwrk[,c(", string5,
                     "'N','CN')],by=list(vett),FUN=sum)", sep = "")
  eval(parse(text = statement))
  statement <- paste("colnames(strwrk2) <- c('gruppo',", string6,
                     "'Nt','COSTt')", sep = "")
  eval(parse(text = statement))
  strwrk <- merge(strwrk, strwrk2)
  rm(strwrk2)
  for (k1 in varloop) {
    statement <- paste("strwrk$diff", k1, " <- strwrk$N * ((1/strwrk$N)*strwrk$TM",
                       k1, " - (1/strwrk$Nt)*strwrk$TM", k1, "t)**2", sep = "")
    eval(parse(text = statement))
  }
  statement <- paste("strwrkagg <- aggregate(strwrk[,c(", string5,
                     string9, string8, "'N')],by=list(strwrk$gruppo),FUN=sum)",
                     sep = "")
  eval(parse(text = statement))
  for (k1 in varloop) {
    statement <- paste("M", k1, " <- round((strwrkagg$TM",
                       k1, " / strwrkagg$N),digits=4)", sep = "")
    eval(parse(text = statement))
    statement <- paste("S", k1, " <- round(sqrt((1/strwrkagg$N)*(strwrkagg$TVAR",
                       k1, " + strwrkagg$diff", k1, ")),digits=4)", sep = "")
    eval(parse(text = statement))
  }
  strato <- strwrkagg$Group.1
  N <- strwrkagg$N
  DOM1 <- c(rep(dominio, nrow(strwrkagg)))
  COST <- aggregate(strwrk$CN, by = list(strwrk$gruppo), FUN = sum)/aggregate(strwrk$N,
                                                                              by = list(strwrk$gruppo), FUN = sum)
  CENS <- c(rep(censiti, nrow(strwrkagg)))
  statement <- paste("strcor <- data.frame(strato,", string10,
                     string11, "N,DOM1,COST$x,CENS)", sep = "")
  eval(parse(text = statement))
  colnames(strcor) <- toupper(colnames(strcor))
  colnames(strcor)[colnames(strcor) == "COST.X"] <- c("COST")
  return(strcor)
}

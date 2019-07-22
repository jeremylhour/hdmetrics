#' ATE and ATT functions
#'
#' Computes ATE and ATT from data and nuisance parameters as in DoubleML p.35 (Chernozhukov et al., 2017)
#' Last edited: 11/01/2018
#' 
#' @param y vector of the dependent variable
#' @param D vector of treatment assignment
#' @param g1 estimate of E(Y_1 \vert X)
#' @param g0 estimate of E(Y_0 \vert X)
#' @param m propensity score estimate P(D=1 \vert X)
#' @param selec indicator of common support
#' 
#' All vectors above must have the same length.
#' 
#' @author Jeremy LHour

ATE <- function(Y,D,g1,g0,m,CS){
  if(missing(CS)){
    n = length(Y)
    CS = rep(TRUE,n)
  }
  gamma = (g1[CS]-g0[CS]) + D[CS]*(Y[CS]-g1[CS])/m[CS] - (1-D[CS])*(Y[CS]-g0[CS])/(1-m[CS])
  return(mean(gamma))
}

SE.ATE <- function(Y,D,g1,g0,m,CS){
  if(missing(CS)) CS = rep(TRUE,n)
  n = length(Y[CS])
  gamma = (g1[CS]-g0[CS]) + D[CS]*(Y[CS]-g1[CS])/m[CS] - (1-D[CS])*(Y[CS]-g0[CS])/(1-m[CS])
  return(sd(gamma)/sqrt(n))
}

ATT <- function(Y,D,g0,m,CS){
  if(missing(CS)) CS = rep(TRUE,length(Y))
  pi = mean(D[CS])
  gamma = D[CS]*(Y[CS]-g0[CS])/pi - (1-D[CS])*m[CS]*(Y[CS]-g0[CS])/(pi*(1-m[CS]))
  return(mean(gamma))
}

SE.ATT <- function(Y,D,g0,m,CS){
  if(missing(CS)) CS = rep(TRUE,length(Y))
  n = length(Y[CS])
  pi = mean(D[CS])
  gamma = D[CS]*(Y[CS]-g0[CS])/pi - (1-D[CS])*m[CS]*(Y[CS]-g0[CS])/(pi*(1-m[CS]))
  return(list(SE = sd(gamma)/sqrt(n),
              gamma = gamma))
}
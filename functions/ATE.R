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
#' 
#' All vectors above must have the same length.
#' 
#' @author Jeremy LHour

ATE <- function(Y,D,g1,g0,m){
  gamma = (g1-g0) + D*(Y-g1)/m - (1-D)*(Y-g0)/(1-m)
  return(mean(gamma))
}

SE.ATE <- function(Y,D,g1,g0,m){
  n = length(Y)
  gamma = (g1-g0) + D*(Y-g1)/m - (1-D)*(Y-g0)/(1-m)
  return(sd(gamma)/sqrt(n))
}

ATT <- function(Y,D,g0,m){
  pi = mean(D)
  gamma = D*(Y-g0)/pi - (1-D)*m*(Y-g0)/(pi*(1-m))
  return(mean(gamma))
}

SE.ATT <- function(Y,D,g0,m){
  n = length(Y)
  pi = mean(D)
  gamma = D*(Y-g0)/pi - (1-D)*m*(Y-g0)/(pi*(1-m))
  return(list(SE = sd(gamma)/sqrt(n),
              gamma = gamma))
}
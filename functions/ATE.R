#' ATE and ATT
#'
#' Computes ATE and ATT from data and nuisance parameters as in DoubleML p.35
#' Last edited: 05/01/2018
#' 
#' @param y vector of the dependent variable
#' @param D vector of treatment assignment
#' @param g1 estimate of E(Y_1 \vert X)
#' @param g0 estimate of E(Y_0 \vert X)
#' @param m propensity score estimate P(D=1 \vert X)
#' 
#' @author Jeremy LHour

ATE <- function(Y,D,g1,g0,m){
  psi = (g1-g0) + D*(Y-g1)/m - (1-D)*(Y-g0)/(1-m)
  return(mean(psi))
}

SE.ATE <- function(Y,D,g1,g0,m){
  n = length(Y)
  psi = (g1-g0) + D*(Y-g1)/m - (1-D)*(Y-g0)/(1-m)
  return(sd(psi)/sqrt(n))
}

ATT <- function(Y,D,g0,m){
  pi = mean(D)
  psi = D*(Y-g0)/pi - (1-D)*m*(Y-g0)/(pi*(1-m))
  return(mean(psi))
}

SE.ATT <- function(Y,D,g0,m){
  n = length(Y)
  pi = mean(D)
  psi = D*(Y-g0)/pi - (1-D)*m*(Y-g0)/(pi*(1-m))
  return(sd(psi)/(pi*sqrt(n)))
}
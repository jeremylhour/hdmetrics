#' DGP for Monte-Carlo Experiment
#' 
#' DataSim edition with Normal residuals and Progit propscore
#' Dummy variables are also included in covariates
#' Needs MASS library to run.
#' 
#' Last edited: 06/16/2017
#' 
#' @param n Sample size.
#' @param p Number of Covariates.
#' @param Ry Desired R-squared for outcome equation.
#' @param Rd Desired R-squared for treatment equation.
#' @param Intercept If TRUE returns an intercept also.
#' @param rho Correlation coefficients between two adjacent covariates.
#' 
#' @return Returns y,d and X necessary data to run Monte Carlo simulations. Also returns paramaters value for a possible check.
#' 
#' @author Jeremy Lhour

DataSim <- function(n=2000,p=50,Ry=.5,Rd=.2,Intercept=T,rho=.5,a=.5){
  
  ### Covariate variance matrix
  Sigma = matrix(0,nrow=p, ncol=p)
  
  for(k in 1:p){
    for(j in 1:p){
      Sigma[k,j] = rho^abs(k-j)
    }
  }
  
  ### Treatment variable coefficient
  b = rep(0,p)
  for(j in 1:abs(p/2)){
    b[j] = 1*(-1)^(j) / j^2
  }
  
  ### Outcome equation coefficients
  gamma = b
  for(j in (abs(p/2)+1):p){
    gamma[j] = (-1)^(j+1) / (p-j+1)^2
  }
  
  ### Adjustment to match R.squared
  c = sqrt((1/t(gamma)%*%Sigma%*%gamma)*(Rd/(1-Rd)))
  gamma = c*gamma
  
  c = sqrt((1/t(b)%*%Sigma%*%b)*(Ry/(1-Ry)))
  b = c*b
  
  ### Generate covariates
  X = mvrnorm(n = n, mu=rep(0,p), Sigma)
  
  # all even-indexed covariates are dummies
  even = 1:p %% 2 == 0
  X[,even] = ifelse(X[,even] > 0,1,0)
  
  # logit = function(x) 1/(1+exp(-x))
  d = as.numeric(runif(n) < pnorm(X%*%gamma))
  
  ### Treatment effect
  y = a*d + X%*%b + rnorm(n)
  
  if(Intercept) X = cbind(rep(1,n),X)
  
  return(list(X=X,
              y=y,
              d=d,
              b=b,
              g=gamma))
}
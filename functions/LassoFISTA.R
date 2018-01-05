#' Computes Lasso solution with FISTA
#'
#' Computes the Lasso solution using the FISTA method of Beck and Teboulle (2014).
#' The objective function is given as the mean squared plus lambda times the L1-norm.
#' Penalty-loadings for each coefficients are not allowed because the algorithm does not converge
#' in that case.
#' Last edited: 19 avril 2016.
#' 
#' @param betaInit Starting value for the coefficient. SHould be a vector of dim ncol(X). 
#' @param y vector of the dependent variable, normalizing it is a good idea.
#' @param X Matrix of Covariates, quicker if normalized.
#' @param W Vector of weights for each observation in the Least Square Objective.
#' @param nopen Set of indices of variables that should not be penalized.
#' @param lambda Overall penalty parameter.
#' @param tol Stopping criteion: difference between the value of objective function at two iterations.
#' @param maxIter Maximal number of iterations for the algorithm.
#' @param trace if TRUE print algorithm info.
#' 
#' @return beta argmin of the function.
#' @return value Value of the function at argmin.
#' @return loss Value of mean squared error.
#' @return l1norm l1-norm of beta.
#' @return nbIter Number of iterations necessary for convergence.
#' @return ConvergenceFISTA 0 if convergence, -555 if not.
#' 
#' @author Jeremy LHour


LassoFISTA <- function(betaInit=rep(0,ncol(X)),y,X,W=rep(1,nrow(X)),
                        nopen=NULL,lambda,
                        tol=1e-8,maxIter=1000,trace=F){
  ### Observation weighting
  W = as.vector(W)
  y = sqrt(W)*y
  X = sweep(X,MARGIN=1,sqrt(W),`*`)
  
  ### Set Algo. Values
  eta = 1/max(2*eigen(t(X)%*%X)$values/nrow(X))
  theta = 1
  thetaO = theta
  beta = betaInit
  v = beta
  cv = 0
  
  k = 0
  repeat{
    k = k+1
    
    thetaO = theta
    theta = (1+sqrt(1+4*thetaO^2))/2
    delta = (1-thetaO)/theta
    
    betaO = beta
    beta = prox(v - eta*LeastSqgrad(v,y,X), lambda*eta,nopen)
    
    v = (1-delta)*beta + delta*betaO
    
    # Show objective function value
    if(trace==T & k%%100 == 0){ print(paste("Objective Func. Value at iteration",k,":",LassoObj(beta,y,X,lambda,nopen))) }
    
    # Break if diverges
    if(is.na(LassoObj(beta,y,X,lambda,nopen) - LassoObj(betaO,y,X,lambda,nopen))){
      cv = -555
      print("LassoFISTA did not converge")
      break
    } else if(sum(abs(LassoObj(beta,y,X,lambda,nopen)-LassoObj(betaO,y,X,lambda,nopen))) < tol || k > maxIter) break
    
  }
  
  if(k > maxIter){
    print("Max. number of iterations reach in Lasso minimization.")
    cv = -555
  } 
  
  return(list(beta=as.vector(beta),
              value=LassoObj(beta,y,X,lambda,nopen),
              loss=LeastSq(beta,y,X),
              l1norm=sum(abs(beta)),
              nbIter=k,
              convergenceFISTA=cv))
}


#################################
#################################
### Define auxiliary functions###
#################################
#################################

prox <- function(x,lambda,nopen){
  y = (abs(x)-lambda)*(abs(x)-lambda > 0) * sign(x)
  y[nopen] = x[nopen] # Do not penalize these variables
  return(y)
}

LeastSq <- function(mu,y,X){
  X = as.matrix(X)
  return(mean((y - X%*%mu)^2))
}

LeastSqgrad <- function(mu,y,X){
  X = as.matrix(X)
  df = as.vector(-2*(t(y - X%*%mu)%*%X) / nrow(X))
  return(df)
}

LassoObj <- function(beta,y,X,lambda,nopen){
  if(length(nopen)>0){
    f = LeastSq(beta,y,X) + lambda*sum(abs(beta[-nopen]))
  } else {
    f = LeastSq(beta,y,X) + lambda*sum(abs(beta))
  }
  return(f)
}
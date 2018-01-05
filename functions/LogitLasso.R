#' Function to compute beta, LogitLasso estimate
#' 
#' Logit Lasso as in Van de Geer (2008). Uses the lbfgs function to perform L1-penalised minimization.
#'  A constant must be included as the first column in X.
#'  Last edited: 11 janvier 2016.
#' 
#' @param d Treatment indicator.
#' @param X Matrix of covariates.
#' @param c Constant for the overall penalty level.
#' @param maxIterPen Maximal number of iterations for penalty estimation.
#' @param trace if TRUE print convergence info.
#' @param PostLasso if TRUE computes the Post-Lasso solution.
#' @param maxIter Maximal number of iteration for Lasso program.
#' 
#' @return SHat Set of indices of non-zero elements in estimated beta.
#' @return betaLasso Lasso solution.
#' @return betaPL Post-Lasso solution.
#' @return lambda Overall penalty level.
#' @return psi Covariate-specific penalty loadings.
#' @return nbIter Number of iterations for penalty level estimation.
#' @return convergence 0 if convergence, -999 if not.
#' 
#' @author Jeremy Lhour


LogitLasso <- function(d,X,c=1.1,
                             maxIterPen=100,PostLasso=F,trace=F){
  ### Load necessary packages
  library("lbfgs")
 
  ### Setting
  d <- as.vector(d)
  d_tilde <- 2*d-1
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  ### First step: Lasso
  
  # Overall penalty level
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  
  # Penalty loadings: get a preliminary estimator
  beta <- c(log(sum(d)/(sum(1-d))),rep(0,p-1))
  Psi <- diag(as.vector(sqrt( t(1/(1+exp(d*(X%*%beta)))^2) %*% X^2 / n )))
  
  # Estimation parameters
  v <- .01 # Stopping rule
  k <- 0
  
  # Lasso estimate
  repeat{
    k <- k+1
    LassoEstim <- lbfgs(Logitloss, Logitlossgrad, Psi%*%beta, d=d_tilde, X=X%*%solve(Psi),
                        orthantwise_c=lambda,orthantwise_start=1,
                        invisible=1)
    beta <- solve(Psi) %*% LassoEstim$par
    
    # Update penalty loadings
    PrePsi <- Psi
    Psi <- diag(as.vector(sqrt( t(1/(1+exp(d*(X%*%beta)))^2) %*% X^2 / n )))
    
    if(trace & k%%5==0) print(paste("Max. pen. loading diff at Lasso Iteration nb.",k,":",max(abs(diag(Psi-PrePsi)))))
    if(k > maxIterPen || max(abs(diag(Psi-PrePsi))) < v) break
  }
  
  betaLasso <- beta 
  SHat <- union(1,which(betaLasso != 0)) #Always put a constant
  
  
  ### Second step: Post-Lasso
  betaPL <- rep(0,p)
  if(PostLasso==T){
    # PLtest <- glm(d ~ X[,SHat]-1, family = "binomial")
    PL <- lbfgs(Logitloss, Logitlossgrad, beta[SHat], d=d_tilde, X=X[,SHat],
                orthantwise_c=0, invisible=1)
    betaPL[SHat] <- PL$par
  }
  
  
  if(k > maxIterPen){
    print("Penalty estimation did not converge.")
    cvg=-999
  } else {
    cvg=0
  }
  
  return(list(SHat=SHat,
              betaPL=betaPL,
              betaLasso=c(betaLasso),
              lambda=lambda, 
              nbIter=k, 
              convergence=cvg,
              psi=diag(Psi)))
  
  
}

##################################
##################################
##################################
### Define Auxiliary Functions ###
##################################
##################################
##################################

Logitloss <- function(beta,d,X){
  X <- as.matrix(X)
  f <- log(1+exp(-d*(X%*%beta)))                                                                      
  return(mean(f,na.rm=T))
}

Logitlossgrad <- function(beta,d,X){
  X <- as.matrix(X)
  f <- d / (1+exp(d*(X%*%beta)))
  f <- -as.vector(f)*X
  return(as.vector(apply(f,2,mean,na.rm=T)))
}
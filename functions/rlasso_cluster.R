globalVariables(c("post", "intercept", "penalty", "control", "error", "n", "select.Z" , "select.X", "aes", "element_blank", "scale_x_discrete", "model.part", "all.categories", "X"))

#' rlasso: Function for Lasso estimation under homoscedastic and heteroscedastic non-Gaussian
#' disturbances
#'
#' The function estimates the coefficients of a Lasso regression with
#' data-driven penalty under homoscedasticity and heteroscedasticity with non-Gaussian noise and X-dependent or X-independent design. The
#' method of the data-driven penalty can be chosen. The object which is
#' returned is of the S3 class \code{rlasso}.
#'
#' The function estimates the coefficients of a Lasso regression with
#' data-driven penalty under homoscedasticity / heteroscedasticity and non-Gaussian noise. The options \code{homoscedastic} is a logical with \code{FALSE} by default.
#' Moreover, for the calculation of the penalty parameter it can be chosen, if the penalization parameter depends on the  design matrix (\code{X.dependent.lambda=TRUE}) or \code{independent} (default, \code{X.dependent.lambda=FALSE}).
#' The default value of the constant \code{c} is \code{1.1} in the post-Lasso case and \code{0.5} in the Lasso case. 
#'  A \emph{special} option is to set \code{homoscedastic} to \code{none} and to supply a values \code{lambda.start}. Then this value is used as penalty parameter with independent design and heteroscedastic errors to weight the regressors.
#' For details of the
#' implementation of the Algorithm for estimation of the data-driven penalty,
#' in particular the regressor-independent loadings, we refer to Appendix A in
#' Belloni et al. (2012). When the option "none" is chosen for \code{homoscedastic} (together with
#' \code{lambda.start}), lambda is set to \code{lambda.start} and the
#' regressor-independent loadings und heteroscedasticity are used. The options "X-dependent" and
#' "X-independent" under homoscedasticity are described in Belloni et al. (2013). 
# \code{lambda.start} can be component-specific. When used with one of the
# other option, the values are used as starting values.
#'
#' The option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of
#' the model with the selected variables.
#'
#' @aliases rlasso
#' @param post logical. If \code{TRUE}, post-Lasso estimation is conducted.
#' @param intercept logical. If \code{TRUE}, intercept is included which is not
#' penalized.
#' @param model logical. If \code{TRUE} (default), model matrix is returned.
#' @param penalty list with options for the calculation of the penalty. 
#' \itemize{
#' \item{\code{c} and \code{gamma}}{ constants for the penalty with default \code{c=1.1} and \code{gamma=0.1}}
#' \item{\code{homoscedastic}}{ logical, if homoscedastic errors are considered (default \code{FALSE}). Option \code{none} is described below.}
#' \item{\code{X.dependent.lambda}}{ logical,  \code{TRUE}, if the penalization parameter depends on the the design of the matrix \code{x}. \code{FALSE}, if independent of the design matrix  (default).}
#' \item{\code{numSim}}{ number of simulations for the dependent methods, default=5000}
#' \item{\code{lambda.start}}{ initial penalization value, compulsory for method "none"}
#' }
#' @param control list with control values.
#' \code{numIter} number of iterations for the algorithm for
#' the estimation of the variance and data-driven penalty, ie. loadings,
#' \code{tol} tolerance for improvement of the estimated variances.
#'\code{threshold} is applied to the final estimated lasso
#' coefficients. Absolute values below the threshold are set to zero.
#' @param ... further arguments (only for consistent defintion of methods)
#' @return \code{rlasso} returns an object of class \code{rlasso}. An object of
#' class "rlasso" is a list containing at least the following components:
#' \item{coefficients}{parameter estimates}
#' \item{beta}{parameter estimates (named vector of coefficients without intercept)}
#' \item{intercept}{value of the intercept}
#' \item{index}{index of selected
#' variables (logical vector)} \item{lambda}{data-driven penalty term for each
#' variable, product of lambda0 (the penalization parameter) and the loadings}
#' \item{lambda0}{penalty term} \item{loadings}{loading for each regressor}
#' \item{residuals}{residuals, response minus fitted values} \item{sigma}{root of the variance of
#' the residuals} \item{iter}{number of iterations} \item{call}{function call}
#' \item{options}{options}
#' \item{model}{model matrix (if \code{model = TRUE} in function call)}
#' @references A. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012).
#' Sparse models and methods for optimal instruments with an application to
#' eminent domain. \emph{Econometrica} 80 (6), 2369-2429.
#' @references A. Belloni, V. Chernozhukov and C. Hansen (2013). Inference for
#' high-dimensional sparse econometric models. In Advances in Economics and
#' Econometrics: 10th World Congress, Vol. 3: Econometrics, Cambirdge
#' University Press: Cambridge, 245-295.
#' @examples 
#' set.seed(1)
#' n = 100 #sample size
#' p = 100 # number of variables
#' s = 3 # nubmer of variables with non-zero coefficients
#' X = Xnames = matrix(rnorm(n*p), ncol=p)
#' colnames(Xnames) <- paste("V", 1:p, sep="")
#' beta = c(rep(5,s), rep(0,p-s))
#' Y = X%*%beta + rnorm(n)
#' reg.lasso <- rlasso(Y~Xnames)
#' Xnew = matrix(rnorm(n*p), ncol=p)  # new X
#' colnames(Xnew) <- paste("V", 1:p, sep="")
#' Ynew =  Xnew%*%beta + rnorm(n)  #new Y
#' yhat = predict(reg.lasso, newdata = Xnew)
#' @export
#' @rdname rlasso
rlasso_cluster <- function(x, ...) {
  UseMethod("rlasso_cluster") # definition generic function
  }
#' @param formula an object of class "formula" (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted in the form
#' \code{y~x}
#' @param data an optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model. If
#' not found in data, the variables are taken from environment(formula),
#' typically the environment from which \code{rlasso} is called.
#' @rdname rlasso
#' @export
rlasso_cluster.formula <- function(formula, data = NULL, post = TRUE, intercept = TRUE, model = TRUE, 
                           penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL, c = 1.1, gamma = .1/log(n)),
                          control = list(numIter = 15, tol = 10^-5, threshold = NULL), ...) {
  cl <- match.call()
  #if (missing(data))  data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 1
  y <- model.response(mf, "numeric")
  n <- length(y)
  x <- model.matrix(mt, mf)[,-1, drop=FALSE]
  if (missing(data)) {
    if (is.call(formula[[3]])) { 
    #colnames(x) <- sub(format(formula[[3]]), "", colnames(x))
    colnames(x) <- sub(re.escape(format(formula[[3]])), "", colnames(x))
    } else {
      colnames(x) <- sub(re.escape(formula[[3]]), "", colnames(x))  
    }
  }
  est <- rlasso_cluster(x, y, post = post, intercept = intercept, penalty=penalty, model=model, 
               control = control)
  est$call <- cl
  return(est)
}

#' @rdname rlasso
#' @export
rlasso_cluster.character <- function(x, data = NULL, post = TRUE, intercept = TRUE, model = TRUE, 
                           penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL, c = 1.1, gamma = .1/log(n)),
                           control = list(numIter = 15, tol = 10^-5, threshold = NULL), ...) {
  formula <- as.formula(x)
  res <- rlasso_cluster.formula(formula, data = data, post = post, intercept = intercept, model = model, 
                        penalty = penalty, control = control, ...)
}


#' @rdname rlasso
#' @export
#' @param y dependent variable (vector, matrix or object can be coerced to matrix)
#' @param x regressors (vector, matrix or object can be coerced to matrix)
rlasso_cluster.default <- function(x, y, T_t=1, post = TRUE, intercept = TRUE, model = TRUE,
                           penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL, c = 1.1, gamma = .1/log(n)),
                           control = list(numIter = 15, tol = 10^-5, threshold = NULL),...) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if (is.null(colnames(x)))
    colnames(x) <- paste("V", 1:p, sep = "")
  ind.names <- 1:p
  # set options to default values if missing
  if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  if (!exists("gamma", where = penalty))  penalty$gamma = 0.1/log(n)
  
  if (penalty$homoscedastic=="none" & !exists("lambda.start", where=penalty)) stop("lambda.start must be provided!")
  # checking input numIter, tol
  if (!exists("numIter", where = control)) {
    control$numIter = 15
  }
  
  if (!exists("tol", where = control)) {
    control$tol = 10^-5
  }
  
  #if (post==FALSE & (!exists("c", where = penalty) | is.na(match("penalty", names(as.list(match.call)))))) {
  if (post==FALSE & (!exists("c", where = penalty))) {  
    penalty$c = 0.5
  }
  
  # Intercept handling and scaling
  if (intercept) {
    meanx <- colMeans(x)
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- y - mu
  } else {
    meanx <- rep(0, p)
    mu <- 0
  }
  
  normx <- sqrt(apply(x, 2, var))
  ind <- rep(FALSE, p) #
  
  # variables with low variation are taken out, because normalization is not reliable
  # eps <- 10^-9  # precision for scaling
  #ind <- which(normx < eps)
  #if (length(ind) != 0) {
  #  x <- x[, -ind]
  #  normx <- normx[-ind]
  #  ind.names <- ind.names[-ind]
  #  p <- dim(x)[2]
  #  if (!is.null(penalty$lambda.start)) {
  #    penalty$lambda.start <- penalty$lambda.start[-ind]
  #  }
  #}
  
  #
  
  XX <- crossprod(x)
  Xy <- crossprod(x, y)
  
  startingval <- init_values(x,y)$residuals
  pen <- lambdaCalculation_cluster(penalty = penalty, y = startingval, x = x,T_t=T_t)
  lambda <- pen$lambda
  Ups0 <- Ups1 <- pen$Ups0
  lambda0 <- pen$lambda0
  
  mm <- 1
  s0 <- sqrt(var(y))
  
  while (mm <= control$numIter) {
    # calculation parameters
    #coefTemp <- LassoShooting.fit(x, y, lambda, XX = XX, Xy = Xy)$coefficients
    #xn <- t(t(x)/as.vector(Ups1))
    if (mm==1 && post) {
      coefTemp <- LassoShooting.fit(x, y, lambda/2, XX = XX, Xy = Xy)$coefficients
      #lasso.reg <- glmnet::glmnet(xn, y, family = c("gaussian"), alpha = 1,
      #                            lambda = lambda0/(2*n)/2, standardize = FALSE, intercept = FALSE)
      #lasso.reg <- glmnet::glmnet(x, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n)/2, standardize = FALSE, intercept = FALSE,  penalty.factor = Ups1)
    } else {
      coefTemp <- LassoShooting.fit(x, y, lambda, XX = XX, Xy = Xy)$coefficients
      #lasso.reg <- glmnet::glmnet(xn, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n), standardize = FALSE, intercept = FALSE)
      #lasso.reg <- glmnet::glmnet(x, y, family = c("gaussian"), alpha = 1,
      #                           lambda = lambda0/(2*n), standardize = FALSE, intercept = FALSE,  penalty.factor = Ups1)
    }
    #coefTemp <- as.vector(lasso.reg$beta)
    #names(coefTemp) <- colnames(x)
    coefTemp[is.na(coefTemp)] <- 0
    ind1 <- (abs(coefTemp) > 0)
    x1 <- as.matrix(x[, ind1, drop = FALSE])
    if (dim(x1)[2] == 0) {
      if (intercept) {
        intercept.value <- mean(y + mu)
        coef <- rep(0,p+1)
        names(coef) <- c("intercept", colnames(x)) #c("intercept", names(coefTemp))
      } else {
        intercept.value <- mean(y)
        coef <- rep(0,p)
        names(coef) <- colnames(x) #names(coefTemp)
      }
      est <- list(coefficients = coef, beta=rep(0,p), intercept=intercept.value, index = rep(FALSE, p),
                  lambda = lambda, lambda0 = lambda0, loadings = Ups0, residuals = y -
                    mean(y), sigma = var(y), iter = mm, call = match.call(),
                  options = list(post = post, intercept = intercept, ind.scale=ind, 
                                 control = control, mu = mu, meanx = meanx))
      if (model) {
        est$model <- x
      } else {
        est$model <- NULL
      }
      est$tss <- est$rss <- sum((y - mean(y))^2)
      est$dev <- y - mean(y)
      class(est) <- "rlasso"
      return(est)
    }
    
    # refinement variance estimation
    if (post) {
      reg <- lm(y ~ -1 + x1)
      coefT <- coef(reg)
      coefT[is.na(coefT)] <- 0
      e1 <- y - x1 %*% coefT
      coefTemp[ind1] <- coefT
    }
    if (!post) {
      e1 <- y - x1 %*% coefTemp[ind1]
    }
    s1 <- sqrt(var(e1))
    
    # homoscedatic and X-independent
    if (penalty$homoscedastic == TRUE && penalty$X.dependent.lambda == FALSE) {
      Ups1 <- s1*normx
      lambda <- rep(pen$lambda0 * s1, p)
    }
    # homoscedatic and X-dependent
    if (penalty$homoscedastic == TRUE && penalty$X.dependent.lambda == TRUE) {
      Ups1 <- s1*normx
      lambda <- rep(pen$lambda0 * s1, p)
    }
    # heteroscedastic and X-independent
    if (penalty$homoscedastic == FALSE && penalty$X.dependent.lambda == FALSE) {
#       Ups1 <- 1/sqrt(n) * sqrt(t(t(e1^2) %*% (x^2)))
      
      N= floor(n/T_t)
      e1_t = matrix(e1,N,T_t)
      #Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
      Ups1 <- NULL
      for( k in 1:p){
        x_t = matrix(x[,k],N,T_t)
        x_t = x_t*e1_t
        Ups1 <-c(Ups1, sum(rowSums(x_t)^2)  )
      }
      Ups1 <- 1/sqrt(n) *sqrt(Ups1)
      
      lambda <- pen$lambda0 * Ups1
    }
    
    # heteroscedastic and X-dependent
    if (penalty$homoscedastic == FALSE && penalty$X.dependent.lambda == TRUE) {
      lc <- lambdaCalculation_cluster(penalty, y=e1, x=x,T_t=T_t)
      Ups1 <- lc$Ups0
      lambda <- lc$lambda
    }
    
    
    
    # none
    if (penalty$homoscedastic == "none") {
      if (is.null(penalty$lambda.start)) stop("Argument lambda.start required!")
      Ups1 <- 1/sqrt(n) * sqrt(t(t(e1^2) %*% (x^2)))
      lambda <- pen$lambda0 * Ups1
    }
    
    mm <- mm + 1
    if (abs(s0 - s1) < control$tol) {
      break
    }
    s0 <- s1
  }
  
  if (dim(x1)[2] == 0) {
    coefTemp = NULL
    ind1 <- rep(0, p)
  }
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp) < control$threshold] <- 0
  ind1 <- as.vector(ind1)
  coefTemp <- as.vector(as.vector(coefTemp))
  names(coefTemp) <- names(ind1) <- colnames(x)
  if (intercept) {
    if (is.null(mu)) mu <-0
    if (is.null(meanx))  meanx <-  rep(0, length(coefTemp))  #<- 0
    if (sum(ind)==0) {
      intercept.value <- mu - sum(meanx*coefTemp)
    } else {
      intercept.value <- mu - sum(meanx*coefTemp) #sum(meanx[-ind]*coefTemp)
    }
  } else {
    intercept.value <- NA
  }
  
  #if (intercept) {
  #  e1 <- y - x1 %*% coefTemp[ind1] - intercept.value 
  #} else {
  #  e1 <- y - x1 %*% coefTemp[ind1]
  #}
  if (intercept) {
    beta <- c(intercept.value, coefTemp)
    names(beta)[1] <- "(Intercept)"
  } else {
    beta <- coefTemp
  }
  
  s1 <- sqrt(var(e1))
  est <- list(coefficients = beta, beta=coefTemp, intercept=intercept.value, index = ind1, lambda = lambda,
              lambda0 = lambda0, loadings = Ups1, residuals = as.vector(e1), sigma = s1,
              iter = mm, call = match.call(), options = list(post = post, intercept = intercept,
                                                             control = control, penalty = penalty, ind.scale=ind,
                                                             mu = mu, meanx = meanx), model=model)
  if (model) {
    x <- scale(x, -meanx, FALSE)
    est$model <- x
  } else {
    est$model <- NULL
  }
  est$tss <- sum((y - mean(y))^2)
  est$rss <- sum(est$residuals^2)
  est$dev <- y - mean(y)
  class(est) <- "rlasso"
  return(est)
}


################ function lambdaCalculation

#' Function for Calculation of the penalty parameter
#'
#' This function implements different methods for calculation of the penalization parameter \eqn{\lambda}. Further details can be found under \link{rlasso}.
#'
#' @param penalty list with options for the calculation of the penalty. 
#' \itemize{
#' \item{\code{c} and \code{gamma}}{ constants for the penalty with default \code{c=1.1} and \code{gamma=0.1}}
#' \item{\code{homoscedastic}}{ logical, if homoscedastic errors are considered (default \code{FALSE}). Option \code{none} is described below.}
#' \item{\code{X.dependent.lambda}}{ if \code{independent} or \code{dependent} design matrix \code{X} is assumed for calculation of the parameter \eqn{\lambda}}
#' \item{\code{numSim}}{ number of simulations for the X-dependent methods}
#' \item{\code{lambda.start}}{ initial penalization value, compulsory for method "none"}
#' }
#' @param x matrix of regressor variables
#' @param y residual which is used for calculation of the variance or the data-dependent loadings
#' @return The functions returns a list with the penalty \code{lambda} which is the product of \code{lambda0} and \code{Ups0}. \code{Ups0}
#' denotes either the variance (\code{independent} case) or the data-dependent loadings for the regressors. \code{method} gives the selected method for the calculation.
#' @export


lambdaCalculation_cluster <- function(penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL, c = 1.1, gamma = 0.1),
                              y = NULL, x = NULL, T_t=1) {
  checkmate::checkChoice(penalty$X.dependent.lambda, c(TRUE, FALSE, NULL))
  checkmate::checkChoice(penalty$homoscedastic, c(TRUE, FALSE, "none"))
  if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  if (!exists("c", where = penalty) & penalty$homoscedastic!="none") {
    penalty$c = 1.1
  }
  if (!exists("gamma", where = penalty) & penalty$homoscedastic!="none") {
    penalty$gamma = 0.1
  }


  # homoscedastic and X-independent
  if (penalty$homoscedastic==TRUE && penalty$X.dependent.lambda == FALSE) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 *
                                                                    p))
    Ups0 <- sqrt(var(y))
    lambda <- rep(lambda0 * Ups0, p)
  }

  # homoscedastic and X-dependent
  if (penalty$homoscedastic==TRUE && penalty$X.dependent.lambda == TRUE) {
    if (!exists("numSim", where = penalty)) {
      penalty$numSim = 5000
    }
    p <- dim(x)[2]
    n <- dim(x)[1]
    R <- penalty$numSim
    sim <- vector("numeric", length = R)
    for (l in 1:R) {
      g <- matrix(rep(rnorm(n), each = p), ncol = p, byrow = TRUE)
      sim[l] <- n * max(2 * colMeans(x * g))
    }
    lambda0 <- penalty$c * quantile(sim, probs = 1 - penalty$gamma)
    Ups0 <- sqrt(var(y))
    lambda <- rep(lambda0 * Ups0, p)
  }

  # heteroscedastic and X-independent (was "standard")
  if (penalty$homoscedastic==FALSE && penalty$X.dependent.lambda == FALSE) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    #lambda0 <- 2*penalty$c*sqrt(n)*sqrt(2*log(2*p*log(n)/penalty$gamma))
    lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 * p * 1))  # 1=num endogenous variables
    
    N= floor(n/T_t)
    y_t = matrix(y,N,T_t)
    #Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    Ups0 <- NULL
    for( k in 1:p){
      x_t = matrix(x[,k],N,T_t)
      x_t = x_t*y_t
      Ups0 <-c(Ups0, sum(rowSums(x_t)^2)  )
    }
    Ups0 <- 1/sqrt(n) *sqrt(Ups0)
    
    lambda <- lambda0 * Ups0
  }
  
  # heteroscedastic and X-dependent
  if (penalty$homoscedastic==FALSE && penalty$X.dependent.lambda == TRUE) {
    if (!exists("numSim", where = penalty)) {
      penalty$numSim = 5000
    }
    p <- dim(x)[2]
    n <- dim(x)[1]
    R <- penalty$numSim
    sim <- vector("numeric", length = R)
    lasso.x.y <- rlasso(y ~ x)
    eh <- lasso.x.y$residuals
    ehat <- matrix(rep(eh, each = p), ncol = p, byrow = TRUE) # might be improved by initial estimator or passed through
    for (l in 1:R) {
      g <- matrix(rep(rnorm(n), each = p), ncol = p, byrow = TRUE)
      sim[l] <- n * max(2 * colMeans(x * ehat* g)) # sqrt(n) or n??
    }
    lambda0 <- penalty$c * quantile(sim, probs = 1 - penalty$gamma)
    Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    lambda <- lambda0 * Ups0
  }
  
  
  if (!is.null(penalty$lambda.start)) {
    p <- dim(x)[2]
    if (length(penalty$lambda.start) == 1) {
      lambda.start <- rep(penalty$lambda.start, p)
    }
    lambda <- as.matrix(penalty$lambda.start)
  }

  if (penalty$homoscedastic == "none") {
    if (is.null(penalty$lambda.start) | !exists("lambda.start", where = penalty))
      stop("For method \"none\" lambda.start must be provided")
    n <- dim(x)[1]
    lambda0 <- penalty$lambda.start
    Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    lambda <- lambda0 * Ups0
  }

  return(list(lambda0 = lambda0, lambda = lambda, Ups0 = Ups0, method = penalty))
}





################# Methods for Lasso

#' Methods for S3 object \code{rlasso}
#'
#' Objects of class \code{rlasso} are constructed by \code{rlasso}.
#' \code{print.rlasso} prints and displays some information about fitted \code{rlasso} objects.
#' \code{summary.rlasso} summarizes information of a fitted \code{rlasso} object.
#' \code{predict.rlasso} predicts values based on a \code{rlasso} object.
#' \code{model.matrix.rlasso} constructs the model matrix of a \code{rlasso} object.
#'
#' @param object an object of class \code{rlasso}
#' @param x an object of class \code{rlasso}
#' @param all logical, indicates if coefficients of all variables (TRUE) should be displayed or only the non-zero ones (FALSE)
#' @param digits significant digits in printout
#' @param newdata new data set for prediction. An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are returned.
#' @param ... arguments passed to the print function and other methods
#' @rdname methods.rlasso
#' @aliases methods.rlasso print.rlasso predict.rlasso model.matrix.rlasso
#' @export

print.rlasso <- function(x, all=TRUE ,digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    if (all) {
      cat("Coefficients:\n")
      print.default(format(coef(x), digits = digits), print.gap = 2L,
                    quote = FALSE)
    } else {
      if (x$options$intercept) {
      print.default(format(coef(x)[c(TRUE,x$index)], digits = digits), print.gap = 2L,
                    quote = FALSE)
      } else {
        print.default(format(x$beta[x$index], digits = digits), print.gap = 2L,
                      quote = FALSE)
      }
    }
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#' @rdname methods.rlasso
#' @export

summary.rlasso <- function(object, all=TRUE, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nPost-Lasso Estimation: ",  paste(deparse(object$options$post), sep = "\n", collapse = "\n"), "\n", sep = " ")
  coefs <- object$coefficients
  p <- length(object$beta)
  num.selected <- sum(abs(object$beta)>0)
  n <- length(object$residuals)
  cat("\nTotal number of variables:", p)
  cat("\nNumber of selected variables:", num.selected, "\n", sep=" ")
  resid <- object$residuals
  cat("\nResiduals: \n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(apply(t(resid), 1L, quantile), dimnames = list(nam, dimnames(resid)[[2L]]))
  print(drop(t(rq)), digits = digits)
  cat("\n")
  if (all) {
    coefm <- matrix(NA, length(coefs), 1)
    coefm[,1] <- coefs
    colnames(coefm) <- "Estimate"
    rownames(coefm) <- names(coefs)
    printCoefmat(coefm, digits = digits, na.print = "NA")
  } else {
    coefs <- coefs[abs(coefs)>0]
    coefm <- matrix(NA, length(coefs), 1)
    coefm[,1] <- coefs
    colnames(coefm) <- "Estimate"
    rownames(coefm) <- names(coefs)
    printCoefmat(coefm, digits = digits, na.print = "NA")
  }
  cat("\nResidual standard error:", format(signif(object$sigma, digits)))
  cat("\n")
  
  if (object$options$intercept) {
    df.int <- 1
  } else {
    df.int <- 0
  }
  
  object$r.squared <- 1 - object$rss/object$tss
  object$adj.r.squared <- 1 - (1-object$r.squared)*((n-df.int)/(n-num.selected-df.int))
  cat("Multiple R-squared: ", formatC(object$r.squared, digits = digits))
  cat("\nAdjusted R-squared: ", formatC(object$adj.r.squared, digits = digits))
  
  if (!is.null(object$model)) {
    object$supscore <- sqrt(n)*max(abs(colMeans(object$model*as.vector(object$dev))))
    R <- 500
    stat <- vector("numeric", length=R)
    for (i in 1:R) {
      g <- rnorm(n)
      dev.g <- as.vector(g*object$dev)
      mat <- object$model*dev.g
      stat[i] <- sqrt(n)*max(abs(colMeans(mat)))
    }
    #quantstat <- quantile(stat, probs=c(1-alpha))
    object$pvalue <- sum(stat>object$supscore)/R
    cat("\nJoint significance test:\n", "the sup score statistic for joint significance test is", formatC(object$supscore, digits = digits), "with a p-value of", 
        formatC(object$pvalue, digits = digits))
  }
  cat("\n")
  invisible(object)
}

#' @rdname methods.rlasso
#' @export

# model.matrix.rlasso <- function(object, ...) {
#   if (is.call(object$call[[2]])) {
#     if(is.null(object$call$data)){
#       X <- model.frame(object$call[[2]])
#       mt <- attr(X, "terms")
#       attr(mt, "intercept") <- 0
#       mm <- model.matrix(mt, X)
#     } else {
#       dataev <- eval(object$call$data)
#       mm <- as.matrix(dataev[,names(object$beta)])
#     }
#   } else {
#     mm <- eval(object$call[[2]])
#   }
#   return(mm)
# }
model.matrix.rlasso <- function (object, ...) 
{
  if(is.null(object$model)){
    if (!is.null(object$call$x)){
      mm <- as.matrix(eval(object$call$x))
    } else {
      if(!is.null(object$call$data)){
        mm <- model.matrix(eval(object$call$formula), data = eval(object$call$data))[, -1, drop = FALSE]
      } else {
        mm <- model.matrix(eval(object$call$formula))[, -1, drop = FALSE]
      }
    }
  } else {
    mm <- object$model
  }
  return(mm)
}



#' @rdname methods.rlasso
#' @export

# predict.rlasso <- function (object, newdata = NULL, ...){
#   if (missing(newdata) || is.null(newdata)) {
#     if (is.matrix(model.matrix(object))) {
#       X <- model.matrix(object)
#     }
#     if (class(model.matrix(object))=="formula") {
#       form <- eval(model.matrix(object))
#       X <- model.matrix(form)[,-1, drop=FALSE]
#     }
#     if(sum(object$options$ind.scale)!=0) {
#       X <- X[,-object$options$ind.scale]
#     }
#   }
#   else {
#     varcoef <- names(object$beta)
#     if(all(is.element(varcoef,colnames(newdata)))){
#       X <- as.matrix(newdata[,varcoef])
#     } else {
#       if (is.matrix(newdata)) {
#         X <- as.matrix(newdata)
#       } else {
#         #X <- as.matrix(newdata)
#         formula <- eval(object$call[[2]])
#         X <- model.matrix(formula, data=newdata)[,-1, drop=FALSE]
#       }
#       if(sum(object$options$ind.scale)!=0) {
#         X <- X[,-object$options$ind.scale]
#       }
#       stopifnot(ncol(X)==length(object$beta))
#     }
#   }
#   n <- dim(X)[1] #length(object$residuals)
#   beta <- object$beta
#   
#   if (object$options[["intercept"]]) {
#     yhat <- X %*% beta + object$intercept
#     if (dim(X)[2]==0) yhat <- rep(object$intercept, n)
#   }
#   if (!object$options[["intercept"]]) {
#     yhat <- X %*% beta
#     if (dim(X)[2]==0) yhat <- rep(0, n)
#   }
#   return(yhat)
# }

# predict.rlasso2 <- function (object, newdata = NULL, ...) 
# {
#   if (missing(newdata) || is.null(newdata)) {
#     X <- model.matrix(object)
#     if (sum(object$options$ind.scale) != 0) {
#       X <- X[, -object$options$ind.scale]
#     }
#   } else {
#     varcoef <- names(object$beta)
#     if (all(is.element(varcoef, colnames(newdata)))){
#       X <- as.matrix(newdata[, varcoef])
#     } else {
#       if(!is.null(object$call$x)){
#         stop("newdata does not contain the correct variables")
#       } else {
#         mod.frame <- as.data.frame(cbind(rep(1, nrow(newdata)), newdata))
#         colnames(mod.frame)[1] <- as.character(eval(object$call$formula)[[2]])
#         X <- try(model.matrix(eval(object$call$formula), data = mod.frame)[, -1, drop = FALSE])
#         if(inherits(X, "try-error")){
#           stop("newdata does not contain the variables specified in formula")
#         }
#       }
#     } 
#     if (sum(object$options$ind.scale) != 0) {
#       X <- X[, -object$options$ind.scale]
#     }
#   }
#   n <- dim(X)[1]
#   beta <- object$beta
#   if (object$options[["intercept"]]) {
#     yhat <- X %*% beta + object$intercept
#     if (dim(X)[2] == 0) 
#       yhat <- rep(object$intercept, n)
#   }
#   if (!object$options[["intercept"]]) {
#     yhat <- X %*% beta
#     if (dim(X)[2] == 0) 
#       yhat <- rep(0, n)
#   }
#   return(yhat)
# }


predict.rlasso <- function (object, newdata = NULL, ...) 
{
  mf <- match.call(expand.dots = TRUE)
  m <- match("newx", names(mf), 0L)
  if (m!=0L) stop("Please use argument \"newdata\" instead of \"newx\" to provide data for prediction.")
  k <- length(object$beta)
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
    if (sum(object$options$ind.scale) != 0) {
      X <- X[, -object$options$ind.scale]
    }
  } else {
    if (is.null(colnames(newdata))) {
      X <- as.matrix(newdata)
      if (dim(X)[2] != k) {
        stop("No variable names provided in newdata and number of parameters does not fit!")
      } else {
        #message("No variable names provided in newdata. Prediction relies on right ordering of the variables.")
      }
    } else {
      varcoef <- names(object$beta)
      if (all(is.element(varcoef, colnames(newdata)))){
        X <- as.matrix(newdata[, varcoef])
      } else {
        mod.frame <- as.data.frame(cbind(rep(1, nrow(newdata)), newdata))
        colnames(mod.frame)[1] <- as.character(eval(object$call$formula)[[2]])
        X <- try(model.matrix(eval(object$call$formula), data = mod.frame)[, -1, drop = FALSE])
        if(inherits(X, "try-error")){
          stop("newdata does not contain the variables specified in formula")
        } 
      }
    }
  }
    if (sum(object$options$ind.scale) != 0) {
      X <- X[, -object$options$ind.scale]
    }
  
  n <- dim(X)[1]
  beta <- object$beta
  if (object$options[["intercept"]]) {
    yhat <- X %*% beta + object$intercept
    if (dim(X)[2] == 0) 
      yhat <- rep(object$intercept, n)
  }
  if (!object$options[["intercept"]]) {
    yhat <- X %*% beta
    if (dim(X)[2] == 0) 
      yhat <- rep(0, n)
  }
  return(yhat)
}













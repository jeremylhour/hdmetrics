format.perc <- function(probs, digits) paste(format(100 * probs, trim = TRUE, 
                                                    scientific = FALSE, digits = digits), "%")

# function for calculation of the errors after choosing the five
# variables with the highest correlation

init_values <- function(X, y, number = 5, intercept = TRUE) {
  suppressWarnings(corr <- abs(cor(y, X)))
  kx <- dim(X)[2]
  index <- order(corr, decreasing = T)[1:min(number, kx)]
  coefficients <- rep(0, kx)
  if (intercept == TRUE) {
    reg <- lm(y ~ X[, index, drop = FALSE])
    coefficients[index] <- coef(reg)[-1]
  } else {
    reg <- lm(y ~ -1 + X[, index, drop = FALSE])
    coefficients[index] <- coef(reg)
  }
  coefficients[is.na( coefficients)] <- 0
  res <- list(residuals = reg$residuals, coefficients = coefficients)
  return(res)
}


# function for escaping [ and other special characters in formula calls

re.escape <- function(strings){
  vals <- c("\\\\", "\\[", "\\]", "\\(", "\\)", 
            "\\{", "\\}", "\\^", "\\$","\\*", 
            "\\+", "\\?", "\\.", "\\|")
  replace.vals <- paste0("\\\\", vals)
  for(i in seq_along(vals)){
    strings <- gsub(vals[i], replace.vals[i], strings)
  }
  strings
}


######################## functions for extracting x,y,d,z from formula element

##### Define the function f.formula
f.formula <- function(formula, data, all.categories = FALSE){
  cl <- match.call()
  if (missing(data))  data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  ### if 'formula' and 'data' arguments are in reverse order: "switch", i.e. restore intended order/re-name
  if( any(c("formula", "Formula", "character") %in% is(data))  & 
      any(c("data.frame", "matrix") %in% is(formula)) ){
    mf$formula <- data
    mf$data <- formula
    formula <- mf$formula
    data <- mf$data
  }
  ### 
  

  matrix.flag <- FALSE
  if(is.matrix(data)){
    warning("'data' is a matrix and has been converted to a data.frame.")
    data <- as.data.frame(data)
    matrix.flag <- TRUE
  }
  
  ### Add (potentially "updated") 'data' to model.frame
  mf$data <- data
  
  formula <- Formula::as.Formula(formula)
  if(!identical(x=length(formula), y=c(1L, 2L))){
    if(length(formula)[1] < 1L) stop("The LHS of formula must include a response variable.", call. = FALSE)
    if(length(formula)[1] > 1L) stop("Multiple responses not implemented. The LHS of formula must include exactly one response variable.", call. = FALSE)
    if(length(formula)[2] != 2L) stop("RHS of formula is of improper length. formula must be specified as y ~ x + d | x + z, where y is response, x are exogenous regressors (on both sides of '|'), d are endogenous regressors (treatments), and z are instruments. Use of '.' to the right of '|' is permitted.", call. = FALSE)
  }
  
  f1 <- formula(formula, rhs = 1)            # "y ~ x + d" - part
  f2 <- formula(formula, lhs = 0, rhs = 2)   # "  ~ x + z" or " ~ . - d + z"
  formula <- Formula::as.Formula(f1, update(formula(formula, lhs = 0, rhs = 1), f2))  # update (only relevant if a `.` to the right of |)
  
  
  # Warning if (some) variable names are given as 'data[, 1]' or 'data[, "y"] instead of simply 'y'.
  if(any(grepl("[", as.character(formula), fixed=TRUE))){
    if(matrix.flag == TRUE){
      stop("Use of '[...]' for variable names, as in 'data[, 1]' or data[, \"y\"], is not possible if 'data' is a matrix. \n '[...]' works for most cases if 'data' is a data.frame, but its use is still discouraged.")
    }else{
      warning("Use of '[...]' for variable names, as in 'data[, 1]' or data[, \"y\"], is discouraged. It will work for most cases, especially when used consistently. Checking if the LHS variable is also included on the RHS is not available and setting 'all.categories' to TRUE has no effect.", call. = FALSE)      
    }
  }
  # Warning if (some) variable names are given as 'data$x' instead of simply 'x'.
  if(any(grepl("$", as.character(formula), fixed=TRUE))){
    if(matrix.flag == TRUE){
      stop("Use of '$' for variable names, as in 'data$y ~ data$x1 + data$d1 | data$x1 + data$z1', is not possible if 'data' is a matrix. \n '$' works for most cases if 'data' is a data.frame, but its use is still discouraged.")      
    }else{
      warning("Use of '$' for variable names, as in 'data$y ~ data$x1 + data$d1 | data$x1 + data$z1', is discouraged. It will work for most cases, especially when used consistently. Checking if the LHS variable is also included on the RHS is not available and setting 'all.categories' to TRUE has no effect.", call. = FALSE)  
    }
    
  }  
  # Check: were there are any variables from LHS used on the RHS (to the left or to the right of | )?
  if(  any(    setdiff(all.vars(formula(formula, lhs=1, rhs=0)), c( deparse(substitute(data)) )) %in%  # Note: remove the name of the data set to avoid wrong errors (in case of variable names given as data$y)
               all.vars(formula(formula, lhs=0, rhs=1:2)) )  ){
    if(missing(data)){
      stop("Argument 'data' is missing.", call. = FALSE)
    }else{
      stop("LHS and RHS of formula must not include the same variable(s).", call. = FALSE) 
    }
  }
  
  
  # Add updated formula to model.frame
  mf$formula <- formula
  #####
  
  
  mf[[1L]] <- as.name("model.frame")
  # Anm: as in ivreg(), whereas in lasso.formula:  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  
  #####
  ## Extract the response variable
  # - only one response variable allowed; stopped above with error() if more than one response variable
  # - extract -after transformation (if any)- response as a single column matrix
  
  Y <- as.matrix(model.part(formula, data = mf, lhs = 1, rhs=0, terms=FALSE, drop = FALSE))  
  # Anm: in ivreg and in lasso.formula:  Y <- model.response(mf, "numeric")
  # -> model.part() works even for multiple responses where model.response()  will fail
  #####
  
  
  
  #####
  ## Extract X, D, and Z
  # Reminder: formula is akin to  y ~ x + d | x + z
  # where x are exogenous, d endogenous, and z instruments (used only for the first stage)
  # the second part (x + z) contains exogenous variables and instruments
  # -> it has already been updated (see above)
  exin <- formula(formula, lhs=0, rhs=2)
  
  # the first part  (x + d) contains exogenous and endogenous variables
  exen <- formula(formula, lhs=0, rhs=1)
  
  # names of the components in X, D, and Z
  x <- intersect(attr(terms(exen), "term.labels"), attr(terms(exin), "term.labels"))
  d <- setdiff  (attr(terms(exen), "term.labels"), attr(terms(exin), "term.labels"))   # {x, d} \ {x, z} = {d}
  z <- setdiff  (attr(terms(exin), "term.labels"), attr(terms(exen), "term.labels"))   # {x, z} \ {x, d} = {z}
  
  if (is.environment(data)) data <- mf
  
  if (length(d)==0) {
    D <- NULL
  } else {
  conlist.D <- sapply(data, is.factor)                                                  # factor variables in data
  form.D <- as.formula(paste(" ~ ", paste(d, collapse= " + ", sep=""), sep=""))  # formula for D
  conlist.D[which(!(names(conlist.D) %in% attr(terms(form.D), "term.labels")))] <- FALSE  # TRUE for all factor variables in D
  D <- model.matrix(form.D, data = mf, contrasts.arg = lapply(data[, conlist.D, drop = FALSE], contrasts, contrasts = !all.categories))
  D <- D[, -which(colnames(D) == "(Intercept)"), drop = FALSE]
  }
  
  if (length(z)==0) {
    Z <- NULL
  } else {
  conlist.Z <- sapply(data, is.factor)                                                  # factor variables in data
  form.Z <- as.formula(paste(" ~ ", paste(z, collapse= " + ", sep=""), sep=""))  # formula for Z
  conlist.Z[which(!(names(conlist.Z) %in% attr(terms(form.Z), "term.labels")))] <- FALSE  # TRUE for all factor variables in Z
  Z <- model.matrix(form.Z, data = mf, contrasts.arg = lapply(data[, conlist.Z, drop = FALSE], contrasts, contrasts = !all.categories))
  Z <- Z[, -which(colnames(Z) == "(Intercept)"), drop = FALSE]
  }
  
  if (length(x)==0) {
    X <- NULL
  } else {
  conlist.X <- sapply(data, is.factor)                                                  # factor variables in data
  form.X <- as.formula(paste(" ~ ", paste(x, collapse= " + ", sep=""), sep=""))  # formula for X
  conlist.X[which(!(names(conlist.X) %in% attr(terms(form.X), "term.labels")))] <- FALSE  # TRUE for all factor variables in X
  X <- model.matrix(form.X, data = mf, contrasts.arg = lapply(data[, conlist.X, drop = FALSE], contrasts, contrasts = !all.categories))
  X <- X[, -which(colnames(X) == "(Intercept)"), drop = FALSE]
  }
  
  attr(D, "assign") <- NULL ; attr(Z, "assign") <- NULL ; attr(X, "assign") <- NULL
  #####
  
  
  ######
  #Extract intercepts for the x + d part (intercept2) and for the x + z part (intercept1): TRUE or FALSE.
  if(attr(terms(exen), "intercept") == 1L){
    intercept2 <- TRUE
  }else{
    intercept2 <- FALSE
  }
  if(attr(terms(exin), "intercept") == 1L){
    intercept1 <- TRUE
  }else{
    intercept1 <- FALSE
  }
  
  Y <- as.vector(Y)

  
  return(list(Y = Y, X = X, D = D, Z = Z, intercept1 = intercept1, intercept2 = intercept2))
  ######
}


############################## help function for tracing the position of selected variables specified by a one-sided formula

include <- function(I, data) {
  
I.orig <- I # save original I

if(!is.null(I)){
  # error message if I is not a formula/Formula or numeric/integer or logical:
  if(!any(c("formula", "Formula", "numeric", "integer", "logical") %in% class(I))){
    stop("'I' must be either an object of class formula/Formula or a numeric/integer or logical vector.")
  }
  if(any(c("numeric", "integer","logical") %in% class(I))){
    
    # convert the x + d part of the rhs of formula to character
    f1.no.y <- as.character(formula(formula, rhs = 1, lhs = 0))[[2]]
    # f1.no.y is to be split by the '+'-operator
    # if '-'-operator is used (e.g., y ~ -1 + x1 + ...), the '+' needs to be added beforehand:
    f1.no.y <- gsub(pattern = "-", replacement ="\\+ -", f1.no.y,
                    ignore.case = FALSE, perl = FALSE,
                    fixed = FALSE, useBytes = FALSE)
    # convert into character and splitting by '+'-operator:
    f1.no.y.string <- as.vector(unlist(strsplit(f1.no.y, "\\+ ")))
    # eliminate empty entries (i.e.: "")
    f1.no.y.string <- f1.no.y.string[which(f1.no.y.string != "")]
    # the elements in f1.no.y.string can now be indexed by numeric/integer or logical 'I'
    
    
    
    if(class(I.orig) %in% c("numeric", "integer")){
      # check: for numeric/integer I, no 'non-existent' elements can be selected
      # i.e. if 'formula' has k elements between '~' and '|', only 1, 2, ..., k are allowed
      # but -5, 0, k+1, 1.5, etc. are not allowed        
      if(!all(I.orig %in% 1:length(f1.no.y.string))){
        index.in <- I.orig %in% 1:length(f1.no.y.string) # logical index vector (problematic entries = FALSE)
        stop(paste("'formula' was determined to have ", length(f1.no.y.string),
                   " elements between '~' and '|'. The numeric/integer 'I' that was provided tried to select a non-existent element (e.g., 0, a negative number, a non-integer number, or a number greater than ", length(f1.no.y.string), "). The problematic entries in I were: ",
                   paste(I.orig[!index.in], collapse = ",", sep=" "), ".", sep=""))
      }
    }
    # check: for logical I, it has to have the appropriate length:
    if(class(I.orig) == "logical" & length(I.orig) != length(f1.no.y.string)){
      stop(paste("If 'I' is logical, its length has to be the same as the number of 'elements' (",
                 length(f1.no.y.string),
                 ") in the 'x + d'-part of 'formula'.", sep=""))
    }
    
    # I is now a formula as implied by the user-specified I, i.e. the 'formula' ''indexed'' by the user-specified I
    I <- as.formula(paste("~ ", paste(f1.no.y.string[I], collapse = " + "), sep=""))
  }
  
  
  # error message if I is too long.
  # formula-class: length must be 2
  # Formula-class: length must be c(0, 1)
  if("formula" %in% class(I)){
    if(!("Formula" %in% class(I))){
      if(length(I) > 2L){ stop("'I' must be a one-sided formula, e.g. '~ x1 + x2'.") }
    }else{
      if(length(I)[1] > 0L){ stop("'I' must be a one-sided formula: '~ x1 + x2' is ok, 'y ~ x1 + x2' is not.") }
      if(length(I)[2] > 1L){ stop("'I' must be a one-sided formula, '~ x1 + x2' is ok, '~ x1 + x2 | z1 + z2' is not.") }
    }
  }
  
  
  # names of the components in X that are selected via the argument I:
  i <- attr(terms(as.formula(I)), "term.labels")
  
  # error: if formula is like y ~ 1 + x1 + ..., or y ~ 0 + x1 + ..., or y ~ -1 + x1 + ...,
  # and if I is ~ 1, or c(1), or c(TRUE, FALSE, FALSE, ...),
  # then I did not selected any terms/variables, but only the intercept. 
  # In this case, the just constructed 'i' is 'character(0)', and operating on such an 'i'
  # would yield an uninformative error.
  # Instead so an informative error is produced in the next lines of code.
  if(length(i) == 0L){
    stop("Argument 'I' did not select any (non-intercept) terms from 'formula'. Most likely reason: 'formula' explicitly addressed the intercept (i.e., -1, 0, 1) and this is also the only component from 'formula' that is selected via 'I'.")
  }
  
  
  # construction of X.I: the model matrix X with only those columns which are selected via I:
  # (construction of X.I is analog to that of X, D, and Z)
  conlist.X.I <- sapply(data, is.factor)                                                  # factor variables in data
  form.X.I <- as.formula(paste(" ~ ", paste(i, collapse= " + ", sep=""), sep=""))  # formula for X.I
  conlist.X.I[which(!(names(conlist.X.I) %in% attr(terms(form.X.I), "term.labels")))] <- FALSE  # TRUE for all factor variables in X
  X.I <- model.matrix(form.X.I, data = data, contrasts.arg = lapply(data[, conlist.X.I, drop = FALSE], contrasts, contrasts = !all.categories))
  X.I <- X.I[, -which(colnames(X.I) == "(Intercept)"), drop = FALSE]
  attr(X.I, "assign") <- NULL
  
  # check: are all columns (columns names in X.I) also in X? if not: stop with error.
  if(!all(colnames(X.I) %in% colnames(X))){
    stop(paste("'I' resulted in column names that are not in X. Please only use terms in I as they appear in the 'x'-part of formula and no 'd'- or 'z'-terms. Also note that any x-term must be present both to the left and to the right of '|' in 'formula'. \n Encountered column names that result from I but are not in X: ",
               paste(setdiff(colnames(X.I), colnames(X)), collapse=", "),
               ".\n", sep=""))
  }
  
  I.c <- which(colnames(X) %in% colnames(X.I))
  
  
  # check: for each column in X.I: check if the column in X with the same name has the same values as the column in X.I
  # (a security check)
  equal.values.invariables <- function(number, X.matrix = X, X.I.matrix = X.I){
    col.name <- colnames(X.matrix)[number]  # column name of the column in 'X.matrix' given by 'number'
    # as all.equal returns TRUE or a message, and not simply TRUE/FALSE, use: isTRUE(all.equal(...))
    return( isTRUE(all.equal(X.matrix[, col.name, drop = FALSE], X.I.matrix[, col.name, drop = FALSE], check.attributes = FALSE)))
  } # function for comparison
  list.equal.values <- sapply(I.c, FUN = equal.values.invariables, X.matrix = X, X.I.matrix = X.I) # FALSE for those columns with non-matching entries
  # if there are any columns with non-matching entries: stop with error message that includes the colnames of the problematic columns:
  if(any(list.equal.values == FALSE)){
    stop( paste("'I' resulted in at least one column that has the same name in X and X.I (the design matrix created from 'I') but has different values in X and X.I. One possible scenario: a factor variable, say 'x5' (levels: A, B, C, D), is included in 'formula' so that dummies x5A, x5B, x5C, x5D are created for design matrix X. If there is already a (different) variable named, say 'x5B', in 'data' and it is included in 'I', 'x5B' in X.I will be the column from 'data' and not the one created by model.matrix() for X.\n The problematic column names are: ", 
                paste(colnames(X)[I.c[!list.equal.values]], collapse = ",", sep=", "), ".\n", sep="") )
  }
  
  
  
}

# generate I.c = NULL in case no I is provided, so that I.c = NULL can be passed on to workhorse function
if(is.null(I)){ I.c <- NULL }

return(I.c)
}

#######################################################################################

############################## help function for tracing the position of selected variables specified by a one-sided formula

# include <- function(I, data) {
#   
#   I.orig <- I # save original I
#   
#   if(!is.null(I)){
#     # error message if I is not a formula/Formula or numeric/integer or logical:
#     if(!any(c("formula", "Formula", "numeric", "integer", "logical") %in% class(I))){
#       stop("'I' must be either an object of class formula/Formula or a numeric/integer or logical vector.")
#     }
#     if(any(c("numeric", "integer","logical") %in% class(I))){
#       
#       # convert the x + d part of the rhs of formula to character
#       f1.no.y <- as.character(formula(formula, rhs = 1, lhs = 0))[[2]]
#       # f1.no.y is to be split by the '+'-operator
#       # if '-'-operator is used (e.g., y ~ -1 + x1 + ...), the '+' needs to be added beforehand:
#       f1.no.y <- gsub(pattern = "-", replacement ="\\+ -", f1.no.y,
#                       ignore.case = FALSE, perl = FALSE,
#                       fixed = FALSE, useBytes = FALSE)
#       # convert into character and splitting by '+'-operator:
#       f1.no.y.string <- as.vector(unlist(strsplit(f1.no.y, "\\+ ")))
#       # eliminate empty entries (i.e.: "")
#       f1.no.y.string <- f1.no.y.string[which(f1.no.y.string != "")]
#       # the elements in f1.no.y.string can now be indexed by numeric/integer or logical 'I'
#       
#       
#       
#       if(class(I.orig) %in% c("numeric", "integer")){
#         # check: for numeric/integer I, no 'non-existent' elements can be selected
#         # i.e. if 'formula' has k elements between '~' and '|', only 1, 2, ..., k are allowed
#         # but -5, 0, k+1, 1.5, etc. are not allowed        
#         if(!all(I.orig %in% 1:length(f1.no.y.string))){
#           index.in <- I.orig %in% 1:length(f1.no.y.string) # logical index vector (problematic entries = FALSE)
#           stop(paste("'formula' was determined to have ", length(f1.no.y.string),
#                      " elements between '~' and '|'. The numeric/integer 'I' that was provided tried to select a non-existent element (e.g., 0, a negative number, a non-integer number, or a number greater than ", length(f1.no.y.string), "). The problematic entries in I were: ",
#                      paste(I.orig[!index.in], collapse = ",", sep=" "), ".", sep=""))
#         }
#       }
#       # check: for logical I, it has to have the appropriate length:
#       if(class(I.orig) == "logical" & length(I.orig) != length(f1.no.y.string)){
#         stop(paste("If 'I' is logical, its length has to be the same as the number of 'elements' (",
#                    length(f1.no.y.string),
#                    ") in the 'x + d'-part of 'formula'.", sep=""))
#       }
#       
#       # I is now a formula as implied by the user-specified I, i.e. the 'formula' ''indexed'' by the user-specified I
#       I <- as.formula(paste("~ ", paste(f1.no.y.string[I], collapse = " + "), sep=""))
#     }
#     
#     
#     # error message if I is too long.
#     # formula-class: length must be 2
#     # Formula-class: length must be c(0, 1)
#     if("formula" %in% class(I)){
#       if(!("Formula" %in% class(I))){
#         if(length(I) > 2L){ stop("'I' must be a one-sided formula, e.g. '~ x1 + x2'.") }
#       }else{
#         if(length(I)[1] > 0L){ stop("'I' must be a one-sided formula: '~ x1 + x2' is ok, 'y ~ x1 + x2' is not.") }
#         if(length(I)[2] > 1L){ stop("'I' must be a one-sided formula, '~ x1 + x2' is ok, '~ x1 + x2 | z1 + z2' is not.") }
#       }
#     }
#     
#     
#     # names of the components in X that are selected via the argument I:
#     i <- attr(terms(as.formula(I)), "term.labels")
#     
#     # error: if formula is like y ~ 1 + x1 + ..., or y ~ 0 + x1 + ..., or y ~ -1 + x1 + ...,
#     # and if I is ~ 1, or c(1), or c(TRUE, FALSE, FALSE, ...),
#     # then I did not selected any terms/variables, but only the intercept. 
#     # In this case, the just constructed 'i' is 'character(0)', and operating on such an 'i'
#     # would yield an uninformative error.
#     # Instead so an informative error is produced in the next lines of code.
#     if(length(i) == 0L){
#       stop("Argument 'I' did not select any (non-intercept) terms from 'formula'. Most likely reason: 'formula' explicitly addressed the intercept (i.e., -1, 0, 1) and this is also the only component from 'formula' that is selected via 'I'.")
#     }
#     
#     
#     # construction of X.I: the model matrix X with only those columns which are selected via I:
#     # (construction of X.I is analog to that of X, D, and Z)
#     conlist.X.I <- sapply(data, is.factor)                                                  # factor variables in data
#     form.X.I <- as.formula(paste(" ~ ", paste(i, collapse= " + ", sep=""), sep=""))  # formula for X.I
#     conlist.X.I[which(!(names(conlist.X.I) %in% attr(terms(form.X.I), "term.labels")))] <- FALSE  # TRUE for all factor variables in X
#     X.I <- model.matrix(form.X.I, data = data, contrasts.arg = lapply(data[, conlist.X.I, drop = FALSE], contrasts, contrasts = !all.categories))
#     X.I <- X.I[, -which(colnames(X.I) == "(Intercept)"), drop = FALSE]
#     attr(X.I, "assign") <- NULL
#     
#     
#     I.c <- which(colnames(X) %in% colnames)
#   
#   # generate I.c = NULL in case no I is provided, so that I.c = NULL can be passed on to workhorse function
#   if(is.null(I)){ I.c <- NULL }
#   
#   return(I.c)
#   }
# }  

###########
  check_variables <- function(formula, varnames) {
    
    if (is.null(formula)) return(NULL)
    
    if("formula" %in% class(formula)){
      if(!("Formula" %in% class(formula))){
        if(length(formula) > 2L){ stop("'I' must be a one-sided formula, e.g. '~ x1 + x2'.") }
      }else{
        if(length(formula)[1] > 0L){ stop("'I' must be a one-sided formula: '~ x1 + x2' is ok, 'y ~ x1 + x2' is not.") }
        if(length(formula)[2] > 1L){ stop("'I' must be a one-sided formula, '~ x1 + x2' is ok, '~ x1 + x2 | z1 + z2' is not.") }
      }
    }
    
    i <- attr(terms(as.formula(formula)), "term.labels")
    I.c <- match(i, varnames)
    return(I.c)
  }


################## function for checking of binary variable
  
check_binary <- function(x) {
  
  if (!setequal(unique(x), c(0,1))) {
        stop("Treatment variable and Instrumental Variable should be binary (0/1)!")
  }
  return(invisible(NULL))
}

################## Function for constructing Instrumenatl Variables following BLP (1995)
constructIV <- function(firmid, cdid, id, X) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  names <- colnames(X)
  if (is.null(names)) names <- paste("V", 1:p, sep="")
  sum.other <-  matrix(NA, nrow = n, ncol = p)
  colnames(sum.other) <- paste("sum.other.", names, sep="")
  sum.rival <-  matrix(NA, nrow = n, ncol = p)
  colnames(sum.rival) <- paste("sum.rival.", names, sep="")
  
  for (i in 1:n) {
    for (j in 1:p) {
      other_ind=(firmid==firmid[i] & cdid==cdid[i] & id!=id[i])
      rival_ind=(firmid!=firmid[i] & cdid==cdid[i])
      
      sum.other[i,j]=sum(X[other_ind==1,j])
      sum.rival[i,j]=sum(X[rival_ind==1,j]) 
    }
  }
  return(cbind(sum.other, sum.rival))
}


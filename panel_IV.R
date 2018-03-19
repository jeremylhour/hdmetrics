###########################
## performances of cluster lasso in IV DML


library("mnormt")
library(hdm)
set.seed(2024)
T_t = 10
T_off = 5

n_seq = c(50,100,150,200)
j=3
n = n_seq[j]
### GENERATE DATA
rho_eps = 0.8
rho_nu = 0.5
rho_u = 0.8
rho_z = 0.8
alpha = 0.5

## Two different numbers of instruments
p_z = n*(T_t-2)
# p_z = n*(T_t+2)

## Designs for the coefficients vectors \pi 
## Design 1 
s= floor(0.5*n^{1/3})
Pi = (-1)^(0:(p_z-1))*(  1/sqrt(s)*matrix(1,p_z,1)   )
Pi[s:p_z,1] <- Pi[s:p_z,1]*sqrt(s)*(s:p_z)^(-2)


### individual heterogeneity
Sigma <- diag(rep(1,n))
for (i in 1:n){
  for (j in 1:n){
    Sigma[i,j] <- (1/2)^{abs(i-j)}
  }
}
Sigma <-Sigma *4/T
e_i <- rmnorm(1, mean = rep(0, nrow( Sigma)), varcov =  Sigma)
f_i = e_i

######  instruments
Sigma <- diag(rep(1,p_z))
for (i in 1:n){
  for (j in 1:n){
    Sigma[i,j] <- (1/2)^{abs(i-j)}
  }
}
var <- rmnorm(n*T_t, mean = rep(0, nrow(Sigma )), varcov = Sigma )
z_all = NULL
z_all_oracle = NULL
for (index in 1:p_z){
  phi =  matrix(var[,index],n,T_t)
  z <- matrix(e_i/(1-rho_z)+ sqrt(1/(1-rho_z^2))*phi[,1],n,1)
  for(i in 2:T_t){
    z  <- cbind(z ,rho_z*z[,i-1]+phi[,i])
  } 
  z_all <- cbind(z_all,z)
}
dim(z_all)

source("C:/Users/gaillac/Documents/GitHub/hdmetrics/functions/rlasso_cluster.R")
source("C:/Users/gaillac/Documents/GitHub/hdmetrics/functions/help_functions.R")


out = NULL
out1 = NULL
out2 = NULL
out3 = NULL
Ups1_all=NULL
MC= 1000
for( m in 1:MC){
  ####### errors
  sigmas <- matrix(c(1,rho_nu,rho_nu,1), 2,2  )
  var <- rmnorm(n*(T_t+T_off ), mean = rep(0, nrow(sigmas)), varcov = sigmas)
  nu_1 <- matrix(var[,1],n,(T_t+T_off ))
  nu_2 <- matrix( var[,2],n,(T_t+T_off))
  
  eps <-matrix(nu_1[,1],n,1)
  us <-matrix(nu_2[,1],n,1)
  for(i in 2:(T_t+T_off)){
    eps <- cbind(eps,rho_eps*eps[,i-1]+nu_1[,i])
    us <- cbind(us,rho_u*us[,i-1]+nu_2[,i])
  } 
  i=4
  acf(eps[i,T_off:(T_t+T_off)])
  acf(eps[i,1:(T_t+T_off)])
  acf(us[i,T_off:(T_t+T_off)])
  acf(us[i,1:(T_t+T_off)])
  
  eps <- eps[,T_off:(T_t+T_off)]
  us <- us[,T_off:(T_t+T_off)]
  
  d = NULL
  y = NULL
  index=1
  for (index in 1:T_t){
    sel_z = seq(index,index+T_t*(p_z-1), by=T_t)
    d_curr =  z_all[,sel_z]%*%Pi  + f_i + us[,index]
    d = cbind(d,  d_curr    )
    y = cbind(y,  alpha*d_curr  + e_i + eps[,index]     )
  }
  
  iv_oracle = NULL
  for (index in 1:T_t){
    sel_z = seq(index,index+T_t*(p_z-1), by=T_t)
    d_curr =  z_all[,sel_z]%*%Pi  + f_i 
    iv_oracle = cbind(iv_oracle,  d_curr    )
  }
  ##### estimation 
  
  ## within transformation
  
  y_dot = y - apply(y,1,mean)
  d_dot = d - apply(d,1,mean)
  iv_oracle_dot = iv_oracle - apply(iv_oracle ,1,mean)
  z_all_dot = z_all
  for (index in 1:p_z){
    z_all_dot[,((index-1)*T_t+1):(index*T_t)] <- z_all[,((index-1)*T_t+1):(index*T_t)]  - apply( z_all_dot[,((index-1)*T_t+1):(index*T_t)] ,1,mean)
  }
  
  ##  post-cluster-Lasso regression of d_dot on z_all_dot
  
  
  d_dot_long = matrix(d_dot,n*T_t,1) ## t=1, then t=2, ... 
  y_dot_long = matrix(y_dot,n*T_t,1)
  iv_oracle_dot_long =  matrix(iv_oracle_dot,n*T_t,1)
  z_all_dot_long = NULL
  for( i in 1:p_z){
    z_all_dot_long =cbind(z_all_dot_long, matrix(z_all_dot[,((i-1)*T_t+1):(i*T_t)],n*T_t,1)        ) 
  }
  dim(z_all_dot_long)
  
  
  
  #### with non-cluster lasso
  ## decomposed
  # fit.lasso <-rlasso(y=d_dot_long,x=z_all_dot_long)
  # d_dot_hat_long <-d_dot_long - fit.lasso$res 
  # fit.lm <- tsls(x=NULL,d=d_dot_long,y=y_dot_long,z=d_dot_hat_long)
  ## direct way
  fit.lasso <-rlassoIV(x=NULL,d=d_dot_long,y=y_dot_long,z=z_all_dot_long, select.X=FALSE, select.Z=TRUE)
  out <- rbind(out,c(fit.lasso$coefficients[1], fit.lasso$se[1],(fit.lasso$coefficients[1]-alpha)/fit.lasso$se[1]))
  
  # bias_1 = alpha - fit.lm $coefficients[1]
  # se_1 = fit.lm$se[1]
  ##### TSLS 
  fit.tsls <- tsls(x=NULL,d=d_dot_long,y=y_dot_long,z=z_all_dot_long)
  ## Attention: should use HAC covariance matrix for se. 
  out1 <- rbind(out1,c(fit.tsls$coefficients[1], fit.tsls$se[1],(fit.tsls$coefficients[1]-alpha)/fit.tsls$se[1]))
  
  ## oracle
  fit.oracle <-  tsls(x=NULL,d=d_dot_long,y=y_dot_long,z=iv_oracle_dot_long)  
  out2 <- rbind(out2,c(fit.oracle $coefficients[1], fit.oracle $se[1],(fit.oracle $coefficients[1]-alpha)/fit.oracle $se[1]))
  
  ## with cluster lasso
  ## decomposed
  fit.lasso <-rlasso_cluster(y=d_dot_long,x=z_all_dot_long,  T_t= T_t)
  
  # s1 <- sqrt(var(fit.lasso$res))
  e1_t = matrix(fit.lasso$res,n,T_t)
  #Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
  Ups1 <- NULL
  for( k in 1:p_z){
    x_t = matrix(z_all_dot_long[,k],n,T_t)
    x_t = x_t*e1_t
    Ups1 <-c(Ups1, sum(rowSums(x_t^2))/sum(rowSums(x_t)^2)  )
  }
  Ups1 <-  T_t*min(Ups1 )
  Ups1_all <- c(Ups1_all,Ups1)
  
  
  d_dot_hat_long <-d_dot_long - fit.lasso$res 
  fit.lm <- tsls(x=NULL,d=d_dot_long,y=y_dot_long,z=d_dot_hat_long)
  out3 <- rbind(out3,c(fit.lm$coefficients[1], fit.lm$se[1],(fit.lm$coefficients[1]-alpha)/fit.lm$se[1]))
}


Results = matrix(ncol=4, nrow=dim(out1)[1])
Results[,1] = out1[1:dim(out1)[1],1]
Results[,2] = out2[1:dim(out1)[1],1]
Results[,3] = out[1:dim(out1)[1],1]
Results[,4] = out3[1:dim(out1)[1],1]
### COMPUTE BIAS AND RMSE
StatDisplay = data.frame()
StatDisplay[1:4,"bias"] = apply(Results-alpha,2,mean)
StatDisplay[1:4,"RMSE"] = sqrt(apply((Results-alpha)^2,2,mean))
StatDisplay[1:4,"MAD"] = sqrt(apply((Results-alpha)^2,2,median))
row.names(StatDisplay) = c("Naive","Oracle", "LASSO,noncluster","LASSO,cluster")
print(StatDisplay)







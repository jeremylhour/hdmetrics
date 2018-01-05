### Double ML Heterogeneous Treatment Effect



### Set working directory
setwd("//ulysse/users/JL.HOUR/1A_These/B. ENSAE Classes/Cours3A/hdmetrics")

rm(list=ls())
set.seed(03101989)


### 0. Settings

### Load packages
library("ggplot2")
library("gridExtra")
library("MASS")

### Load user-defined functions
source("functions/LassoFISTA.R")
source("functions/LogitLasso.R")
source("functions/ATE.R")

### Pennsylvania Bonus Data
Penn = as.data.frame(read.table("packages/DMLonGitHub-master/penn_jae.dat", header=T ))

########################### Sample Construction ######################

index = (Penn$tg==0) | (Penn$tg==4) # Only keep treatment=4 and controls
data = Penn[index,]
data$tg[data$tg==4] = 1
data[,"dep1"] = data$dep==1
data[,"dep2"] = data$dep==2
data$inuidur1 = log(data$inuidur1)


Y = data[,"inuidur1"] # Outcome
D = data[,"tg"] # Treatment
X = as.matrix(data[,c("female","black","othrace","dep1","dep2","q2","q3","q4","q5","q6","agelt35","agegt54","durable","lusd","husd")]) # Covariates

K = 5
n = length(Y)

X = cbind(rep(1,n),X)
split = runif(n)
cvgroup = as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = T))  

### 1. Lasso-style
theta = vector(length=K)
for(k in 1:K){
  Ik = cvgroup==k # Separate the sample
  NIk = cvgroup!=k
  
  # A. Propensity Score
  treatfit = LogitLasso(D[NIk],X[NIk,],c=.2,PostLasso=T)
  mhat = 1/(1+exp(-X[Ik,] %*% treatfit$betaPL))
  
  # B. Selection on Y_0
  g = .1/log(max(ncol(X),sum(Ik)))
  lambdastar = 2.2*qnorm(1-.5*g/ncol(Ik))/sqrt(sum(Ik)) # Lasso penalty level

  outcomefit0 = LassoFISTA(y=Y[D[NIk]==0],X=X[D[NIk]==0,],nopen=c(1),lambda=0*lambdastar) # Do not penalize the constant

  
}
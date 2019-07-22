### Double ML: Lalonde dataset
### Jeremy LHour
### 24/04/2018


### Set working directory
setwd("//ulysse/users/JL.HOUR/1A_These/B. ENSAE Classes/High-Dimensional Econometrics/hdmetrics")

rm(list=ls())
set.seed(12071990)

### Packages
library("hdm")
library("causalsens")
library("randomForest")

### Load user-defined functions
source("functions/ATE.R")

##########################
##########################
### 0. Lalonde dataset ###
##########################
##########################

### Experimental Estimate
data("lalonde.exp")
summary(lm(re78 ~ treat, data=lalonde.exp))


### min-max scale
mMscale <- function(X){
  X = as.matrix(X)
  mins = apply(X,2,min)
  maxs = apply(X,2,max)
  return(scale(X, center=mins, scale=maxs-mins))
}

# Load data
data("lalonde.psid")
  
D = lalonde.psid[,"treat"]
Y = lalonde.psid[,"re78"]

X.Main = data.frame("Cons"=1,
                      lalonde.psid[,c("age","education","married","black","hispanic","re74","re75","nodegree")],
                      "NoIncome74"=as.numeric(lalonde.psid[,"re74"]==0),
                      "NoIncome75"=as.numeric(lalonde.psid[,"re75"]==0)
  )
  
  

X.Main[,c("age","education","re74","re75")] = mMscale(X.Main[,c("age","education","re74","re75")])
  
X.ContInterac = data.frame(
    "AgexMarried"=lalonde.psid[,"age"]*lalonde.psid[,"married"],
    "AgexNodegree"=lalonde.psid[,"age"]*lalonde.psid[,"nodegree"],
    "AgexBlack"=lalonde.psid[,"age"]*lalonde.psid[,"black"],
    "AgexHispanic"=lalonde.psid[,"age"]*lalonde.psid[,"hispanic"],
    "AgexNoIncome74"=lalonde.psid[,"age"]*(lalonde.psid[,"re74"]==0),
    "AgexNoIncome75"=lalonde.psid[,"age"]*(lalonde.psid[,"re75"]==0),
    "EducxMarried"=lalonde.psid[,"education"]*lalonde.psid[,"married"],
    "EdcuxNodegree"=lalonde.psid[,"education"]*lalonde.psid[,"nodegree"],
    "EdcuxBlack"=lalonde.psid[,"education"]*lalonde.psid[,"black"],
    "EdcuxHispanic"=lalonde.psid[,"education"]*lalonde.psid[,"hispanic"],
    "EdcuxNoIncome74"=lalonde.psid[,"education"]*(lalonde.psid[,"re74"]==0),
    "EdcuxNoIncome75"=lalonde.psid[,"education"]*(lalonde.psid[,"re75"]==0),
    "Income74xMarried"=lalonde.psid[,"re74"]*lalonde.psid[,"married"],
    "Income74xNoDegree"=lalonde.psid[,"re74"]*lalonde.psid[,"nodegree"],
    "Income74xBlack"=lalonde.psid[,"re74"]*lalonde.psid[,"black"],
    "Income74xHispanic"=lalonde.psid[,"re74"]*lalonde.psid[,"hispanic"],
    "Income74xNoIncome75"=lalonde.psid[,"re74"]*(lalonde.psid[,"re75"]==0),
    "Income75xMarried"=lalonde.psid[,"re75"]*lalonde.psid[,"married"],
    "Income75xNodegree"=lalonde.psid[,"re75"]*lalonde.psid[,"nodegree"],
    "Income75xBlack"=lalonde.psid[,"re75"]*lalonde.psid[,"black"],
    "Income75xHispanic"=lalonde.psid[,"re75"]*lalonde.psid[,"hispanic"],
    "Income75xNoIncome74"=lalonde.psid[,"re75"]*(lalonde.psid[,"re74"]==0)
  )
  
X.ContInterac = mMscale(X.ContInterac)
  
X.DumInterac = data.frame(
    "MarriedxNodegree"=lalonde.psid[,"married"]*lalonde.psid[,"nodegree"],
    "MarriedxBlack"=lalonde.psid[,"married"]*lalonde.psid[,"black"],
    "MarriedxHispanic"=lalonde.psid[,"married"]*lalonde.psid[,"hispanic"],
    "MarriedxNoIncome74"=lalonde.psid[,"married"]*(lalonde.psid[,"re74"]==0),
    "MarriedxNoIncome75"=lalonde.psid[,"married"]*(lalonde.psid[,"re75"]==0),
    "NodegreexBlack"=lalonde.psid[,"nodegree"]*lalonde.psid[,"black"],
    "NodegreexHispanic"=lalonde.psid[,"nodegree"]*lalonde.psid[,"hispanic"],
    "NodegreexNoIncome74"=lalonde.psid[,"nodegree"]*(lalonde.psid[,"re74"]==0),
    "NodegreexNoIncome75"=lalonde.psid[,"nodegree"]*(lalonde.psid[,"re75"]==0),
    "BlackxNoIncome74"=lalonde.psid[,"black"]*(lalonde.psid[,"re74"]==0),
    "BlackxNoIncome75"=lalonde.psid[,"black"]*(lalonde.psid[,"re75"]==0),
    "HispanicxNoIncome74"=lalonde.psid[,"hispanic"]*(lalonde.psid[,"re74"]==0),
    "HispanicxNoIncome75"=lalonde.psid[,"hispanic"]*(lalonde.psid[,"re75"]==0),
    "NoIncome74x75"=(lalonde.psid[,"re74"]==0)*(lalonde.psid[,"re75"]==0)
  )
  
X.Poly = poly(as.matrix(lalonde.psid[,c("age","education","re74","re75")]),degree=5)
X.Poly = mMscale(X.Poly)
X.Age  = model.matrix(~as.factor(lalonde.psid[,"age"]) - 1)[,-1]
  
X = cbind(X.Main,X.ContInterac,X.DumInterac,X.Poly)
X = as.matrix(X)

K = 5
n = length(Y)


###############################################
###############################################
### 1. Simple DML w Lasso and Random Forest ###
###############################################
###############################################

split = runif(n)
cvgroup = as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = T))

thetaL = vector(length=K); gammaL = vector(length=n)
thetaRF = vector(length=K); gammaRF = vector(length=n)

for(k in 1:K){
  print((paste("Iteration ",k)))
  Ik = cvgroup==k # Split the sample
  NIk = cvgroup!=k
  
  # A. Propensity Score
  treatfit = rlassologit(D[NIk] ~ X[NIk,]) # Logit Lasso
  mhatL = predict(treatfit,newdata=X[Ik,])
  IndL = (mhatL>0) & (mhatL<1)
  
  optRF = tuneRF(X[NIk,],as.factor(D[NIk]),stepFactor=1.5, improve=0.05, nodesize=5, doBest=T, plot=F, trace=F)
  min = optRF$mtry
  treatForest = randomForest(X[NIk,],as.factor(D[NIk]), mtry=min)
  mhatRF = predict(treatForest, X[Ik,], type="prob")[,2]
  IndRF = (mhatRF>0) & (mhatRF<1)
  
  # B. E(Y_0 \vert X)
  outcomefit0 = rlasso(Y[D[NIk]==0] ~ X[D[NIk]==0,]) # Lasso
  g0hatL = predict(outcomefit0,newdata = X[Ik,])
  
  optRF = tuneRF(X[D[NIk]==0,],Y[D[NIk]==0],stepFactor=1.5, improve=0.05, nodesize=5, doBest=T, plot=F, trace=F)
  min = optRF$mtry
  outcomeForest0 = randomForest(X[D[NIk]==0,],Y[D[NIk]==0], mtry=min)
  g0hatRF = predict(outcomeForest0, X[Ik,], type="response")
  
  # C. ATT
  thetaL[k] = ATT(Y[Ik],D[Ik],g0hatL,mhatL,IndL)
  gammaL[Ik] = SE.ATT(Y[Ik],D[Ik],g0hatL,mhatL,IndL)$gamma
  
  thetaRF[k] = ATT(Y[Ik],D[Ik],g0hatRF,mhatRF,IndRF)
  gammaRF[Ik] = SE.ATT(Y[Ik],D[Ik],g0hatRF,mhatRF,IndL)$gamma
}

### Results

## A. Lasso
SE.thetaL = sd(gammaL)/sqrt(n)

print("ATT, Point estimate")
print(mean(thetaL))
print("Standard Error")
print(SE.thetaL)
print(".95 Confidence Interval")
print(paste("[",mean(thetaL)+qnorm(.025)*SE.thetaL,";",mean(thetaL)+qnorm(.975)*SE.thetaL,"]"))

### B. Random Forest
SE.thetaRF = sd(gammaRF, na.rm=T)/sqrt(n)

print("ATT, Point estimate")
print(mean(thetaRF))
print("Standard Error")
print(SE.thetaRF)
print(".95 Confidence Interval")
print(paste("[",mean(thetaRF)+qnorm(.025)*SE.thetaRF,";",mean(thetaRF)+qnorm(.975)*SE.thetaRF,"]"))


########################################################
########################################################
### 2. Median to mitigate impact of sample-splitting ###
########################################################
########################################################

# cf. DML paper, definition 3.3, p. 30
# Here the partition of the sample is recomputed S times
S = 20

thetaL = matrix(nrow=K,ncol=S); gammaL = matrix(nrow=n,ncol=S)
thetaRF = matrix(nrow=K,ncol=S); gammaRF = matrix(nrow=n,ncol=S)

for(s in 1:S){
  split = runif(n)
  cvgroup = as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = T))
  
  for(k in 1:K){
    Ik = cvgroup==k # Split the sample
    NIk = cvgroup!=k
    
    # A. Propensity Score
    treatfit = rlassologit(D[NIk] ~ X[NIk,]) # Logit Lasso
    mhatL = predict(treatfit,newdata=X[Ik,])
    IndL = (mhatL>0) & (mhatL<1)
    
    optRF = tuneRF(X[NIk,],as.factor(D[NIk]),stepFactor=1.5, improve=0.05, nodesize=5, doBest=T, plot=F, trace=F)
    min = optRF$mtry
    treatForest = randomForest(X[NIk,],as.factor(D[NIk]), mtry=min)
    mhatRF = predict(treatForest, X[Ik,], type="prob")[,2]
    IndRF = (mhatRF>0) & (mhatRF<1)
    
    # B. E(Y_0 \vert X)
    outcomefit0 = rlasso(Y[D[NIk]==0] ~ X[D[NIk]==0,]) # Lasso
    g0hatL = predict(outcomefit0,newdata = X[Ik,])
    
    optRF = tuneRF(X[D[NIk]==0,],Y[D[NIk]==0],stepFactor=1.5, improve=0.05, nodesize=5, doBest=T, plot=F, trace=F)
    min = optRF$mtry
    outcomeForest0 = randomForest(X[D[NIk]==0,],Y[D[NIk]==0], mtry=min)
    g0hatRF = predict(outcomeForest0, X[Ik,], type="response")
    
    # C. ATT
    nIk = Ik
    nIk[Ik] = IndL
    
    thetaL[k,s] = ATT(Y[Ik],D[Ik],g0hatL,mhatL,IndL)
    gammaL[nIk,s] = SE.ATT(Y[Ik],D[Ik],g0hatL,mhatL,IndL)$gamma
    
    nIk = Ik
    nIk[Ik] = IndRF
    
    thetaRF[k,s] = ATT(Y[Ik],D[Ik],g0hatRF,mhatRF,IndRF)
    gammaRF[nIk,s] = SE.ATT(Y[Ik],D[Ik],g0hatRF,mhatRF,IndRF)$gamma
  }
}

### Results based on the median

## A. Lasso
thetaL.DML = apply(thetaL,2,mean)
thetaL.med = median(thetaL.DML)
SE.thetaL = apply(gammaL,2,sd,na.rm=T)/sqrt(n)

print("ATT, Point estimate")
print(thetaL.med)
print("Standard Error")
print(sqrt(median(SE.thetaL^2 +  (thetaL.DML-thetaL.med)^2)))

### B. Random Forest
thetaRF.DML = apply(thetaRF,2,mean)
thetaRF.med = median(thetaRF.DML)
SE.thetaRF = apply(gammaRF,2,sd,na.rm=T)/sqrt(n)

print("ATT, Point estimate")
print(thetaRF.med)
print("Standard Error")
print(sqrt(median(SE.thetaRF^2 +  (thetaRF.DML-thetaRF.med)^2)))
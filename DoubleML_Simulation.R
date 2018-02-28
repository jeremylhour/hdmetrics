### Double Selection and sample-splitting, simulations
### Used in Section 2.2 of the HDMetrics class.
### Jeremy L Hour
### 04/01/2018
### Edited: 28/02/2018

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
source("functions/DataSim.R") 
source("functions/LassoFISTA.R")

### Simulation parameters
R = 10000 # nb simulations
n = 200 # sample size
p = 150 # nb variables
K = 5 # nb folds
a = .5 # ATT

split = runif(n)
cvgroup = as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = T))  

g = .1/log(max(p,n))
lambda = 2.2*qnorm(1-.5*g/p)/sqrt(n) # Lasso penalty level

Results = matrix(ncol=3, nrow=R)
stdev = vector(length=R)
t_start = Sys.time()
pb = txtProgressBar(style = 3)

for(r in 1:R){
  ### GENERATE DATA
  data = DataSim(n=n,p=p,Ry=.3,Rd=.7)
  X = data$X; y = data$y; d = data$d
  
  ### METHOD 1: Naive selection
  lassoselec = LassoFISTA(y=y,X=cbind(d,X),nopen=c(1,2),lambda=lambda) # Do not penalize the constant
  Snaive = which(lassoselec$beta != 0)
  Snaive = Snaive[!Snaive %in% c(1,2)] # delete intercept and treatment variable
  if(length(Snaive)==0){
    naivefit = lm(y ~ d)
  } else {
    naivefit = lm(y ~ d + X[,Snaive])
  }
  
  
  ### METHOD 2: Double-Selection, no sample-splitting
  # A. Selection on Treatment
  treatfit = LassoFISTA(y=d,X=X,nopen=c(1),lambda=.15*lambda) # Do not penalize the constant
  Sd = which(treatfit$beta != 0)
  Sd = Sd[!Sd == 1] # delete intercept
  
  # B. Selection on Outcome
  outcomefit = LassoFISTA(y=y,X=X,nopen=c(1),lambda=lambda) # Do not penalize the constant
  Sy = which(outcomefit$beta != 0)
  Sy = Sy[!Sy == 1] # delete intercept
  
  # C. Compute Post-Double-Selection
  Shat = union(Sy,Sd)
  if(length(Shat)==0){
    DSfit = lm(y ~ d)
  } else {
    DSfit = lm(y ~ d + X[,Shat])
  }
  
  # D. Compute sd
  if(length(Shat)==0){
    treatfit = lm(d ~ 1)
  } else {
    treatfit = lm(d ~ X[,Shat])
  }
  sigmaNum = sum(treatfit$residuals^2*DSfit$residuals^2) /(n - length(Shat) - 1)
  sigmaDenom = sum(treatfit$residuals^2) / n
  stdev[r] = sqrt( sigmaNum / sigmaDenom^2) / sqrt(n)
  
  
  ### METHOD 3: Double Selection with Sample Splitting
  theta = vector(length=K)
  for(k in 1:K){
    Ik = cvgroup==k # Separate the sample
    NIk = cvgroup!=k
    
    # 0. Adjust Lasso penalty level
    gstar = .1/log(max(p,sum(Ik)))
    lambdastar = 2.2*qnorm(1-.5*g/p)/sqrt(sum(Ik)) # Lasso penalty level
    
    # Abis. Selection on Treatment
    treatfit = LassoFISTA(y=d[NIk],X=X[NIk,],nopen=c(1),lambda=.15*lambdastar) # Do not penalize the constant
    Sd = which(treatfit$beta != 0)
    Sd = Sd[!Sd == 1] # delete intercept
    
    # Bbis. Selection on Outcome
    outcomefit = LassoFISTA(y=y[NIk],X=X[NIk,],nopen=c(1),lambda=lambdastar) # Do not penalize the constant
    Sy = which(outcomefit$beta != 0)
    Sy = Sy[!Sy == 1] # delete intercept
    
    # Cbis. Compute Post-Double-Selection
    Shat = union(Sy,Sd)
    if(length(Shat)==0){
      outcomePL = lm(y[NIk] ~ 1)
      treatPL = lm(d[NIk] ~ 1)
    } else {
      outcomePL = lm(y[NIk] ~ X[NIk,Shat])
      treatPL = lm(d[NIk] ~ X[NIk,Shat])
    }
    
    # D. Target param on left-out sample
    ytilde = y[Ik] - cbind(rep(1,sum(Ik)),X[Ik,Shat])%*%coef(outcomePL)
    dtilde = d[Ik] - cbind(rep(1,sum(Ik)),X[Ik,Shat])%*%coef(treatPL)
    Ikfit = lm(ytilde ~ dtilde)
    
    theta[k] = Ikfit$coef['dtilde']
  }
  
  
  ### COLLECTING RESULTS
  Results[r,] = c(naivefit$coef['d'],
                  DSfit$coef['d'],
                  mean(theta))
  
  setTxtProgressBar(pb, r/R)
}

close(pb)
print(Sys.time()-t_start)

### COMPUTE BIAS AND RMSE
StatDisplay = data.frame()
StatDisplay[1:3,"bias"] = apply(Results-a,2,mean)
StatDisplay[1:3,"RMSE"] = sqrt(apply((Results-a)^2,2,mean))
row.names(StatDisplay) = c("Naive","Immunized","Immunized, Cross-fitted")
print(StatDisplay)

### DRAW CHARTS
id = c(mapply(function(x) rep(x,R),1:3))
val = c(Results)-a
data_res = data.frame(val = val, model = id)

M = max(abs(quantile(Results,.01,na.rm=T)),abs(quantile(Results,.99,na.rm=T)))
lb = -1.1*M; ub = 1.1*M

get.plot <- function(data,modelS,title="A Title",s){
  plot_res <- ggplot(subset(data, (model==modelS)), aes(x=val)) + 
    geom_histogram(binwidth = .02, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
    ggtitle(title) + 
    stat_function(fun = dnorm, args=list(mean=0, sd=s), colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  return(plot_res)
} # plot func

pdf("plots/Immunized.pdf",width=14,height=4)
grid.arrange(get.plot(data_res,1,"Naive Post-Selec", mean(stdev)), get.plot(data_res,2,"Double-Selec", mean(stdev)), get.plot(data_res,3,"Double-Selec, Cross-fitting", mean(stdev)), ncol=3)
dev.off()


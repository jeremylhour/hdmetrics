### Double Selection Procedure
### Jeremy L Hour
### 06/12/2017

### Set working directory
setwd("//ulysse/users/JL.HOUR/1A_These/B. ENSAE Classes/Cours3A/Code")

rm(list=ls())
set.seed(30031987)


### 0. Settings

### Load packages
library("ggplot2")
library("gridExtra")
library("MASS")

### Load user-defined functions
source("functions/DataSim.R") 
source("functions/LassoFISTA.R")


### 1. Simulations
R = 500
n = 200
S = cbind(1:(n/2),(n/2+1):n)
p = 50
g = .1/log(max(p,n))
lambda = 2.2*qnorm(1-.5*g/p)/sqrt(n)

Results = matrix(ncol=3, nrow=R)
t_start = Sys.time()
pb = txtProgressBar(style = 3)

for(r in 1:R){
  ### 1. Generate data
  data = DataSim(n=n,p=p,Ry=.1,Rd=.5)
  X = data$X; y = data$y; d = data$d
  
  ### 2. Selection on outcome
  outcomefit = LassoFISTA(y=y,X=X,nopen=c(1),lambda=lambda) # Do not penalize the constant
  Sy = which(outcomefit$beta != 0)
  
  ### 2. Selection on treatment
  treatfit = LassoFISTA(y=d,X=X,nopen=c(1),lambda=.15*lambda) # Do not penalize the constant
  Sd = which(treatfit$beta != 0)
  
  ### 3. Compute naive plug-in
  naivefit = lm(y ~ d + X[,Sy]-1)
  
  ### 4. Compute Post-Double-Selection
  Shat = union(Sy,Sd)
  DSfit = lm(y ~ d + X[,Shat]-1)
  
  ### 4.5 Compute standard error (not useful here)
  treatfit = lm(d ~ X[,Shat]-1)
  sigmaNum = sum(treatfit$residuals^2*DSfit$residuals^2) /(n - length(Shat) - 1)
  sigmaDenom = sum(treatfit$residuals^2) / n
  DSsigma = sqrt( sigmaNum / sigmaDenom^2) / sqrt(n)
  
  ### 5. Double ML
  theta = vector(length=2) # ici pb car le 'n' n'est pas le meme
  for(i in 1:2){
    ### Compute nuisance param on Sample 1
    outcomefit = LassoFISTA(y=y[S[,i]],X=X[S[,i],],nopen=c(1),lambda=lambda)
    Sy = which(outcomefit$beta != 0)
    
    treatfit = LassoFISTA(y=d[S[,i]],X=X[S[,i],],nopen=c(1),lambda=.15*lambda)
    Sd = which(treatfit$beta != 0)
    
    Shat = union(Sy,Sd)
    outcomePL = lm(y[S[,i]] ~ X[S[,i],Shat]-1)
    treatPL = lm(d[S[,i]] ~ X[S[,i],Shat]-1)
    
    ### Compute param of interest on Sample 2
    ytilde = y[S[,3-i]] - as.matrix(X[S[,3-i],Shat])%*%coef(outcomePL)
    dtilde = d[S[,3-i]] - as.matrix(X[S[,3-i],Shat])%*%coef(treatPL)
    oosfit = lm(ytilde ~ dtilde)
    
    theta[i] = oosfit$coef['dtilde']
  }
  
  ### 6. Third step: ATT estimation
  Results[r,] = c(naivefit$coef['d'],
                  DSfit$coef['d'],
                  mean(theta))
  
  setTxtProgressBar(pb, r/R)
}

close(pb)
print(Sys.time()-t_start)

# Draw the charts
id = c(mapply(function(x) rep(x,R),1:3))
val = c(Results)
data_res = data.frame(val = val, model = id)

M = max(abs(quantile(Results,.01,na.rm=T)),abs(quantile(Results,.99,na.rm=T)))
lb = -1.1*M
ub = 1.1*M
sdBCH = DSsigma

### Function for plot
get.plot <- function(data,modelS,title="A Title",sdBCH){
  plot_res <- ggplot(subset(data, (model==modelS)), aes(x=val)) + 
    geom_histogram(binwidth = .02, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
    ggtitle(title) + 
    stat_function(fun = dnorm, args=list(mean=.5, sd=sdBCH), colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  return(plot_res)
}


grid.arrange(get.plot(data_res,1,"Naive Plug-In", sdBCH), get.plot(data_res,2,"Immunized", sdBCH), ncol=2)


### Compute bias and RMSE
StatDisplay = data.frame()
StatDisplay[1:3,"bias"] = apply(Results-.5,2,mean)
StatDisplay[1:3,"RMSE"] = sqrt(apply((Results-.5)^2,2,mean))
row.names(StatDisplay) = c("NaivePlugIn","Immunized","DoubleML")
print(StatDisplay)
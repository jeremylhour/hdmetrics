###### # application IV eminent domain

###### Simulations from  Chern. Hansen. Spindler. 
### 0. Settings

### Load packages
library("ggplot2")
library("gridExtra")
library("MASS")
library("mnormt")
library(hdm)
library(AER)
### Simulation parameters
set.seed(13571113)
p_x = 200 ## number of controls
p_z = 150 ## number of instruments 
n = 202 ## total sample size
K = 2 # nb folds

#### Splitting decision rules
split = runif(n)
cvgroup = as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = T))  

##number of MC replications
MC =1000
Results = matrix(ncol=3, nrow=MC)

out = NULL
out1 = NULL
out2 = NULL
out3=NULL
outKK1 = NULL
outKK2 = NULL
for (kk in 1:MC){
  
  ### GENERATE DATA
  means <- c(0,0,0,0)
  Sigma <- matrix(0,p_x,p_x)
  for (i in 1:p_x){
    for (j in 1:p_x){
      Sigma[i,j] <- (1/2)^{abs(i-j)}
    }
  }
  
  nu <- 4/9 + sum((1:p_x)^(-2))
  beta <- matrix(0,1,p_x)
  beta[1,1:4] <-  1/(9*nu)
  beta[1,5:p_x] <-  1/nu*(5:p_x)^(-2)
  
  delta <- matrix(3*(1:p_z)^(-2),p_z,1)
  Pi_m <- cbind(diag(1,p_z,p_z),matrix(0,p_z,(p_x-p_z)))
  
  sigmas <- matrix(0,2+p_z+p_x, 2+p_z+p_x)
  sigmas[1:2,1:2] <- matrix(c(1,0.6,0.6,1), 2,2  )
  sigmas[3:(2+p_z),3:(2+p_z)] <- diag(1,p_z,p_z)
  sigmas[(3+p_z):(2+p_z+p_x),(3+p_z):(2+p_z+p_x)] <-Sigma
  
  
  var <- rmnorm(n, mean = rep(0, nrow(sigmas)), varcov = sigmas)
  dim(var )
  eps <- var[,1]
  us <- var[,2]
  zetas <- var[,3:(2+p_z)]
  x <- var[,(3+p_z):(2+p_z+p_x)]
  gamma = beta
  tau = 1.5
  z <- Pi_m%*%t(x) + 0.125*t(zetas)
  d <- x%*%t(gamma) +  t(z)%*%delta + us
  y <- tau *d + x%*%t(beta) + 2*eps
  z <- t(z)
  
  ### METHOD 0bis: oracle
  zO =z%*%delta
  xO =x%*%t(beta) 
  ivfit.lasso = tsls(y=y,d=d, x=xO, z=zO)
  out3 <- rbind(out3,c(ivfit.lasso$coef[1], ivfit.lasso$se[1],ivfit.lasso$coef[1]/ivfit.lasso$se[1]))
  
  ### METHOD 1: Double-Selection, no sample-splitting
  ## Do LASSO of D on X to obtain gamma
  rD_x = rlasso(d ~ x + z)
  ## Do LASSO of Y on X to obtain theta, and extract residuals
  rY_x = rlasso(y ~ x)
  rY = rY_x$res
  ## Build D_hat from estimated gamma and delta
  dhat <- (d- rD_x$res) 
  ## regress d_hat on X to get nu
  rD = rlasso(dhat ~ x)
  ## extract the residuals
  rD_res =rD$res
  rD1 <- d - (dhat -rD$res ) 
  ## Do TSLS of the residuals of Y/X on residuals of D/X using residuals of Dhat/X as instruments
  ivfit.lasso = tsls(y=rY,d=rD1, x=NULL, z=rD_res, intercept = FALSE)
  out <- rbind(out,c(ivfit.lasso$coef[1], ivfit.lasso$se[1],ivfit.lasso$coef[1]/ivfit.lasso$se[1]))
  
  ### Build in function to do all this.... 
 # ivfit.lasso2 = rlassoIV(x,d,y,z, select.X=TRUE, select.Z=TRUE)
 #  out1 <- rbind(out1,c(ivfit.lasso2$coef, ivfit.lasso2$se,ivfit.lasso2$coef/ivfit.lasso2$se ))

  ### METHOD 0: selection, alternative (Non-orthogonal)
  ## select all the controls selected by the two Lasso
  sel = (abs(rD_x$coefficients[2:(dim(x)[2]+1)])> 10^(-6))*1 + (rY_x$coefficients[2:(dim(x)[2]+1)]> 10^(-6))*1
  sel[sel ==2] <- 1 
  sel_z = (rD_x$coefficients[(2+dim(x)[2]):(1+dim(x)[2]+dim(z)[2])] > 10^(-6))*1 
  ## Do TSLS 
  x_sel = x[,sel==1]
  z_sel = z[,sel_z==1]
  if(sum(sel)>0 & sum(sel_z)>0){
     ivfit.lm = ivreg(y ~ d  + x_sel| z_sel + x_sel)
   }else if (sum(sel)==0 & sum(sel_z)>0){
     ivfit.lm = ivreg(y ~ d  | z_sel)
   }
  se <-  coef(summary(ivfit.lm))[2, "Std. Error"]
  out1 <- rbind(out1,c(ivfit.lm$coef["d"],  se ,ivfit.lm$coef["d"]/se))
  
 
  
  
 ### METHOD 2: Double Selection with Sample Splitting
 outK = matrix(ncol=3, nrow=K)
 for(k in 1:K){
   Ik = cvgroup==k # Separate the sample
   NIk = cvgroup!=k
   ind <- matrix(1,dim(d)[1],1)
   ind_x <- matrix(1,dim(x[Ik,])[1],1)
   ## Do LASSO of D on X to obtain gamma
   rD_x = rlasso(d[NIk,] ~ x[NIk,] + z[NIk,])
   ## Do LASSO of Y on X to obtain theta, and extract residuals
   rY = rlasso(y[NIk,] ~ x[NIk,])
   ## Build D_hat from estimated gamma and delta
   dhat <-  cbind( ind ,x,z)%*%rD_x$coefficients
   ## regress d_hat on X to get nu
   rD = rlasso(dhat[NIk,] ~ x[NIk,])
   ## extract the residuals
   rD_res = dhat[Ik,] - cbind( ind_x ,x[Ik,])%*%rD$coefficients
   rD1 <- d[Ik,] - cbind( ind_x ,x[Ik,])%*%rD$coefficients
   rY <-  y[Ik,] -  cbind( ind_x ,x[Ik,])%*%rY$coefficients
   ## Do TSLS 
   ivfit.lasso = tsls(y=rY,d=rD1, x=NULL, z=rD_res, intercept = FALSE)
   outK[k,] <- c(ivfit.lasso$coef[1], ivfit.lasso$se[1],ivfit.lasso$coef[1]/ivfit.lasso$se[1])
  }
  outK1 <- outK
  coef1 <- median( outK1[,1])
  outK1[,2] <- outK1[,2] +(  outK1[,1] -   coef1  )^2
  outKK1 <-  rbind(outKK1,c(coef1, median(  outK1[,2]),coef1/median(  outK1[,2])))
 
  outK1 <- outK
  coef1 <- mean( outK1[,1])
  outK1[,2] <- outK1[,2] +(  outK1[,1] -   coef1  )^2
  outKK2 <-  rbind(outKK2,c(coef1, mean(  outK1[,2]),coef1/mean(  outK1[,2])))
 
 
}



hist(out[,3]-tau/out[,2],100, prob=T, col=4)
xseq = seq(-3,3,length.out=100)
lines(xseq,dnorm(xseq),col=2, lwd=2)

hist(out1[,3]-tau/out1[,2],100, prob=T, col=4)
xseq = seq(-3,3,length.out=100)
lines(xseq,dnorm(xseq),col=2, lwd=2)

hist(out3[,3]-tau/out3[,2],100, prob=T, col=4)
xseq = seq(-3,3,length.out=100)
lines(xseq,dnorm(xseq),col=2, lwd=2)

# hist(out[,3],100, prob=T, col=4)
# xseq = seq(-3,3,length.out=100)
# lines(xseq,dnorm(xseq),col=2, lwd=2)
# 
# hist(out1[,3],100, prob=T, col=4)
# xseq = seq(-3,3,length.out=100)
# lines(xseq,dnorm(xseq),col=2, lwd=2)
# 
# hist(out3[,3],100, prob=T, col=4)
# xseq = seq(-3,3,length.out=100)
# lines(xseq,dnorm(xseq),col=2, lwd=2)
# 
# 
# hist(outKK1[,3],100, prob=T, col=4)
# xseq = seq(-3,3,length.out=100)
# lines(xseq,dnorm(xseq),col=2, lwd=2)
# 
# hist(outKK2[,3],100, prob=T, col=4)
# xseq = seq(-3,3,length.out=100)
# lines(xseq,dnorm(xseq),col=2, lwd=2)

Results = matrix(ncol=5, nrow=dim(outKK2)[1])
Results[,2] = out[1:dim(outKK2)[1],1]
Results[,1] = out1[1:dim(outKK2)[1],1]
Results[,3] = out3[1:dim(outKK2)[1],1]
Results[,4] = outKK1[1:dim(outKK2)[1],1]
Results[,5] = outKK2[1:dim(outKK2)[1],1]
### COMPUTE BIAS AND RMSE
StatDisplay = data.frame()
StatDisplay[1:5,"bias"] = apply(Results-tau,2,mean)
StatDisplay[1:5,"RMSE"] = sqrt(apply((Results-tau)^2,2,mean))
StatDisplay[1:5,"MAD"] = sqrt(apply((Results-tau)^2,2,median))
row.names(StatDisplay) = c("Naive","Immunized","Oracle", "Cross-fitted med.","Cross-fitted mean.")
print(StatDisplay)


library(xtable)
xtable(StatDisplay)

Results[,1] = out1[1:dim(outKK2)[1],3] - tau/out1[1:dim(outKK2)[1],2] 
Results[,2] = out[1:dim(outKK2)[1],3]- tau/out[1:dim(outKK2)[1],2] 
Results[,3] = out3[1:dim(outKK2)[1],3]- tau/out3[1:dim(outKK2)[1],2] 
Results[,4] = outKK1[1:dim(outKK2)[1],3]- tau/outKK1[1:dim(outKK2)[1],2] 
Results[,5] = outKK2[1:dim(outKK2)[1],3]- tau/outKK2[1:dim(outKK2)[1],2] 



Results_s <- Results
### DRAW CHARTS
id = c(mapply(function(x) rep(x,dim(outKK2)[1]),1:5))
val = c(Results)
data_res = data.frame(val = val, model = id)
length(id)

M = max(abs(quantile(Results,.05,na.rm=T)),abs(quantile(Results,.95,na.rm=T)))
lb = -4; ub = 4

get.plot <- function(data,modelS,title="A Title",s){
  plot_res <- ggplot(subset(data, (model==modelS)), aes(x=val)) + 
    geom_histogram(binwidth = .02, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
    ggtitle(title) + 
    stat_function(fun = dnorm, args=list(mean=0, sd=s), colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  return(plot_res)
} # plot func

# pdf("plots/Immunized.pdf",width=14,height=4)
stdev=1
grid.arrange(get.plot(data_res,1,"Naive", 1), get.plot(data_res,2,"Immunized",1), get.plot(data_res,3,"Oracle",1), ncol=3)

grid.arrange(get.plot(data_res,4,"Cross-fitted med.",1),get.plot(data_res,5,"Cross-fitted mean", 1), ncol=2)


###### ##### 



library(BLPestimatoR)



library(xtable)


table(out)





###### 
install.packages("Rmosek", type="source")


library(devtools)
install_github("cran/Rmosek")
install_github("cran/hdm")

libra
install.packages('C:/Users/gaillac/Downloads/hdm_0.2.3.tar.gz', repos = NULL, type="source")
install.packages('C:/Users/gaillac/Downloads/Rmosek_1.2.5.1.tar.gz', repos = NULL, type="source")

install.packages("gpuR")

install_github("cdeterman/gpuR")
devtools::install_github("cdeterman/RViennaCL")
library("gpuR")
# verify you have valid GPUs
detectGPUs()
set.seed(123)
gpuA <- gpuMatrix(rnorm(16), nrow=4, ncol=4)
gpuB <- gpuA %*% gpuA

















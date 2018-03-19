rm(list=ls())
library(hdm)
library(foreign)

## import dataset
path = "C:/Users/gaillac/Documents/GitHub/hdmetrics/data/NEW7080.dta"
data <- read.dta(path )
head(data)
## rename columns
colnames(data) <- c("AGE","AGEQ","v3","EDUC","ENOCENT","ESOCENT","v7","v8","LWKLYWGE","MARRIED",
                    "MIDATL","MT","NEWENG","v14","v15", "CENSUS","v17","QOB","RACE","SMSA","SOATL","v22","v23",
                    "WNOCENT","WSOCENT","v26","YOB")
# LWKLYWGE: Log weekly earnings
data1 <- subset( data,CENSUS==80 & v17 <=50 & YOB <=40)
head(data1)
rm(data)

### generate dummines for Year of birth, Quarter of birth, and interactions with the state
YOB.f = factor(data1[,"YOB"] )
v7.f = factor(data1[,"v17"] )
x  = model.matrix(~YOB.f + v7.f +YOB.f*v7.f)
QOB.f = factor(data1[,"QOB"] )
z  = model.matrix(~QOB.f +YOB.f + v7.f  + QOB.f*v7.f + QOB.f*YOB.f)

z  = z[,-1]
head(z)
cont <- colnames(x)
iv <- colnames(z)

####### Baseline 2SLS with 3 IV
Result = matrix(ncol=2, nrow=3)
form <- paste("LWKLYWGE", paste(c(cont,"EDUC"), collapse=" + "), sep=" ~ ")
fit.tsls.b <- tsls(x=NULL,d=as.matrix(data1[,"EDUC"]),y=as.matrix(data1[,"LWKLYWGE"]),z=z[,1:3])
summary(fit.tsls.b)
Result[1,] <-  c(fit.tsls.b$coefficients["d1",], fit.tsls.b $se["d1"])

####### Baseline 2SLS with 180 IV
fit.tsls.b <- tsls(x=NULL,d=as.matrix(data1[,"EDUC"]),y=as.matrix(data1[,"LWKLYWGE"]),z=z)
summary(fit.tsls.b)
Result[2,] <-  c(fit.tsls.b$coefficients["d1",], fit.tsls.b $se["d1"])

## Baseline 2SLS Selection
d = as.matrix(data1[,"EDUC"])
y = as.matrix(data1[,"LWKLYWGE"])
rD_x = rlasso(d ~ x + z)
## Do LASSO of Y on X to obtain theta, and extract residuals
rY_x = rlasso(y ~ x)
rY = rY_x$res
rm(rY_x)
## Build D_hat from estimated gamma and delta
dhat <- (d- rD_x$res) 
## regress d_hat on X to get nu
rD = rlasso(dhat ~ x)
## extract the residuals
rD_res =rD$res
rD1 <- d - (dhat -rD$res ) 
rm(rD_x)
rm(rD)
## Do TSLS 
ivfit.lasso = tsls(y=rY,d=rD1, x=NULL, z=rD_res, intercept = FALSE)
summary(ivfit.lasso)
# fit.lasso.b <-rlassoIV(x=x,d=as.matrix(data1[,"EDUC"]),y=as.matrix(data1[,"LWKLYWGE"]),z=z, select.X=TRUE, select.Z=TRUE)
Result[3,] <- c(ivfit.lasso$coefficients,ivfit.lasso$se)

library(ivmodel)
##### for Fuller corrected estimates 
iv1 = ivmodel(Y=y, D=d, Z=z, X=z)
summary(iv1)
Fuller(iv1)


iv2 = ivmodel(Y=y, D=d, Z=z[,1:3], X=z[,1:3])
summary(iv2)
Fuller(iv2)
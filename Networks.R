### Network simulations
### Jeremy L Hour
### 26/01/2018


### Set working directory
setwd("//ulysse/users/JL.HOUR/1A_These/B. ENSAE Classes/Cours3A/hdmetrics")

rm(list=ls())
set.seed(12071990)

### Load packages
library("MASS")

### 1. Simple example from Bramoulle
n = 10
beta = .5
eta = .2
gamma = .3
W = matrix(0, nrow=n, ncol=n)
for(i in 2:n) W[i,i-1] = 1

X = rnorm(n,0,sd=.7)

### Residual variance-covariance matrix
rho = .5
Sigma = matrix(0,nrow=n, ncol=n)

for(k in 1:n){
  for(j in 1:n){
    Sigma[k,j] = rho^abs(k-j)
  }
}
eps = mvrnorm(n = 1, mu=rep(0,n), Sigma)

M = solve(diag(n)-beta*W)
y = M %*% (eta*diag(n) + gamma*W) %*% X + M %*% eps

### Regression
peery = W%*%y
peerX = W%*%X
simplereg = lm(y ~ peery + peerX + X)
summary(simplereg)
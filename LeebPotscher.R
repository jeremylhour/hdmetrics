### Post-selection estimator distribution
### From Leeb and Potscher (2006)
### 29 decembre 2017
### JL

setwd("/Users/jeremylhour/Documents/R/HD_Econometrics")

### 0. User defined func
Delta = function(a,b){
  pnorm(a+b)-pnorm(a-b)
}

dPS = function(x,beta=.5,rho=.7,n=100){
  # variances are normalized to 1
  c = sqrt(log(n))
  probaR = Delta(sqrt(n)*beta,c)
  y = probaR*dnorm(x,mean=-rho*sqrt(n)*beta,sd=sqrt(1-rho^2)) +
    (1-Delta((sqrt(n)*beta+rho*x)/sqrt(1-rho^2),c/sqrt(1-rho^2)))*dnorm(x)
  return(list(y=y,
              probaR=probaR))
}

### 1. Simple graph
xset = seq(-3,3,by=.1)
yset = dPS(xset,beta=0,rho=.9)
print(yset$probaR)

plot(xset,yset$y,type="line")

### 1. Graph for handout
## G1
Myset = cbind(dPS(xset,beta=.5)$y,
              dPS(xset,beta=.3)$y,
              dPS(xset,beta=.2)$y,
              dPS(xset,beta=.1)$y)

pdf("PSDensity-rho7.pdf", width=8,height=6)
matplot(xset,Myset,type="l",
        xlab="",ylab="",lwd=c(2,2,2,2))
legend(1.8,.47, c(".5", ".3", ".2",".1"),
       col=seq_len(4), lty=seq_len(4), lwd=c(2,2,2,2))
dev.off()

## G2
Myset = cbind(dPS(xset,beta=.5,rho=.9)$y,
              dPS(xset,beta=.3,rho=.9)$y,
              dPS(xset,beta=.2,rho=.9)$y,
              dPS(xset,beta=.1,rho=.9)$y)

pdf("PSDensity-rho9.pdf", width=8,height=6)
matplot(xset,Myset,type="l",
        xlab="",ylab="",lwd=c(2,2,2,2))
legend(1.8,.77, c(".5", ".3", ".2",".1"),
       col=seq_len(4), lty=seq_len(4), lwd=c(2,2,2,2))
dev.off()

## G2
Myset = cbind(dPS(xset,beta=.5,rho=.4)$y,
              dPS(xset,beta=.3,rho=.4)$y,
              dPS(xset,beta=.2,rho=.4)$y,
              dPS(xset,beta=.1,rho=.4)$y)

pdf("PSDensity-rho4.pdf", width=8,height=6)
matplot(xset,Myset,type="l",
        xlab="",ylab="",lwd=c(2,2,2,2))
legend(1.8,.37, c(".5", ".3", ".2",".1"),
       col=seq_len(4), lty=seq_len(4), lwd=c(2,2,2,2))
dev.off()

## Test @b=0
Myset = cbind(dPS(xset,beta=0,rho=.1,n=100)$y,
              dPS(xset,beta=0,rho=.3,n=100)$y,
              dPS(xset,beta=0,rho=.5,n=100)$y,
              dPS(xset,beta=0,rho=.7,n=100)$y)

matplot(xset,Myset,type="l",
        xlab="",ylab="",lwd=c(2,2,2,2))
legend(1.8,.47, c(".1", ".3", ".5",".7"),
       col=seq_len(4), lty=seq_len(4), lwd=c(2,2,2,2))

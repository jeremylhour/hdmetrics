### Fisher example
### Jeremy L Hour
### 22/11/2017
### MODIFIED: 22/03/2019

setwd("//ulysse/users/JL.HOUR/1A_These/B. ENSAE Classes/Cours3A/")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("ggplot2")
library("gridExtra")

### Use-defined Functions
theta.iter = function(D,Y,Dstar,C=0){
  ChangeTreat = (Dstar==1-D)
  Ystar = Y
  Ystar[ChangeTreat] = Y -(2*Dstar-1)*C
  return(mean(Ystar[Dstar==1]) - mean(Ystar[Dstar==0]))
}

permutation.iter = function(D,Y,C=0){
  return(theta.iter(D,Y,sample(D),C))
}

compute.pval = function(theta.obs, theta.sim,C=0){
  B = length(theta.sim) 
  return((1 + sum(abs(theta.sim-C)>abs(theta.obs-C)))/(B+1))
}


### 1. Simulations

set.seed(12071990)

# param
N = 200
pi = .2
B = 10000

### Case 1
mu = 0.75

d = ifelse(runif(N) < pi,1,0)
y = d*mu + rnorm(N)

theta.hat = mean(y[d==1]) - mean(y[d==0])
print(theta.hat) # Theta


Dpermut = replicate(B, sample(d)) # each column is a random permutation of d

theta.reshuffled = mapply(function(b) theta.iter(d,y,Dpermut[,b],C=0), b=1:B)
p.val1 = compute.pval(theta.hat,theta.reshuffled,C=0)

titer = data.frame(val=theta.reshuffled)
p1 = ggplot(titer, aes(x=val)) + 
  geom_histogram(binwidth = .01, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(name="Estimated ATE") +
  ggtitle(expression(tau*"=.5")) + 
  geom_vline(xintercept = theta.hat, colour="darkorchid3", size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")


C.val = seq(-3,3,.1)
p.vals1 = vector(length=length(C.val))
for(i in 1:length(C.val)){
  theta.sim = mapply(function(b) theta.iter(d,y,Dpermut[,b],C=C.val[i]), b=1:B)
  p.vals1[i] = compute.pval(theta.hat,theta.sim,C=C.val[i])
}

### Case 2
mu = 0

d = ifelse(runif(N) < pi,1,0)
y = d*mu + rnorm(N)

theta.hat = mean(y[d==1]) - mean(y[d==0])
print(theta.hat) # Theta


Dpermut = replicate(B, sample(d)) # each column is a random permutation of d

theta.reshuffled = mapply(function(b) theta.iter(d,y,Dpermut[,b],C=0), b=1:B)
p.val2 = compute.pval(theta.hat,theta.reshuffled,C=0)

titer = data.frame(val=theta.reshuffled)
p2 = ggplot(titer, aes(x=val)) + 
  geom_histogram(binwidth = .01, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(name="Estimated ATE") +
  ggtitle(expression(tau*"=0")) + 
  geom_vline(xintercept = theta.hat, colour="darkorchid3", size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")


C.val = seq(-3,3,.1)
p.vals2 = vector(length=length(C.val))
for(i in 1:length(C.val)){
  theta.sim = mapply(function(b) theta.iter(d,y,Dpermut[,b],C=C.val[i]), b=1:B)
  p.vals2[i] = compute.pval(theta.hat,theta.sim,C=C.val[i])
}

print(p.val1)
print(p.val2)

# Plot 1
pdf("FisherTest.pdf", width=10, height=5)
grid.arrange(p1, p2, ncol=2)
dev.off()

# Plot 2
plotdata = data.frame(val1=p.vals1,val2=p.vals2,C=C.val)

pdf("p-val.pdf", width=8, height=6)
ggplot(plotdata, aes(x=C.val)) + 
  geom_line(aes(y = val1, colour = "red")) +
  geom_line(aes(y = val2, colour = "blue")) +
  geom_point(aes(y = val1, colour = "red")) +
  geom_point(aes(y = val2, colour = "blue")) +
  geom_hline(yintercept = .05, colour="darkorchid3", size=1)+
  scale_x_continuous(name="C", limits=c(-1.5,2.5)) +
  scale_y_continuous(name="p-value") +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
dev.off()

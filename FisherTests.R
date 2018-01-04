### Fisher example
### Jeremy L Hour
### 22/11/2017

setwd("//ulysse/users/JL.HOUR/1A_These/B. ENSAE Classes/Cours3A/")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("ggplot2")
require(gridExtra)

### Permutation function
permutation.iter = function(d,y){
  dstar = sample(d)
  return(mean(y[dstar==1]) - mean(y[dstar==0]))
}

### Compute p.val
compute.pval = function(y,d,dpermut,C,theta.obs){
  Ypot0 = y - d*C; Ypot1 = y + (1-d)*C;
  theta.reshuffled = mapply(function(r) mean(Ypot1[dpermut[,r]==1]) - mean(Ypot0[dpermut[,r]==0]), 1:ncol(dpermut))
  p.val = mean(abs(theta.reshuffled - C) >= abs(theta.obs-C))
  return(p.val)
}


set.seed(12071990)
### Case 1
n = 200
mu = .5
pi = .2
d = ifelse(runif(n) < pi,1,0)
y = d*mu + rnorm(n)

theta.hat = mean(y[d==1]) - mean(y[d==0])
print(theta.hat) ## Original sample stat

B = 10000
theta.reshuffled = replicate(B, permutation.iter(d,y), simplify="vector")
p.val1 = sum(abs(theta.hat) < abs(theta.reshuffled))/B

titer = data.frame(val=theta.reshuffled)

p1 = ggplot(titer, aes(x=val)) + 
  geom_histogram(binwidth = .01, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(name="Estimated ATE") +
  ggtitle(expression(tau*"=.5")) + 
  geom_vline(xintercept = theta.hat, colour="darkorchid3", size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")

dpermut = replicate(B, sample(d)) # each column is a random permutation of d

C.val = seq(-3,3,.1)
p.vals1 = vector(length=length(C.val))
for(i in 1:length(C.val)){
  p.vals1[i] = compute.pval(y,d,dpermut,C.val[i],theta.hat)
}

### Case 2
mu = 0
y = d*mu + rnorm(n)

theta.hat = mean(y[d==1]) - mean(y[d==0])
print(theta.hat) ## Original sample stat

theta.reshuffled = replicate(B, permutation.iter(d,y), simplify="vector")
p.val2 = sum(abs(theta.hat) < abs(theta.reshuffled))/B

titer = data.frame(val=theta.reshuffled)

p2 = ggplot(titer, aes(x=val)) + 
  geom_histogram(binwidth = .01, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(name="Estimated ATE") +
  ggtitle(expression(tau*"=0")) + 
  geom_vline(xintercept = theta.hat, colour="darkorchid3", size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")

dpermut = replicate(B, sample(d)) # each column is a random permutation of d

p.vals2 = vector(length=length(C.val))
for(i in 1:length(C.val)){
  p.vals2[i] = compute.pval(y,d,dpermut,C.val[i],theta.hat)
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
  geom_point(aes(y = val1, colour = "red",size=1)) +
  geom_line(aes(y = val1, colour = "red")) +
  geom_point(aes(y = val2, colour = "blue",size=1)) +
  geom_line(aes(y = val2, colour = "blue")) +
  geom_hline(yintercept = .05, colour="darkorchid3", size=1)+
  scale_x_continuous(name="C", limits=c(-.75,1.25)) +
  scale_y_continuous(name="p-value") +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
dev.off()

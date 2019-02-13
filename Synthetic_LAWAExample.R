### Legal Arizona Workers Act
### From Bohn, Lofstrom and Raphael (2014) 
### "Did the 2007 Legal Arizona Workers Act Reduce the State's unauthorized immigrant population?"
### Jeremy L Hour
### Last Edited: 15/01/2018

### Set working directory
setwd("//ulysse/users/JL.HOUR/1A_These/B. ENSAE Classes/Cours3A/hdmetrics")

rm(list=ls())


### Load packages
library("Synth")

### 0. Load Data
load("data/LAWAdata.Rda")

# Generate Proportion with HS diploma or less
LAWAdata[,"hispnoncitizenHS"] = LAWAdata[,"hispnoncitizene11545"]+LAWAdata[,"hispnoncitizene21545"]

# Drop States with legislation targeting employment of undocumented immigrants
data = LAWAdata[!LAWAdata$statefip%in%c(28,44,45,49),]

data = data[order(data$statefip,data$year),]
data[,"statefip"] = as.numeric(data$statefip)
StateName = c('Alabama','Alaska','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','District of Columbia','Florida',
              'Georgia','Hawaii','Idaho','Illinois','Indiana','Iowa','Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts','Michigan',
              'Minnesota','Missouri','Montana','Nebraska','Nevada','New Hampshire','New Jersey','New Mexico','New York','North Carolina','North Dakota',
              'Ohio','Oklahoma','Oregon','Pennsylvania','South Dakota','Tennessee','Texas','Vermont','Virginia','Washington','West Virginia',
              'Wisconsin','Wyoming')

controls = unique(data[,"statefip"]) 
controls = controls[controls!=4] # Control States are all states except AZ (4)


### 1. Difference-in-differences

# First model
data[,"TreatPeriod"] =  data[,"year"] > 2006
data[,"TreatUnit"] =  data[,"statefip"] == 4
DiD = lm(hispnoncitizenHS ~ TreatUnit*TreatPeriod, data=data)
summary(DiD)

# More complex model
data[,"Y2007"] =  data[,"year"] == 2007
data[,"Y2008"] =  data[,"year"] == 2008
data[,"Y2009"] =  data[,"year"] == 2009
FlexDiD = lm(hispnoncitizenHS ~ Y2007 + Y2008 + Y2009 + TreatUnit + TreatUnit*Y2007 + TreatUnit*Y2008 + TreatUnit*Y2009,  data=data)
summary(FlexDiD)

# Placebo check A: use a ficticious treatment year
ATET = vector(length=6)
for(FakeYear in 2001:2006){
  databis = data[data[,"year"]<FakeYear+1,]
  databis[,"TreatPeriod"] =  databis[,"year"] > FakeYear-1
  PDiD = lm(hispnoncitizenHS ~ TreatUnit*TreatPeriod, data=databis)
  ATET[FakeYear-2000] = coef(PDiD)[4]
}
print(ATET)

# Placebo check B: use a ficticious treated state
ATET = vector(length=length(controls))
databis = data
i=0
for(FakeTreat in controls){
  i=i+1
  databis[,"TreatUnit"] =  data[,"statefip"] == FakeTreat
  PDiD = lm(hispnoncitizenHS ~ TreatUnit*TreatPeriod, data=databis)
  ATET[i] = coef(PDiD)[4]
}
print(ATET)

# Common Trend Assumption?
Ymat = matrix(data[,"hispnoncitizenHS"],nrow=12)
plotdata = ts(cbind(Ymat[,3],apply(Ymat[,-3],1,mean)),start=c(1998), freq=1)

plot(plotdata, plot.type="single",
     col=c("steelblue","darkorchid"), lwd=2,
     lty=c(1,6),xlab="", ylab="Prop of Hispanic non-citizen",
     ylim=c(0,.15))
lim <- par("usr")
rect(2007, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(2002,.08,
       legend=c("Arizona", "Rest of US"),
       col=c("steelblue","darkorchid"), lwd=2,
       lty=c(1,6))


### 2. Synthetic Controls
# Prepare data
dataprep.out = dataprep(foo = data,dependent = "hispnoncitizenHS",unit.variable = "statefip",time.variable = "year",
    special.predictors = list(
      list("ind1",1998:2000,"mean"),list("ind1",2001:2003,"mean"),list("ind1",2004:2006,"mean"),
      list("ind2",1998:2000,"mean"),list("ind2",2001:2003,"mean"),list("ind2",2004:2006,"mean"),
      list("ind3",1998:2000,"mean"),list("ind3",2001:2003,"mean"),list("ind3",2004:2006,"mean"),
      list("ind4",1998:2000,"mean"),list("ind4",2001:2003,"mean"),list("ind4",2004:2006,"mean"),
      list("ind5",1998:2000,"mean"),list("ind5",2001:2003,"mean"),list("ind5",2004:2006,"mean"),
      list("ind6",1998:2000,"mean"),list("ind6",2001:2003,"mean"),list("ind6",2004:2006,"mean"),
      list("ind8",1998:2000,"mean"),list("ind8",2001:2003,"mean"),list("ind8",2004:2006,"mean"),
      list("ind9",1998:2000,"mean"),list("ind9",2001:2003,"mean"),list("ind9",2004:2006,"mean"),
      list("educ1",1998:2000,"mean"),list("educ1",2001:2003,"mean"),list("educ1",2004:2006,"mean"),
      list("educ2",1998:2000,"mean"),list("educ2",2001:2003,"mean"),list("educ2",2004:2006,"mean"),
      list("educ3",1998:2000,"mean"),list("educ3",2001:2003,"mean"),list("educ3",2004:2006,"mean"),
      list("RATE",1998:2000,"mean"),list("RATE",2001:2003,"mean"),list("RATE",2004:2006,"mean")
    ),
    treatment.identifier = 4,
    controls.identifier = controls,
    time.predictors.prior = c(1998:2006),
    time.optimize.ssr = c(1998:2006),
    time.plot = 1998:2009
  )

# Run synthetic Controls: Example from the paper
synthreg = synth(dataprep.out)

# Composition
Wsynth = ifelse(synthreg$solution.w < 1e-5,0,synthreg$solution.w)
Wsynth = Wsynth / sum(Wsynth)
rownames(Wsynth) = StateName[-3]
print(Wsynth)

gaps = dataprep.out$Y1plot- dataprep.out$Y0plot%*%synthreg$solution.w
print(gaps)

synth.tables = synth.tab(dataprep.res = dataprep.out,synth.res = synthreg)
print(synth.tables)

path.plot(synthreg,dataprep.out,tr.intake=2007)
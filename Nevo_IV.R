#############################################################################################################################
#############################################################################################################################
rm(list=ls())

#### import Dataset (semi-fabricated) replicated Aviv Nevo's dataset
path = "C:/Users/gaillac/Documents/GitHub/hdmetrics/data/cereal_ps3.csv"
data_nevo <- read.csv(path, header=T, sep=";",dec=",")
head(data_nevo)
# data_nevo <- as.matrix(data_nevo)
data_nevo  <- as.data.frame(data_nevo)
## Baseline OLS
baseline = colnames(data_nevo)[12:27]
controls =colnames(data_nevo)[10:11]
# controls =colnames(BLP)[c(8:11)]
Result = matrix(ncol=3, nrow=4)
data_nevo$s <- log(data_nevo$share)
endog = "price"
y_n <- "s"


form <- paste(y_n, paste(c(controls,endog), collapse=" + "), sep=" ~ ")
fit.ols.b <- lm(as.formula(form), data=data_nevo)
sds <- coef(summary(fit.ols.b ))[, "Std. Error"]
## own price elasticities 
elas <-  fit.ols.b$coefficients["price"]*data_nevo[,"price"]*(1-as.numeric(data_nevo[,"share" ]))
nb.ine <- sum(elas >= -1)
Result[1,] <- c(fit.ols.b$coefficients["price"], sds["price"],nb.ine)


## Baseline 2SLS
form <- paste(y_n, paste(c(controls,endog), collapse=" + "), sep=" ~ ")
form <- paste(form, paste(c(baseline), collapse=" + "), sep=" | ")
fit.tsls.b <- tsls(as.formula(form), data=data_nevo)
## own price elasticities 
elas <-  fit.tsls.b$coefficients["price",]*data_nevo[,"price"]*(1-data_nevo[,"price" ])
nb.ine <- sum(elas >= -1)
Result[2,] <- c(fit.tsls.b$coefficients["price",], fit.tsls.b $se["price"],nb.ine)

### 2SLS Estimates With \Double Selection"

## Baseline 2SLS Selection
fit.lasso.b <-rlassoIV(x=as.matrix(data_nevo[,controls]),
                       d=as.matrix(data_nevo[,endog]),y=as.matrix(data_nevo[,y_n]),z=as.matrix(data_nevo[,baseline]), select.X=TRUE, select.Z=TRUE)
## own price elasticities 
elas <- fit.lasso.b$coefficients*data_nevo[,"price"]*(1-data_nevo[,"share" ])
nb.ine <- sum(elas >= -1)
Result[3,] <- c(fit.lasso.b$coefficients,fit.lasso.b$se,nb.ine)


## Augmented 2SLS Selection
cont <-  paste("s ~ ",paste(paste0("(",paste(c(controls), collapse=" + ")),")^3"))
zz <-  paste("s ~ ",paste(paste0("(",paste(c(baseline), collapse=" + ")),")^3"))
xsel <- model.matrix(as.formula(cont), data = data_nevo)
zsel <- model.matrix(as.formula(zz), data = data_nevo)
dim(x)
fit.lasso.aug <-rlassoIV(x=xsel,
                         d=as.matrix(data_nevo[,endog]),y=as.matrix(data_nevo[,y_n]),
                         z=zsel, select.X=TRUE, select.Z=TRUE)
elas <- fit.lasso.aug$coefficients*data_nevo[,"price"]*(1-data_nevo[,"share" ])
nb.ine <- sum(elas >= -1)
Result[4,] <- c(fit.lasso.aug$coefficients,fit.lasso.aug$se,nb.ine)

## own price elasticities 
Result <- as.data.frame(Result)
row.names(Result) = c("Baseline OLS","Baseline 2SLS","Baseline 2SLS Selection","Augmented 2SLS Selection")
colnames(Result) <- c("Price Coefficient","Standard Error","Number Inelastic")
print(Result)
library(xtable)
xtable(Result)












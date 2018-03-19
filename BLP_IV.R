# if (!require(hdm)) install.packages("hdm")
library(hdm)

## load dataset, of import it from /data/ in the Github repository
data(BLP)
BLP <- as.data.frame(BLP)
head(BLP)

## Note: despite what the hdm package's documentation says, price is price - mean(price) (in /1000 1983$).
## If you compare the quantiles of price with Table II from BLP (1995), you see that everything is
## shifted down by 11.761
BLP$BLP.price <- BLP$BLP.price + 11.761
## (check with table II from BLP 1995)
quantile(BLP$BLP.price)
quantile(BLP$BLP.hpwt)
quantile(BLP$BLP.mpd)
quantile(BLP$BLP.y)
quantile(BLP$BLP.mpg)
BLP$BLP.mpd <- BLP$BLP.mpd*10
BLP$BLP.mpg <- BLP$BLP.mpg*10
Result = matrix(ncol=3, nrow=6)
###  Estimates Without Selection

## Baseline OLS IV
baseline = colnames(BLP)[64:73]
## Augmented  OLS IV
augmented = colnames(BLP)[16:63]
controls =colnames(BLP)[c(c(7:8),c(10:11))]
# controls =colnames(BLP)[c(8:11)]
controls_cont = controls[-1]

## log shares
BLP$s <- log(BLP$BLP.share)
endog = "BLP.price"
y_n <- "s"
head(BLP)

form <- paste(y_n, paste(c(controls,endog), collapse=" + "), sep=" ~ ")
fit.ols.b <- lm(as.formula(form), data=BLP)
sds <- coef(summary(fit.ols.b ))[, "Std. Error"]
## own price elasticities 
elas <-  fit.ols.b$coefficients["BLP.price"]*BLP[,"BLP.price"]*(1-BLP[,"BLP.share" ])
nb.ine <- sum(elas >= -1)
Result[1,] <- c(fit.ols.b$coefficients["BLP.price"], sds["BLP.price"],nb.ine)


## Baseline 2SLS
form <- paste(y_n, paste(c(controls,endog), collapse=" + "), sep=" ~ ")
form <- paste(form, paste(c(baseline), collapse=" + "), sep=" | ")
fit.tsls.b <- tsls(as.formula(form), data=BLP)
## own price elasticities 
elas <-  fit.tsls.b$coefficients["BLP.price",]*BLP[,"BLP.price"]*(1-BLP[,"BLP.share" ])
nb.ine <- sum(elas >= -1)
Result[2,] <- c(fit.tsls.b$coefficients["BLP.price",], fit.tsls.b $se["BLP.price"],nb.ine)


## Augmented OLS
cont <-  paste(paste("(",paste(c(controls), collapse=" + ")),")^3")
form <- paste(y_n,paste(cont,"+ BLP.trend + BLP.air + BLP.price"), sep=" ~ ")
fit.ols.b <- lm(as.formula(form), data=BLP)
sds <- coef(summary(fit.ols.b ))[, "Std. Error"]
## own price elasticities 
elas <-  fit.ols.b$coefficients["BLP.price"]*BLP[,"BLP.price"]*(1-BLP[,"BLP.share" ])
nb.ine <- sum(elas >= -1)
Result[3,] <- c(fit.ols.b$coefficients["BLP.price"], sds["BLP.price"],nb.ine)

## Augmented TSLS
cont <-  paste(paste("(",paste(c(controls), collapse=" + ")),")^3")
form <- paste(y_n,paste(cont,"+ BLP.trend + BLP.air + BLP.price"), sep=" ~ ")
form <- paste(form, paste(c(augmented), collapse=" + "), sep=" | ")
fit.tsls.b  <-  tsls(as.formula(form), data=BLP)
## own price elasticities 
elas <-  fit.tsls.b$coefficients["BLP.price",]*BLP[,"BLP.price"]*(1-BLP[,"BLP.share" ])
nb.ine <- sum(elas >= -1)
Result[4,] <- c(fit.tsls.b$coefficients["BLP.price",], fit.tsls.b $se["BLP.price"],nb.ine)

### 2SLS Estimates With \Double Selection"

## Baseline 2SLS Selection
fit.lasso.b <-rlassoIV(x=as.matrix(BLP[,controls]),
                         d=as.matrix(BLP[,endog]),y=as.matrix(BLP[,y_n]),z=as.matrix(BLP[,baseline]), select.X=TRUE, select.Z=TRUE)
## own price elasticities 
elas <- fit.lasso.b$coefficients*BLP[,"BLP.price"]*(1-BLP[,"BLP.share" ])
nb.ine <- sum(elas >= -1)
Result[5,] <- c(fit.lasso.b$coefficients,fit.lasso.b$se,nb.ine)


## Augmented 2SLS Selection
cont <-  paste("s ~ ",paste(paste0("(",paste(c(controls), collapse=" + ")),")^3"))
xsel <- model.matrix(as.formula(cont), data = BLP)
x <- cbind(xsel[,-1],BLP[,c("BLP.air")])
dim(x)
fit.lasso.aug <-rlassoIV(x=x,
                          d=as.matrix(BLP[,endog]),y=as.matrix(BLP[,y_n]),z=as.matrix(BLP[,augmented]), select.X=TRUE, select.Z=TRUE)
## own price elasticities 

elas <- fit.lasso.aug$coefficients*BLP[,"BLP.price"]*(1-BLP[,"BLP.share" ])
nb.ine <- sum(elas >= -1)
Result[6,] <- c(fit.lasso.aug$coefficients,fit.lasso.aug$se,nb.ine)

Result <- as.data.frame(Result)
row.names(Result) = c("Baseline OLS","Baseline 2SLS","Augmented OLS"
                           , "Augmented 2SLS","Baseline 2SLS Selection","Augmented 2SLS Selection")
colnames(Result) <- c("Price Coefficient","Standard Error","Number Inelastic")
print(Result)
library(xtable)
xtable(Result)











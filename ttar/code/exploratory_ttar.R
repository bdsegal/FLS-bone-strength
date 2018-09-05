# Exploratory plots for distribution of total area (ttar)

library(mgcv)

# set computer-specific paths (note: 'path' points to the folder
# containing the data, and is used in 'data_pred.R')
# alter as necessary
if( length(grep("bdsegal",getwd()))>0 ){
  computer <- "C:/Users/bdsegal"
  path <- "Q:/U/bdsegal/Bones/data"
} else if (length(grep("bsegal", getwd()))>0){
  computer <- "/home/bsegal"
  path <- "/home/bsegal/Documents/Research/Bones/data"
} else{
  computer <- "C:/Users/Segal"
  path <- "C:/Users/Segal/Documents/Research/Bones/data"
}

setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))


# functions -------------------------------------------------------------------
gammaOverlay <- function(x, variable_name, b=25){
  mu <- mean(x)
  sigma2 <- var(x)

  alpha <- mu^2/sigma2
  beta <- sigma2/mu

  hist(x,freq=FALSE,breaks=b, main = paste("Histogram of ", variable_name, sep=""), xlab="", ylab="density")
  curve(dgamma(x,shape=alpha,scale=beta),xlim=range(x),add=TRUE,col="red")
  
  readline()
  
 qqplot(qgamma(ppoints(500),shape=alpha,scale=beta),x,pch=19,
    main=paste("Gamma Q-Q plot of ", variable_name, sep=""),
    xlab="theoretical", ylab="observed")
  qqline(x, distribution = function(p) qgamma(p, shape=alpha,scale=beta), col = "red")
}

normalOverlay <- function(x, variable_name, b = 25){
  mu <- mean(x)
  sigma <- sd(x)

  hist(x, freq = FALSE, breaks = b, main = paste("Histogram of ", variable_name, sep=""), xlab="", ylab="density")
  curve(dnorm(x, mean = mu, sd = sigma), xlim = range(x), add = TRUE, col = "red")
  
  readline()
  qqplot(qnorm(ppoints(500), mean = mu, sd = sigma), x, pch = 19,
    main = paste("Normal Q-Q plot of ", variable_name, sep=""), xlab = "theoretical", ylab = "observed")
  qqline(x, distribution = function(p) qnorm(p, mean = mu, sd = sigma), col = "red")
}

gammaOverlay(dataSub$ttar[which(dataSub$targage==8 & dataSub$sex == "Male")], "TTAR")
gammaOverlay(dataSub$ttar[which(dataSub$targage==18 & dataSub$sex == "Male")], "TTAR")

gammaOverlay(dataSub$ttar[which(dataSub$targage==8 & dataSub$sex == "Female")], "TTAR")
gammaOverlay(dataSub$ttar[which(dataSub$targage==18 & dataSub$sex == "Female")], "TTAR")

m1ttarCent <- gamm(ttar ~ 
        te(skelage, birthday, k=c(5,5), bs="cr")+
        te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
        family=Gamma(link=log),
        random=list(pedno=~1, ptno=~1),
        data=dataSub)

# or, after running ttar_modelFits, load m1ttarcent
# load(file.path("ttar","code","m1ttarCent.Rdata"))

dataSub$expectedVals <- exp(predict(m1ttarCent$lme))
dataSub$errors <- with(dataSub, ttar - expectedVals)

m1 <- lm(log(dataSub$errors^2) ~ log(dataSub$expectedVals))
coef(m1)
              # (Intercept) log(dataSub$expectedVals) 
              #   -7.270134                  1.944297 

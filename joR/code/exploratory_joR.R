# Exploratory plots for distribution of jo/R (BSI)

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

gammaOverlay(data$joR[which(data$targage==8 & data$sex == "Male")], "BSI")
gammaOverlay(data$joR[which(data$targage==18 & data$sex == "Male")], "BSI")

gammaOverlay(data$joR[which(data$targage==8 & data$sex == "Female")], "BSI")
gammaOverlay(data$joR[which(data$targage==18 & data$sex == "Female")], "BSI")

m1cent <- gamm(joR ~ 
        te(skelage, birthday, k=c(5,5), bs="cr")+
        te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
        family=Gamma(link=log),
        random=list(pedno=~1, ptno=~1),
        data=dataSub)

# or, after running joR_modelFits, load m1cent
# load(file.path("joR","code","m1cent.Rdata"))

dataSub$expectedVals <- exp(predict(m1cent$lme))
dataSub$errors <- with(dataSub, joR - expectedVals)

residFit <- lm(log(dataSub$errors^2) ~ log(dataSub$expectedVals))
coef(residFit)

normalOverlay(log(data$joR[which(data$targage==8 & data$sex == "Male")]), "log(BSI)")
normalOverlay(log(data$joR[which(data$targage==18 & data$sex == "Male")]), "log(BSI)")

normalOverlay(log(data$joR[which(data$targage==8 & data$sex == "Female")]), "log(BSI)")
normalOverlay(log(data$joR[which(data$targage==18 & data$sex == "Female")]), "log(BSI)")

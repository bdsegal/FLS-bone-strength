# exploring models for body size
# This code was then incorporated into the dataPred.R file
# This script is not used directly in any subsequent data prep

library(mgcv)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cubature)

# set computer-specific paths (note: 'path' points to the folder
# containing the data, and is used in 'dataPred.R')
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

dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

paperPath <- file.path(computer,
  "Dropbox/Research/Bones/final_analysis/plots_report/all")

# load functions for obtaining predictions
source(file.path(dataPrepPath,"predict_functions_2prod_centered.R"))

overlay <- function(outcome,b=25){
# make overlay plots for gamma distribution
  mu <- mean(outcome)
  sigma2 <- var(outcome)

  alpha <- mu^2/sigma2
  beta <- sigma2/mu

  hist(outcome,freq=FALSE,breaks=b)
  curve(dgamma(x,shape=alpha,scale=beta),xlim=range(outcome),add=TRUE,col="red")

  readline()
  
  qqplot(qgamma(ppoints(500),shape=alpha,scale=beta),outcome,pch=19,
    main="Q-Q plot for Gamma",
    xlab="Theoretical",ylab=expression(j[0]))
  qqline(outcome, distribution = function(p) qgamma(p, shape=alpha,scale=beta), col = "red")
}

data <- read.csv(file.path(path,"skeleton_long_INCLUDE.csv"), na.strings=".")

data$sex[which(data$sex==1)] <- "Male"
data$sex[which(data$sex==2)] <- "Female"
data$male <- ifelse(data$sex=="Male",1,0)

# data$jo_over_radius <- with(data, jo/(mc2tmd/2))
data$joR <- with(data, jo/(mc2tmd/2))

data$ptno <- as.factor(data$ptno)
data$pedno <- as.factor(data$pedno)

# fix missing skelage ------------------------------------------------

skelageToFix <- which(with(data, is.na(skelage) & age>=18))
data[skelageToFix,"skelage"] <- 18

# Body size standardization ------------------------------------------

# complete cases
keep <- which(complete.cases(data[,c("skelage","birthday","bodysize","male", "ptno","pedno")]))
dataSub <- data[keep,]

# keep only those above target age 8
dataSub <- dataSub[which(dataSub$targage>=8),]
# sex-specific standardization ----------------------------

# get percent male to make population average predictions
sexTable <- table(dataSub[,"sex"])
percMale <- sexTable[2]/sum(sexTable)

# data for predictions
dataPred <- as.data.frame(cbind(skelage = dataSub$skelage, male=percMale))

b1 <- gamm(bodysize ~ s(skelage, k=25)+
  s(skelage, k=25,  by=male), data=dataSub, family=Gamma(link=log),
  random=list(pedno=~1, ptno=~1))

png(file.path(paperPath, "b1_resid.png"))
plot(b1$lme)
dev.off()

gam.check(b1$gam)

summary(b1$gam)

# mean-variance relationship
dataSub$expectedVals <- exp(predict(b1$lme))
dataSub$errors <- with(dataSub, bodysize - expectedVals)

residFit <- lm(log(errors^2) ~ log(expectedVals), data = dataSub)
coef(residFit)

      # (Intercept) log(expectedVals) 
       # -10.049502          1.737098 

# 95% confidence interval
coef(residFit)[2] + c(-1, 1) * 1.96 * sqrt(vcov(residFit)[2, 2])
# [1] 1.008372 2.465828

overlay(data$bodysize[which(data$targage==8)],b=40)
overlay(data$bodysize[which(data$targage==10)],b=40)
overlay(data$bodysize[which(data$targage==12)],b=40)
overlay(data$bodysize[which(data$targage==14)],b=40)
overlay(data$bodysize[which(data$targage==16)],b=40)
overlay(data$bodysize[which(data$targage==18)],b=40)

skelSeq <- seq(8, 18, .1)
linPredMale <- predict(b1$gam, newdata=data.frame(skelage=skelSeq, male=1))
linPredFemale <- predict(b1$gam, newdata=data.frame(skelage=skelSeq, male=0))
linPredOverall <- predict(b1$gam, newdata=data.frame(skelage=skelSeq, male=percMale))

vc <- VarCorr(b1$lme)

sigmaPed <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtno <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigma <- diag(c(sigmaPed, sigmaPtno)^2)
logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))

lim <- c(sigmaPed, sigmaPtno) * 100

meanMale <- sapply(linPredMale, function(x){
                     hcubature(integrate2dRandomEffects,
                               lowerLimit = -lim,
                               upperLimit = lim, 
                               linPred = x,
                               sigma = sigma,
                               mean = c(0, 0),
                               logdet = logdet,
                               tol=1e-4,
                               vectorInterface = TRUE)$integral})

meanFemale <- sapply(linPredFemale, function(x){
                     hcubature(integrate2dRandomEffects,
                               lowerLimit = -lim,
                               upperLimit = lim, 
                               linPred = x,
                               sigma = sigma,
                               mean = c(0, 0),
                               logdet = logdet,
                               tol=1e-4,
                               vectorInterface = TRUE)$integral})

meanOverall <- sapply(linPredOverall, function(x){
                     hcubature(integrate2dRandomEffects,
                               lowerLimit = -lim,
                               upperLimit = lim, 
                               linPred = x,
                               sigma = sigma,
                               mean = c(0, 0),
                               logdet = logdet,
                               tol=1e-4,
                               vectorInterface = TRUE)$integral})

bodysizePlot <- data.frame(skelage=rep(skelSeq,3),
                           Sex=rep(c("Male","Female","Overall"), each=length(skelSeq)),
                           bodysize=c(meanMale, meanFemale, meanOverall))

dev.new(height = 5, width = 8)
ggplot(aes(x=skelage, y=bodysize, color=Sex, linetype=Sex),
  data=bodysizePlot)+
  # [which(bodysizePlot$Sex %in% c("Male","Female")),])+
  geom_line(size=1)+
  # geom_line(color="black", size=1, linetype="solid", data=bodysizePlot[which(bodysizePlot$Sex=="avg"),])+
  theme_bw(22)+
  labs(x="Skeletal age", y="Body size")+
  scale_color_manual("",values=c("red","blue","black"))+
  scale_linetype_manual("",values=c("dashed","dotdash","solid")) +
  scale_x_continuous(breaks = seq(6, 18, 2))
ggsave(file.path(paperPath, "bodysizeAvg.png"))

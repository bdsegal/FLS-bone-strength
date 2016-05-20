# exploring models for body size
# This code was then incorporated into the dataPred.R file
# This script is not used directly in any subsequent data prep

library(mgcv)

library(dplyr)
library(ggplot2)
library(reshape2)

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

paperPath <- file.path(computer,
  "Dropbox/Research/Bones/final_analysis/joR/paper")

overlay <- function(jo,b=25){
# make overlay plots for gamma distribution
  mu <- mean(jo)
  sigma2 <- var(jo)

  alpha <- mu^2/sigma2
  beta <- sigma2/mu

  hist(jo,freq=FALSE,breaks=b)
  curve(dgamma(x,shape=alpha,scale=beta),xlim=range(jo),add=TRUE,col="red")

  readline()
  
  qqplot(qgamma(ppoints(500),shape=alpha,scale=beta),jo,pch=19,
    main="Q-Q plot for Gamma",
    xlab="Theoretical",ylab=expression(j[0]))
  qqline(jo, distribution = function(p) qgamma(p, shape=alpha,scale=beta), col = "red")
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
# save(b1, file="b1.Rdata")

summary(b1$gam)
errors <- residuals(b1$lme)
expectedVals <- predict(b1$lme)

coef(lm(log(errors^2) ~ log(expectedVals)))
      # (Intercept) log(expectedVals) 
       # -10.049502          1.737098 

summary(lm(log(errors^2) ~ log(expectedVals)))
# 95% confidence interval
1.7371 + c(-1,1)*1.96*0.3718
# [1] 1.008372 2.465828

overlay(data$bodysize[which(data$targage==8)],b=40)
overlay(data$bodysize[which(data$targage==10)],b=40)
overlay(data$bodysize[which(data$targage==12)],b=40)
overlay(data$bodysize[which(data$targage==14)],b=40)
overlay(data$bodysize[which(data$targage==16)],b=40)
overlay(data$bodysize[which(data$targage==18)],b=40)

vc <- VarCorr(b1$lme)
nvc <- nrow(vc)
vc[(nvc-4):nvc,]

predMale <- exp(predict(b1$gam, newdata=data.frame(skelage=seq(8,18,.1), male=1)))
predFemale <- exp(predict(b1$gam, newdata=data.frame(skelage=seq(8,18,.1), male=0)))
predAvg <- exp(predict(b1$gam, newdata=data.frame(skelage=seq(8,18,.1), male=percMale)))

skelSeq <- seq(8,18,.1)
bodysizePlot <- data.frame(
  skelage=rep(skelSeq,3),
  Sex=rep(c("Male","Female","Mean"), each=length(skelSeq)),
  bodysize=c(predMale, predFemale, predAvg)
)

ggplot(aes(x=skelage, y=bodysize, color=Sex, linetype=Sex),
  data=bodysizePlot)+
  # [which(bodysizePlot$Sex %in% c("Male","Female")),])+
  geom_line(size=1)+
  # geom_line(color="black", size=1, linetype="solid", data=bodysizePlot[which(bodysizePlot$Sex=="avg"),])+
  theme_bw(22)+
  labs(x="Skeletal age", y="Body size")+
  scale_color_manual("",values=c("red","blue","black"))+
  scale_linetype_manual("",values=c("dashed","dotdash","solid"))

ggsave(file.path(paperPath, "bodysizeAvg.png"))
# regress body size on skelage, and get expected values and residuals


# dataSub$bodysizeResid <- with(dataSub, bodysize - bodysizePred)
plot(residuals(b1$gam, type="response"))
plot(dataSub$bodysizeCentered)

# # regress squared residuals on skelage, and get expected variance
b2 <- gamm(bodysizeCentered^2 ~ s(skelage, k=25), 
    random=list(pedno=~1, ptno=~1),
    data=dataSub, family=Gamma(link=log))
save(b2, file="b2.Rdata")

dataSub$bodysizeVar <- predict(b2$gam,type="response")

dataSub$bodysizeStandard <- with(dataSub, bodysizeCentered/sqrt(bodysizeVar))

plot(x=dataSub$skelage, y=sqrt(dataSub$bodysizeVar))
plot(x=dataSub$skelage, y=dataSub$bodysizeStandard)

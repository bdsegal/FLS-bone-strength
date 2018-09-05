# Functions for obtaining predicted values for models that use
# 2-way tensor product between skelage and birthday
# with centered body size and a GLM with log link.
# Used to fill matrices in multiple imputations ('impute_phvage_multiple_post.R')

# (might already be loaded)
library(reshape2)

# load fits for predicting sex-specific standardized body size
# (might already be loaded)
load(file.path(dataPrepPath,"b1.Rdata"))

# phvageMeanSex <- tapply(dataSub$phvage, dataSub$sex, mean, na.rm=TRUE)
# phvageCentSex <- phvageMeanSex-mean(dataSub$phvage, na.rm=TRUE)

## multivariaten normal random sampling, from Wood
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig); m <- ncol(L);
  (mu + L%*%matrix(rnorm(m*n),m,n)) 
}

plotTrend1 <- function(fit, skelagePred=seq(8,18, 2), 
  birthdayPred=seq(1930,2000, 0.5), B=1000, cti=FALSE){

  # produces regression surface without random effects
  # with Bayesian credible intervals (conditional on smoothing parameter)
  
  nSkelage <- length(skelagePred)
  nBirthday <- length(birthdayPred)
  
  beta <- coef(fit)
  
  basePlusMaleInd <- grep("skelage",names(beta))
  maleInd <- grep("male",names(beta))
  femaleInd <- basePlusMaleInd[-which(basePlusMaleInd %in% maleInd)]

  keepMale <- c(1,basePlusMaleInd)
  keepFemale <- c(1,femaleInd)

  # male specific, with overall mean for male (for predicting bodysize)
  newd0 <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=percMale
            )
      
  # for getting mean for males, with male=1 (for predicting body size)
  newd0Male <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=1,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Male")])
            )

  # diff between avg male and avg overall body size: 
  # bodysize|male=1 - bodysize|male=overall mean
  newd0Male$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd0Male)) - 
         exp(predict(b1$gam, newdata=newd0))
  
  # same as above for females
  newd0Female <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  newd0Female$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd0Female)) - 
         exp(predict(b1$gam, newdata=newd0))

  # get design matrices for males and females
  X0MaleTemp <- predict(fit, newdata=newd0Male, type= "lpmatrix")
  X0FemaleTemp <- predict(fit, newdata=newd0Female, type= "lpmatrix")
  
  # only keep columns corresponding to males and females, respectively
  X0Male <- X0MaleTemp*0
  X0Male[,keepMale] <- X0MaleTemp[,keepMale]
  
  X0Female <- X0FemaleTemp*0
  X0Female[,keepFemale] <- X0FemaleTemp[,keepFemale]
  
  # simulate from posterior of betas
  br <- rmvn(B, beta, fit$Vp)
  
  # get simulated values of joR
  if(!cti){
    simMale <- exp(X0Male %*% br)
    simFemale <- exp(X0Female %*% br)
  }else{
    simMale <- X0Male %*% br
    simFemale <- X0Female %*% br
  }

  return(list(simMale=simMale, simFemale=simFemale))
}

plotDerivative1 <- function(fit, skelagePred=seq(8,18,2), 
  birthdayPred=seq(1930,2000,0.5), B=10000, eps=1e-7, cti=FALSE){

  # produces partial derivative of regression surface wrt birthday
  # via finite differences without random effects
  # with Bayesian credible intervals (conditional on smoothing parameter)
    
  nSkelage <- length(skelagePred)
  nBirthday <- length(birthdayPred)
  
  beta <- coef(fit) 
  
  basePlusMaleInd <- grep("skelage",names(beta))
  maleInd <- grep("male",names(beta))
  femaleInd <- basePlusMaleInd[-which(basePlusMaleInd %in% maleInd)]

  keepMale <- c(1,basePlusMaleInd)
  keepFemale <- c(1,femaleInd)

  # joR at initial point
  newd0 <- data.frame( skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=percMale
            )
            
  newd0Male <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=1,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Male")])
            )

  newd0Male$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd0Male)) - 
         exp(predict(b1$gam, newdata=newd0))

  newd0Female <- data.frame( skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  newd0Female$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd0Female)) - 
         exp(predict(b1$gam, newdata=newd0))
     
  X0MaleTemp <- predict(fit, newdata=newd0Male, type= "lpmatrix")
  X0FemaleTemp <- predict(fit, newdata=newd0Female, type= "lpmatrix")
  
  X0Male <- X0MaleTemp*0
  X0Male[,keepMale] <- X0MaleTemp[,keepMale]
  
  X0Female <- X0FemaleTemp*0
  X0Female[,keepFemale] <- X0FemaleTemp[,keepFemale]
  
  # joR at initial point + eps increment in birthday
  newd1 <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred+eps, each=nSkelage),
            male=percMale
            )
            
  newd1Male <- data.frame( skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred+eps, each=nSkelage),
            male=1,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Male")])
            )

  newd1Male$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd1Male)) - 
         exp(predict(b1$gam, newdata=newd1))

  newd1Female <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred+eps, each=nSkelage),
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  newd1Female$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd1Female)) - 
         exp(predict(b1$gam, newdata=newd1))
  
  XdMaleTemp <- predict(fit, newdata=newd1Male, type= "lpmatrix")
  XdFemaleTemp <- predict(fit, newdata=newd1Female, type= "lpmatrix")
  
  XdMale <- XdMaleTemp*0
  XdMale[,keepMale] <- XdMaleTemp[,keepMale]
  
  XdFemale <- XdFemaleTemp*0
  XdFemale[,keepFemale] <- XdFemaleTemp[,keepFemale]

  # simulate from posterior for beta
  br <- rmvn(B, beta, fit$Vp)
  
  # use finite differences to approximate partial derivative wrt birth year
  if(!cti){
    simMale <- (exp(XdMale %*% br)-exp(X0Male %*% br))/eps
    simFemale <- (exp(XdFemale %*% br)-exp(X0Female %*% br))/eps
  }else{
    simMale <- (XdMale %*% br - X0Male %*% br)/eps
    simFemale <- (XdFemale %*% br - X0Female %*% br)/eps
  }

  return(list(simMale=simMale, simFemale=simFemale,
    skelagePred=skelagePred, birthdayPred=birthdayPred,
    nSkelage=nSkelage,
    B=B))
}

plotPerChange1 <- function(fit, initialBirthday=1930, endBirthday=1995, 
  skelagePred=seq(8,18,.1), B=10000, cti=FALSE){

  # produces percent change

  nSkelage <- length(skelagePred)
  
  beta <- coef(fit) 
  
  basePlusMaleInd <- grep("skelage",names(beta))
  maleInd <- grep("male",names(beta))
  femaleInd <- basePlusMaleInd[-which(basePlusMaleInd %in% maleInd)]

  keepMale <- c(1,basePlusMaleInd)
  keepFemale <- c(1,femaleInd)

  newd0 <- data.frame(skelage=skelagePred,
            birthday=initialBirthday,
            male=percMale
            )
            
  newd0Male <- data.frame(skelage=skelagePred,
            birthday=initialBirthday,
            male=1,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Male")])
            )

  newd0Male$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd0Male)) - 
         exp(predict(b1$gam, newdata=newd0))

  newd0Female <- data.frame(skelage=skelagePred,
            birthday=initialBirthday,
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  newd0Female$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd0Female)) - 
         exp(predict(b1$gam, newdata=newd0))
     
  X0MaleTemp <- predict(fit, newdata=newd0Male, type= "lpmatrix")
  X0FemaleTemp <- predict(fit, newdata=newd0Female, type= "lpmatrix")
  
  X0Male <- X0MaleTemp*0
  X0Male[,keepMale] <- X0MaleTemp[,keepMale]

  X0Female <- X0FemaleTemp*0  
  X0Female[,keepFemale] <- X0FemaleTemp[,keepFemale]
    
  newd1 <- data.frame(skelage=skelagePred,
            birthday=endBirthday,
            male=percMale
            )
            
  newd1Male <- data.frame(skelage=skelagePred,
            birthday=endBirthday,
            male=1,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Male")])
            )

  newd1Male$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd1Male)) - 
         exp(predict(b1$gam, newdata=newd1))

  newd1Female <- data.frame(skelage=skelagePred,
            birthday=endBirthday,
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  newd1Female$bodysizeCentered <- 
         exp(predict(b1$gam, newdata=newd1Female)) - 
         exp(predict(b1$gam, newdata=newd1))
  
  XdMaleTemp <- predict(fit, newdata=newd1Male, type= "lpmatrix")
  XdFemaleTemp <- predict(fit, newdata=newd1Female, type= "lpmatrix")
  
  XdMale <- XdMaleTemp*0  
  XdMale[,keepMale] <- XdMaleTemp[,keepMale]
  
  XdFemale <- XdFemaleTemp*0  
  XdFemale[,keepFemale] <- XdFemaleTemp[,keepFemale]
 
  # simulate from posterior for beta
  br <- rmvn(B, beta, fit$Vp)
  
  # Calculate percent change
  if (!cti){
    simMale <- (exp((XdMale - X0Male) %*% br)-1)
    simFemale <- (exp((XdFemale - X0Female) %*% br)-1)
  }else{
    simMale <- (XdMale %*% br - X0Male %*% br)/ (X0Male %*% br)
    simFemale <- (XdFemale %*% br - X0Female %*% br) / (X0Female %*% br)
  }

  return(list(simMale=simMale, simFemale=simFemale,
    skelagePred=skelagePred, initialBirthday=initialBirthday,
    endBirthday=endBirthday,
    B=B))
}

sliceFun <- function(out, i, skelage = seq(8,18,2)){
  # slices the 'out' objects created by plotTrend and plotDerivative to
  # produce data for 2D plots for each skelage

  m <- out$z[[i]]$mean
  l <- out$z[[i]]$l
  u <- out$z[[i]]$u

  meanM <- melt(m, value.name="mean")
  lM <- melt(l, value.name="l")
  uM <- melt(u, value.name="u")
  
  meanMsub <- meanM[which(meanM$skelage %in% skelage),]
  lMsub <- lM[which(lM$skelage %in% skelage),]
  uMsub <- uM[which(uM$skelage %in% skelage),]
  
  meanMsub$skelage <- factor(meanMsub$skelage)
  levels(meanMsub$skelage) <- paste("Skeletal age", levels(meanMsub$skelage))
  
  lMsub$skelage <- factor(lMsub$skelage)
  levels(lMsub$skelage) <- paste("Skeletal age", levels(lMsub$skelage))
  
  uMsub$skelage <- factor(uMsub$skelage)
  levels(uMsub$skelage) <- paste("Skeletal age", levels(uMsub$skelage))
  
  slice <- merge(merge(meanMsub, lMsub), uMsub)
  
  return(slice)
}

# Functions for obtaining predicted values for models that use
# 2-way tensor product between skelage and birthday
# with centered body size and a GLM with log link
# Used by 'joR_modelEval_cent.R' and 'ttar_modelEval_cent.R')

# (might already be loaded)
library(reshape2)
library(cubature)

# load fits for predicting sex-specific standardized body size
# (might already be loaded)
load(file=file.path(computer,"Dropbox/Research/Bones/final_analysis", "b1.Rdata"))

# integrate function over two-dimensional random effects
integrate2dRandomEffects <- function (b, linPred, mean, sigma, logdet) {
    distval <- stats::mahalanobis(t(b), center = mean, cov = sigma)
    return(exp(matrix(linPred + colSums(b) -(2 * log(2 * pi) + logdet + distval)/2, ncol = ncol(b))))
}

# integrate derivative over two-dimensional random effects
integrate2dRandomEffectsDer <- function(b, linPred0, linPred1, eps, mean, sigma, logdet) {
    distval <- stats::mahalanobis(t(b), center = mean, cov = sigma)
    return(matrix(((exp(linPred1 + colSums(b)) - exp(linPred0 + colSums(b)))/eps) * 
                    exp(-(2 * log(2 * pi) + logdet + distval)/2), ncol = ncol(b)))
}

## multivariaten normal random sampling, from Wood
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  return(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

plotTrend <- function(fit0,
                      sigmaOutcome = NULL,
                      logdetOutcome = NULL,
                      limOutcome = NULL,
                      skelagePred=seq(8,18, 2),
                      birthdayPred=seq(1930,2000, 0.5),
                      alpha=0.05,
                      B=10000,
                      Bmc = 1000,
                      integrateRandomEffects = "none"){

  if(!integrateRandomEffects %in% c("none", "MC", "cubature")) {
    stop("integrateRandomEffects must be equal to 'none', 'MC', or 'cubature'")
  }

  # produces regression surface with Bayesian credible intervals (conditional on smoothing parameter)
  
  fit <- fit0$gam
  
  nSkelage <- length(skelagePred)
  nBirthday <- length(birthdayPred)
  
  beta <- coef(fit)
  
  basePlusMaleInd <- grep("skelage",names(beta))
  maleInd <- grep("male",names(beta))
  femaleInd <- basePlusMaleInd[-which(basePlusMaleInd %in% maleInd)]

  keepMale <- c(1,basePlusMaleInd)
  keepFemale <- c(1,femaleInd)

  newd0 <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=percMale
            )
            
  newd0Male <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=1,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Male")])
            )

  newd0Female <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  ## center body size
  vc <- VarCorr(b1$lme)

  sigmaPed <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
  sigmaPtno <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

  sigma <- diag(c(sigmaPed, sigmaPtno)^2)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))

  lim <- c(sigmaPed, sigmaPtno) * 100

  linPredOverall0 <- predict(b1$gam, newdata=newd0)
  linPredMale0 <- predict(b1$gam, newdata=newd0Male)
  linPredFemale0 <- predict(b1$gam, newdata=newd0Female)

  meanOverall0 <- sapply(linPredOverall0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanMale0 <- sapply(linPredMale0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanFemale0 <- sapply(linPredFemale0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  newd0Male$bodysizeCentered <- meanMale0 - meanOverall0
  newd0Female$bodysizeCentered <- meanFemale0 - meanOverall0

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
  if (integrateRandomEffects == "none") {
    simMale <- exp(X0Male %*% br)
    simFemale <- exp(X0Female %*% br)

  } else if (integrateRandomEffects == "MC") {

    simMaleLin <- X0Male %*% br
    simFemaleLin <- X0Female %*% br

    # MC approach  
    bSum <- colSums(rmvn(n = Bmc, mu = c(0, 0), sig = sigmaOutcome))

    simMale <- array(0, dim = dim(simMaleLin))
    simFemale <- array(0, dim = dim(simFemaleLin))
    for (iter in 1:Bmc) {
      simMale <- simMale + exp(simMaleLin + bSum[iter]) / Bmc
      simFemale <- simFemale + exp(simFemaleLin + bSum[iter]) / Bmc
    }

  } else if (integrateRandomEffects == "cubature") {

    # Cubature
    simMale <- apply(simMaleLin, 1:2, function(x){
                         hcubature(integrate2dRandomEffects,
                                   lowerLimit = -limOutcome,
                                   upperLimit = limOutcome, 
                                   linPred = x,
                                   sigma = sigmaOutcome,
                                   mean = c(0, 0),
                                   logdet = logdetOutcome,
                                   tol=1e-4,
                                   vectorInterface = TRUE)$integral})

    simFemale <- apply(simFemaleLin, 1:2, function(x){
                         hcubature(integrate2dRandomEffects,
                                   lowerLimit = -limOutcome,
                                   upperLimit = limOutcome, 
                                   linPred = x,
                                   sigma = sigmaOutcome,
                                   mean = c(0, 0),
                                   logdet = logdetOutcome,
                                   tol=1e-4,
                                   vectorInterface = TRUE)$integral})
  }

  # get quantiles of joR
  Hmale <- apply(simMale,1,quantile,c(alpha/2, 1-alpha/2))
  Hfemale <- apply(simFemale,1,quantile,c(alpha/2, 1-alpha/2))
  
  # return matrices of joR at specific skeleton ages and birthdays
  z <- list()
  
  z[[1]] <- list(l=matrix(Hmale[1,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          u=matrix(Hmale[2,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          mean=matrix(apply(simMale,1,mean), nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)))
          
  z[[2]] <- list(l=matrix(Hfemale[1,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          u=matrix(Hfemale[2,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          mean=matrix(apply(simFemale,1,mean), nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)))
            
  names(z) <- c("Male","Female")

  return(list(z=z, skelagePred=skelagePred, birthdayPred=birthdayPred, 
    alpha=alpha))
}
  
plotDerivative <- function(fit0,
                           sigmaOutcome = NULL,
                           logdetOutcome = NULL,
                           limOutcome = NULL,
                           skelagePred=seq(6,18,2),
                           birthdayPred=seq(1930,1995,0.5),
                           alpha=0.05,
                           B=10000,
                           eps=1e-7,
                           Bmc = 1000,
                           integrateRandomEffects = "none"){

  if(!integrateRandomEffects %in% c("none", "MC", "cubature")) {
    stop("integrateRandomEffects must be equal to 'none', 'MC', or 'cubature'")
  }

  # produces partial derivative of regression surface wrt birthday
  # via finite differences without random effects
  # with Bayesian credible intervals (conditional on smoothing parameter)
  
  fit <- fit0$gam
  
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

  newd0Female <- data.frame( skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred, each=nSkelage),
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

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

  newd1Female <- data.frame(skelage=rep(skelagePred, times=nBirthday),
            birthday=rep(birthdayPred+eps, each=nSkelage),
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  ## center body size
  vc <- VarCorr(b1$lme)

  sigmaPed <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
  sigmaPtno <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

  sigma <- diag(c(sigmaPed, sigmaPtno)^2)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))

  lim <- c(sigmaPed, sigmaPtno) * 100

  linPredOverall0 <- predict(b1$gam, newdata=newd0)
  linPredMale0 <- predict(b1$gam, newdata=newd0Male)
  linPredFemale0 <- predict(b1$gam, newdata=newd0Female)

  meanOverall0 <- sapply(linPredOverall0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanMale0 <- sapply(linPredMale0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanFemale0 <- sapply(linPredFemale0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  newd0Male$bodysizeCentered <- meanMale0 - meanOverall0
  newd0Female$bodysizeCentered <- meanFemale0 - meanOverall0

  linPredOverall1 <- predict(b1$gam, newdata=newd1)
  linPredMale1 <- predict(b1$gam, newdata=newd1Male)
  linPredFemale1 <- predict(b1$gam, newdata=newd1Female)

  meanOverall1 <- sapply(linPredOverall1, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanMale1 <- sapply(linPredMale1, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanFemale1 <- sapply(linPredFemale1, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  newd1Male$bodysizeCentered <- meanMale1 - meanOverall1
  newd1Female$bodysizeCentered <- meanFemale1 - meanOverall1


  X0MaleTemp <- predict(fit, newdata=newd0Male, type= "lpmatrix")
  X0FemaleTemp <- predict(fit, newdata=newd0Female, type= "lpmatrix")
  
  X0Male <- X0MaleTemp*0
  X0Male[,keepMale] <- X0MaleTemp[,keepMale]
  
  X0Female <- X0FemaleTemp*0
  X0Female[,keepFemale] <- X0FemaleTemp[,keepFemale]
  
  XdMaleTemp <- predict(fit, newdata=newd1Male, type= "lpmatrix")
  XdFemaleTemp <- predict(fit, newdata=newd1Female, type= "lpmatrix")
  
  XdMale <- XdMaleTemp*0
  XdMale[,keepMale] <- XdMaleTemp[,keepMale]
  
  XdFemale <- XdFemaleTemp*0
  XdFemale[,keepFemale] <- XdFemaleTemp[,keepFemale]

  # simulate from posterior for beta
  br <- rmvn(B, beta, fit$Vp)
  
  # get simulated values of joR
  if (integrateRandomEffects == "none") {
    # use finite differences to approximate partial derivative wrt birth year
    simMale <- (exp(XdMale %*% br) - exp(X0Male %*% br))/eps
    simFemale <- (exp(XdFemale %*% br) - exp(X0Female %*% br))/eps

  } else if (integrateRandomEffects == "MC") {

    simMaleLin0 <- X0Male %*% br
    simFemaleLin0 <- X0Female %*% br

    simMaleLin1 <- XdMale %*% br
    simFemaleLin1 <- XdFemale %*% br

    # MC approach  
    bSum <- colSums(rmvn(n = Bmc, mu = c(0, 0), sig = sigmaOutcome))

    simMale <- array(0, dim = dim(simMaleLin0))
    simFemale <- array(0, dim = dim(simFemaleLin0))

    for (iter in 1:Bmc) {
      simMale <- simMale + (exp(simMaleLin1 + bSum[iter]) - exp(simMaleLin0 + bSum[iter])) / (eps * Bmc)
      simFemale <- simFemale + (exp(simFemaleLin1 + bSum[iter]) - exp(simFemaleLin0 + bSum[iter])) / (eps * Bmc)
    }

  } else if (integrateRandomEffects == "cubature") {
    # Cubature
    for (i in 1:nrow(simMale)) {
      for (j in 1:ncol(simMale)) {
        simMale[i, j] <- hcubature(integrate2dRandomEffectsDer,
                                   lowerLimit = -limOutcome,
                                   upperLimit = limOutcome, 
                                   linPred0 = simMaleLin0[i, j],
                                   linPred1 = simMaleLin1[i, j],
                                   eps = eps,
                                   sigma = sigmaOutcome,
                                   mean = c(0, 0),
                                   logdet = logdetOutcome,
                                   tol=1e-4,
                                   vectorInterface = TRUE)$integral
      }
    }

    for (i in 1:nrow(simFemale)) {
      for (j in 1:ncol(simFemale)) {
        simFemale[i, j] <- hcubature(integrate2dRandomEffectsDer,
                                   lowerLimit = -limOutcome,
                                   upperLimit = limOutcome, 
                                   linPred0 = simFemaleLin0[i, j],
                                   linPred1 = simFemaleLin1[i, j],
                                   eps = eps,
                                   sigma = sigmaOutcome,
                                   mean = c(0, 0),
                                   logdet = logdetOutcome,
                                   tol=1e-4,
                                   vectorInterface = TRUE)$integral
      }
    }
  }

  # get quantiles for credible intervals
  Hmale <- apply(simMale,1,quantile,c(alpha/2, 1-alpha/2))
  Hfemale <- apply(simFemale,1,quantile,c(alpha/2, 1-alpha/2))
  
  z <- list()
  
  z[[1]] <- list(l=matrix(Hmale[1,], nrow=nSkelage, 
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          u=matrix(Hmale[2,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          mean=matrix(apply(simMale,1,mean), nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)))
          
  z[[2]] <- list(l=matrix(Hfemale[1,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          u=matrix(Hfemale[2,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          mean=matrix(apply(simFemale,1,mean), nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)))
          
  names(z) <- c("Male","Female")

  return(list(z=z, skelagePred=skelagePred, birthdayPred=birthdayPred, 
    alpha=alpha))
}

plotPerChange <- function(fit0,
                          initialBirthday=1930,
                          endBirthday=1995,
                          skelagePred=seq(8,18,.1), 
                          alpha=0.01,
                          B=10000,
                          ttar=FALSE){
  # produces percent change

  fit <- fit0$gam
  
  nSkelage <- length(skelagePred)
  
  beta <- coef(fit)  
  
  basePlusMaleInd <- grep("skelage",names(beta))
  maleInd <- grep("male",names(beta))
  femaleInd <- basePlusMaleInd[-which(basePlusMaleInd %in% maleInd)]

  keepMale <- c(1,basePlusMaleInd)
  keepFemale <- c(1,femaleInd)

  # initial timepoint 0
  newd0 <- data.frame(skelage=skelagePred,
            birthday=initialBirthday,
            male=percMale
            )
            
  newd0Male <- data.frame(skelage=skelagePred,
            birthday=initialBirthday,
            male=1,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Male")])
            )

  newd0Female <- data.frame(skelage=skelagePred,
            birthday=initialBirthday,
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  ## center body size
  vc <- VarCorr(b1$lme)

  sigmaPed <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
  sigmaPtno <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

  sigma <- diag(c(sigmaPed, sigmaPtno)^2)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))

  lim <- c(sigmaPed, sigmaPtno) * 100

  linPredOverall0 <- predict(b1$gam, newdata=newd0)
  linPredMale0 <- predict(b1$gam, newdata=newd0Male)
  linPredFemale0 <- predict(b1$gam, newdata=newd0Female)

  meanOverall0 <- sapply(linPredOverall0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanMale0 <- sapply(linPredMale0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanFemale0 <- sapply(linPredFemale0, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  newd0Male$bodysizeCentered <- meanMale0 - meanOverall0
  newd0Female$bodysizeCentered <- meanFemale0 - meanOverall0

  ## Get design matrix for making predictions
  X0MaleTemp <- predict(fit, newdata=newd0Male, type= "lpmatrix")
  X0FemaleTemp <- predict(fit, newdata=newd0Female, type= "lpmatrix")
  
  X0Male <- X0MaleTemp*0
  X0Male[,keepMale] <- X0MaleTemp[,keepMale]

  X0Female <- X0FemaleTemp*0  
  X0Female[,keepFemale] <- X0FemaleTemp[,keepFemale]
  
  # subsequent time point 1
  newd1 <- data.frame(skelage=skelagePred,
            birthday=endBirthday,
            male=percMale
            )
            
  newd1Male <- data.frame(skelage=skelagePred,
            birthday=endBirthday,
            male=1,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Male")])
            )

  newd1Female <- data.frame(skelage=skelagePred,
            birthday=endBirthday,
            male=0,
            phvageCentered=as.numeric(phvageCentSex[which(names(phvageCentSex)=="Female")])
            )

  linPredOverall1 <- predict(b1$gam, newdata=newd1)
  linPredMale1 <- predict(b1$gam, newdata=newd1Male)
  linPredFemale1 <- predict(b1$gam, newdata=newd1Female)

  meanOverall1 <- sapply(linPredOverall1, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanMale1 <- sapply(linPredMale1, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  meanFemale1 <- sapply(linPredFemale1, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  newd1Male$bodysizeCentered <- meanMale1 - meanOverall1
  newd1Female$bodysizeCentered <- meanFemale1 - meanOverall1

  # Get design matrix for predictions
  XdMaleTemp <- predict(fit, newdata=newd1Male, type= "lpmatrix")
  XdFemaleTemp <- predict(fit, newdata=newd1Female, type= "lpmatrix")
  
  XdMale <- XdMaleTemp*0  
  XdMale[,keepMale] <- XdMaleTemp[,keepMale]
  
  XdFemale <- XdFemaleTemp*0  
  XdFemale[,keepFemale] <- XdFemaleTemp[,keepFemale]
 
  # simulate from posterior for beta
  br <- rmvn(B, beta, fit$Vp)
  
  # Calculate percent change
  simMale <- (exp((XdMale - X0Male) %*% br)-1)
  simFemale <- (exp((XdFemale - X0Female) %*% br)-1)

  Hmale <- apply(simMale,1,quantile,c(alpha/2, 1-alpha/2))
  Hfemale <- apply(simFemale,1,quantile,c(alpha/2, 1-alpha/2))
  
  z <- list()
  
  z[[1]] <- data.frame(skelage=skelagePred,
          l=Hmale[1,],
          u=Hmale[2,],
          mean=apply(simMale,1,mean),
          comparison="Male")
          
  z[[2]] <- data.frame(skelage=skelagePred,
          l=Hfemale[1,],
          u=Hfemale[2,],
          mean=apply(simFemale,1,mean),
          comparison="Female")
          
  names(z) <- c("Male","Female")
  
  zAll <- do.call(rbind,z)

  if(!ttar){
  ggAll <- ggplot(aes(x=skelage, y=mean*100), data=zAll)+
    geom_line()+
    geom_line(aes(y=l*100),linetype="dashed")+
    geom_line(aes(y=u*100),linetype="dashed")+
    theme_bw(16)+
    facet_wrap(~comparison)+
    labs(
      # title= paste("Percent change ", initialBirthday, " to ", endBirthday, ": ",(1-alpha)*100, "% credible intervals", sep=""),
      y=expression(paste("Percent change in BSI",sep="")), x="Skeletal age")+
    geom_hline(yintercept=0, color="red")+
    scale_x_continuous(breaks=seq(min(skelagePred), max(skelagePred), 2))+
    scale_y_continuous(breaks=seq(-40,40,5))

  } else if(ttar){
  ggAll <- ggplot(aes(x=skelage, y=mean*100), data=zAll)+
    geom_line()+
    geom_line(aes(y=l*100),linetype="dashed")+
    geom_line(aes(y=u*100),linetype="dashed")+
    theme_bw(16)+
    facet_wrap(~comparison)+
    labs(
      # title= paste("Percent change ", initialBirthday, " to ", endBirthday, ": ",(1-alpha)*100, "% credible intervals", sep=""),
      y="Percent change in TA", x="Skeletal age")+
    geom_hline(yintercept=0, color="red")+
    scale_x_continuous(breaks=seq(min(skelagePred), max(skelagePred), 2))+
    scale_y_continuous(breaks=seq(-40,40,5))
  }
  
  return(list(ggAll=ggAll, z=z, zAll=zAll, skelagePred=skelagePred,
          initialBirthday=initialBirthday,
          endBirthday=endBirthday, 
          simMale = simMale,
          simFemale = simFemale))
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

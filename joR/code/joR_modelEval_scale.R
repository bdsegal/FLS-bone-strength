# Evaluate the results of the joR models,
# which are fit in 'joR_modelFits_scale.R'.
# Produces summary plots of trends, derivatives, and percent change.

library(mgcv)

library(reshape2)
library(ggplot2)

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

# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/joR/code"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))

# load functions for obtaining predictions
# source(file.path(dataPrepPath,"predict_functions_2prod_centered.R"))

# evaluate models -------------------------------------------------------------

load("m1scale.Rdata")  
load("m2scale.Rdata")
load("m3scale.Rdata") # note: fit with optim
load("m4scale.Rdata")
load("m5scale.Rdata")

aic <- data.frame(model=1:5,
                  aic=c(AIC(m1scale$lme),
                        AIC(m2scale$lme),
                        AIC(m3scale$lme),
                        AIC(m4scale$lme),
                        AIC(m5scale$lme)))
round(aic)

data.frame(model = 1:5,
           adjRsq = round(c(summary(m1scale$gam)$r.sq,
                            summary(m2scale$gam)$r.sq,
                            summary(m3scale$gam)$r.sq,
                            summary(m4scale$gam)$r.sq,
                            summary(m5scale$gam)$r.sq), 2))

# variance components
vc <- VarCorr(m1scale$lme)
nvc <- nrow(vc)
vc[(nvc-4):nvc,]

vc <- VarCorr(m2scale$lme)
nvc <- nrow(vc)
vc[(nvc-5):nvc,]

vc <- VarCorr(m3scale$lme)
nvc <- nrow(vc)
vc[(nvc-2):nvc,]

vc <- VarCorr(m4scale$lme)
nvc <- nrow(vc)
vc[(nvc-3):nvc,]

vc <- VarCorr(m5scale$lme)
nvc <- nrow(vc)
vc[(nvc-5):nvc,]

# other diagnostics -----------------------------------------------------------
beta <- cbind(coef(m1scale$gam),
              coef(m2scale$gam),
              coef(m3scale$gam),
              coef(m4scale$gam),
              coef(m5scale$gam))

matplot(beta)
matplot(log(abs(beta),10))

plot(m1scale$gam, scheme=1)
plot(m2scale$gam, scheme=1)
plot(m3scale$gam, scheme=1)
plot(m4scale$gam, scheme=1)
plot(m5scale$gam, scheme=1)

concurvity(m1scale$gam)   
concurvity(m2scale$gam)  
concurvity(m3scale$gam)  
concurvity(m4scale$gam)  
concurvity(m5scale$gam)  

gam.check(m1scale$gam) 
gam.check(m2scale$gam) 
gam.check(m3scale$gam) 
gam.check(m4scale$gam) 
gam.check(m5scale$gam) 

plot(m1scale$lme)
plot(m2scale$lme)
plot(m3scale$lme)
plot(m4scale$lme)
plot(m5scale$lme)

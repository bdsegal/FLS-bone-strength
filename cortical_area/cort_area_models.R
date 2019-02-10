
library(mgcv)
library(ggplot2)
library(nlme)
library(rgl)

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

paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/plots_report/all")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))

# functions -------------------------------------------------------------------
Rsq <- function(model, cortArea = dataSub$cortArea, type = "fixed") {
  # This function computes R^2 for lme objects, both with the fixed
  # effects component only (type = 'fixed') and with the BLUPS (type = 'random')

  if (class(model) != "lme") {
    stop("model must be an lme class object")
  }

  if (! type %in% c("fixed", "random")) {
    stop("type must be either 'fixed' or 'random'")
  }
  
  if (type == "fixed") {
    sse <- sum((cortArea - predict(model, level = 0))^2)
  } else if (type == "random") {
    sse <- sum((cortArea - predict(model))^2)
  }
  sst <- sum((cortArea - mean(cortArea))^2)
  
  return(1 - sse/sst)
}

# fit models ------------------------------------------------------------------
m1 <- lme(cortArea ~ birthday + ttar + sex + skelage + bodysizeCentered, 
           random = ~ 1|pedno/ptno,
           data = dataSub)

m2 <- lme(cortArea ~ birthday + ttar + sex + skelage*bodysizeCentered, 
           random = ~ 1|pedno/ptno,
           data = dataSub)

m3 <- lme(cortArea ~ birthday + ttar + sex*skelage*bodysizeCentered, 
           random = ~ 1|pedno/ptno,
           data = dataSub)

m4 <- lme(cortArea ~ birthday + ttar*sex*skelage*bodysizeCentered, 
           random = ~ 1|pedno/ptno,
           data = dataSub)

m5 <- lme(cortArea ~ birthday*ttar*sex*skelage*bodysizeCentered, 
           random = ~ 1|pedno/ptno,
           data = dataSub)

m6 <- lme(cortArea ~ birthday + ttar + sex:skelage + bodysizeCentered, 
           random = ~ 1|pedno/ptno,
           data = dataSub)

m7 <- lme(cortArea ~ birthday + ttar + sex:skelage, 
           random = ~ 1|pedno/ptno,
           data = dataSub)

plot(m1)
plot(m2)
plot(m3)
png(file.path(paperPath, "ca_residuals.png"))
plot(m4)
dev.off()
plot(m5)
plot(m6)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)
summary(m7)

# variance components
VarCorr(m1)
VarCorr(m2)
VarCorr(m3)
VarCorr(m4)
VarCorr(m5)
VarCorr(m6)
VarCorr(m7)

mse <- c(mean((predict(m1, level = 0) - dataSub$cortArea)^2),
         mean((predict(m2, level = 0) - dataSub$cortArea)^2),
         mean((predict(m3, level = 0) - dataSub$cortArea)^2),
         mean((predict(m4, level = 0) - dataSub$cortArea)^2),
         mean((predict(m5, level = 0) - dataSub$cortArea)^2),
         mean((predict(m6, level = 0) - dataSub$cortArea)^2),
         mean((predict(m7, level = 0) - dataSub$cortArea)^2))

fitTable <- data.frame(model = paste0("m", 1:7),
                       aic = round(c(AIC(m1),
                                     AIC(m2),
                                     AIC(m3),
                                     AIC(m4),
                                     AIC(m5),
                                     AIC(m6),
                                     AIC(m7))),
                       Rsq = signif(c(Rsq(m1),
                                      Rsq(m2),
                                      Rsq(m3),
                                      Rsq(m4),
                                      Rsq(m5),
                                      Rsq(m6),
                                      Rsq(m7)), 2),
                       mse = signif(mse, 2))
fitTable

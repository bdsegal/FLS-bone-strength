library(cubature)

# Prepare data for regressions 
if("dataSub.Rdata" %in% list.files(file.path(path))) {
  load(file.path(path, "dataSub.Rdata"))

} else {

  data <- read.csv(file.path(path,"skeleton_long_INCLUDE.csv"),na.strings=".")

  data$sex[which(data$sex==1)] <- "Male"
  data$sex[which(data$sex==2)] <- "Female"
  data$male <- ifelse(data$sex=="Male",1,0)

  data$joR <- with(data, jo/(mc2tmd/2))

  data$ptno <- as.factor(data$ptno)
  data$pedno <- as.factor(data$pedno)

  # fix missing skelage ------------------------------------------------

  skelageToFix <- which(with(data, is.na(skelage) & age>=18))
  length(skelageToFix)
  # [1] 52

  data[skelageToFix,"skelage"] <- 18

  # keep only those above target age 8
  # (restricting analysis to ages for which BSI was collected across all birthdays)
  dataSub <- data[which(data$targage>=8),]

  # keeping only complete cases makes it easier to add predicted values to the dataset
  keep <- which(complete.cases(dataSub[,c("skelage","birthday","bodysize","male", "ptno",
    "pedno")]))

  # number of observations removed (assumed missing completely at random)
  nrow(dataSub) - length(keep)
  # [1] 37

  # view subjects with missing data
  # remove <- which(!complete.cases(dataSub[,c("skelage","birthday","bodysize","male", "ptno",
  #   "pedno")]))
  # dataSub[remove,]

  dataSub <- dataSub[keep,]

  # center bodysize on mean trajectory ----------------------------
  # Note: this script also includes code for standardizing and scaling bodysize,
  # and for standardizing phv age. In the final analysis, we did not used the
  # standardized or scaled versions.

  if("b1.Rdata" %in% list.files(file.path(computer,"Dropbox/Research/Bones/final_analysis"))) {
    load(file=file.path(computer,"Dropbox/Research/Bones/final_analysis", "b1.Rdata"))
  } else {
    b1 <- gamm(bodysize ~ s(skelage, k=25)+
      s(skelage, k=25,  by=male), 
      family=Gamma(link=log),
      data=dataSub,
      random=list(pedno=~1, ptno=~1))
    save(b1, file=file.path(computer, "Dropbox/Research/Bones/final_analysis", "b1.Rdata"))
  }

  # get percent male to make population average predictions
  # using percent of observations, not subjects
  sexTable <- table(dataSub[,"sex"])
  percMale <- as.numeric(sexTable[which(names(sexTable) == "Male")]/sum(sexTable))

  linPredOverall <- predict(b1$gam, newdata=data.frame(skelage=dataSub$skelage, male=percMale))

  vc <- VarCorr(b1$lme)

  sigmaPed <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
  sigmaPtno <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

  sigma <- diag(c(sigmaPed, sigmaPtno)^2)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))

  lim <- c(sigmaPed, sigmaPtno) * 100

  dataSub$bodysizeMean <- sapply(linPredOverall, function(x){
                       hcubature(integrate2dRandomEffects,
                                 lowerLimit = -lim,
                                 upperLimit = lim, 
                                 linPred = x,
                                 sigma = sigma,
                                 mean = c(0, 0),
                                 logdet = logdet,
                                 tol=1e-4,
                                 vectorInterface = TRUE)$integral})

  dataSub$bodysizeCentered <- with(dataSub, bodysize - bodysizeMean)
  dataSub$bodysizeScaled <- with(dataSub, bodysizeCentered/500)

  dataSub$phvageCentered <- with(dataSub, phvage - mean(phvage, na.rm=TRUE))

  # create cortical area 
  dataSub$R <- dataSub$mc2tmd / 2
  dataSub$r <- with(dataSub, (R^4 - 2 / pi * R * joR)^(1/4))
  dataSub$cortArea <- with(dataSub, pi * (R^2 - r^2))


  save(dataSub, file=file.path(path,"dataSub.Rdata"))
}

# load into global environment
sexTable <- table(dataSub[,"sex"])
percMale <- as.numeric(sexTable[which(names(sexTable) == "Male")]/sum(sexTable))

phvageMeanSex <- tapply(dataSub$phvage, dataSub$sex, mean, na.rm=TRUE)
phvageCentSex <- phvageMeanSex-mean(dataSub$phvage, na.rm=TRUE)

args <- commandArgs(TRUE)

library(mgcv)

path <- file.path(dirname(dirname(getwd())),"data")

# prep data
load(file.path(dirname(dirname(getwd())),"data", "dataSub.Rdata"))

# indicator variable for missing, and indices of missingness
dataSub$phvageMiss <- is.na(dataSub$phvage)
phvageNA <- which(dataSub$phvageMiss)

fitPhv <- lm(phvage ~ birthday*male, data=dataSub)

# get variance of predicted phvage
sigma <- summary(fitPhv)$sigma
bd <- dataSub$birthday[!dataSub$phvageMiss]
mal <- dataSub$male[!dataSub$phvageMiss]
X <- cbind(1, bd, mal, bd*mal)
XtX <- solve(t(X)%*%X)

# variance of predictions
varPred <- rep(NA, nrow(X))
for (i in 1:length(varPred)){
  varPred[i] <- sigma^2 * (1 + t(X[i,]) %*% XtX %*% X[i,])
}

# mean of predicted phvage
mu <- predict(fitPhv, newdata=dataSub[phvageNA,])

# averaging over phvage_miss
M <- 50

phvCent <- matrix(nrow=M, ncol=2)
multImpute <- list()

for (m in 1:M){

  print(m)
  # impute phvage as draws from posterior
  dataSub$phvage[phvageNA] <- rnorm(n=length(phvageNA), mean=mu, sd=sqrt(varPred))

  dataSub$phvageCentered <- with(dataSub, phvage - mean(phvage))
  
  phvageMeanSex <- tapply(dataSub$phvage, dataSub$sex, mean, na.rm=TRUE)
  phvageCentSex <- phvageMeanSex-mean(dataSub$phvage, na.rm=TRUE)

  # adding phvage to model m1cent
  try({
    fit <- gamm(mc2tci ~ 
          te(skelage, birthday, k=c(5,5), bs="cr")+
          te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
          te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
          te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr")+
          te(skelage, birthday, k=c(5,5), by=phvageCentered, bs="cr")+
          te(skelage, birthday, k=c(5,5), by=male*phvageCentered, bs="cr"),
          random=list(pedno=~1, ptno=~1),
          data=dataSub)
          
  multImpute[[m]] <- fit$gam
  phvCent[m,] <- phvageCentSex
  })
  
}

# keep iterations with models that converged
keep <- which(!sapply(multImpute, is.null))
multImpute <- multImpute[keep]
phvCent <- phvCent[keep,]

save(multImpute, phvCent, file=paste("imputeList_cti_",args,".Rdata",sep=""))

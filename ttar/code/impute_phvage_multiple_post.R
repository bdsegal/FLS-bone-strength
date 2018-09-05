# Process results from the multiple imputation models.
# The matrices with posterior estimates are about 4 GB each
# using the bigmemory package.
# First, need to fit models with the 'impute_phvage_multiple.R' file
# in the 'imputations' folder.

library(mgcv)
library(bigmemory)
library(ggplot2)

# set paths and load functions and data ----------------------------------------

# set computer-specific paths (note: 'path' used in dataPred.R)
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

paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/plots_report/all")

setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/ttar"))
resultsPath <- file.path(getwd(), "imputations")
files <- grep("Rdata",list.files(resultsPath), value=TRUE)

dataPrepPath <- dirname(getwd())

# prep data (needed for percMale)
source(file.path(dataPrepPath,"data_prep.R"))

# load functions for obtaining predictions
source(file.path(dataPrepPath,
  "predict_functions_2prod_centered_impute_matrix_fill.R"))

# set values for obtaining posterior distributions and making plots ------------

if(!"numModels.Rdata" %in% list.files()) {
  numModels <- rep(NA, length(files))
  for (i in 1:length(files)){
    print(i)
    load(file.path(resultsPath, files[i]))
    numModels[i] <- length(multImpute)
  }
  save(numModels,file="numModels.Rdata")
} else {
  load("numModels.Rdata")
}
sum(numModels)
# 591

skelagePred <- seq(8,18, 2)
skelagePredPc <- seq(8,18,.1)
birthdayPred <- seq(1930,2000, 0.5)
nSkelage <- length(skelagePred)

B=1000
p <- length(skelagePred)*length(birthdayPred)
n <- sum(numModels)*B

initialBirthday <- 1930
endBirthday <- 2000


alpha <- 0.05

# matrices for storing posterior estimates (too big for RAM)
trendMale <- big.matrix(nrow=n,
                  ncol=p,
                  init=0,
                  backingfile="trendMale.bin",
                  descriptorfile="trendMale.desc")

trendFemale <- big.matrix(nrow=n,
                  ncol=p,
                  init=0,
                  backingfile="trendFemale.bin",
                  descriptorfile="trendFemale.desc")

derivMale <- big.matrix(nrow=n,
                  ncol=p,
                  init=0,
                  backingfile="derivMale.bin",
                  descriptorfile="derivMale.desc")

derivFemale <- big.matrix(nrow=n,
                  ncol=p,
                  init=0,
                  backingfile="derivFemale.bin",
                  descriptorfile="derivFemale.desc")

pcMale <- big.matrix(nrow=n,
                  ncol=length(skelagePredPc),
                  init=0,
                  backingfile="pcMale.bin",
                  descriptorfile="pcMale.desc")

pcFemale <- big.matrix(nrow=n,
                  ncol=length(skelagePredPc),
                  init=0,
                  backingfile="pcFemale.bin",
                  descriptorfile="pcFemale.desc")

# fill in posterior matrices
count <- 1
for (i in 1:length(files)){

  load(file.path(resultsPath, files[i]))
  M <- length(multImpute)
  
  for (m in 1:M){
    print(paste("file ", i, " model ", m, sep=""))

    phvageCentSex <- phvCent[m,]
    names(phvageCentSex) <- c("Female","Male")

    trendTemp <- plotTrend1(multImpute[[m]], B=B,
      skelagePred=skelagePred,
      birthdayPred=birthdayPred)

    trendMale[((count-1)*B +1):(count*B),] <- t(trendTemp$simMale)
    trendFemale[((count-1)*B +1):(count*B),] <- t(trendTemp$simFemale)

    derivTemp <- plotDerivative1(multImpute[[m]], B=B,
      skelagePred=skelagePred,
      birthdayPred=birthdayPred)

    derivMale[((count-1)*B +1):(count*B),] <- t(derivTemp$simMale)
    derivFemale[((count-1)*B +1):(count*B),] <- t(derivTemp$simFemale)

    pcTemp <- plotPerChange1(multImpute[[m]], B=B, skelagePred=skelagePredPc,
        initialBirthday=initialBirthday, endBirthday=endBirthday)

    pcMale[((count-1)*B +1):(count*B),] <- t(pcTemp$simMale)
    pcFemale[((count-1)*B +1):(count*B),] <- t(pcTemp$simFemale)

    count <- count + 1
  }
}

# get quantiles for trend ------------------------------------------------------
trendMale <- attach.big.matrix("trendMale.desc")
trendFemale <- attach.big.matrix("trendFemale.desc")

Hmale <- matrix(nrow=2, ncol=ncol(trendMale))
Hfemale <- matrix(nrow=2, ncol=ncol(trendFemale))
meanMale <- rep(NA, ncol(trendMale))
meanFemale <- rep(NA, ncol(trendFemale))

for (i in 1:ncol(trendMale)){
  print(i)
  Hmale[,i] <- quantile(trendMale[,i], probs=c(alpha/2, 1-alpha/2))
  Hfemale[,i] <- quantile(trendFemale[,i], probs=c(alpha/2, 1-alpha/2))
  meanMale[i] <- mean(trendMale[,i])
  meanFemale[i] <- mean(trendFemale[,i])
}

z <- list()
  
  z[[1]] <- list(l=matrix(Hmale[1,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          u=matrix(Hmale[2,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          mean=matrix(meanMale, nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)))
          
  z[[2]] <- list(l=matrix(Hfemale[1,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          u=matrix(Hfemale[2,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          mean=matrix(meanFemale, nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)))
            
  names(z) <- c("Male","Female")

out <- list(z=z, skelagePred=skelagePred, birthdayPred=birthdayPred,
  alpha=alpha)

# 2d plots: mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="Total area", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,70))
ggsave(file.path(paperPath, "m1ttarCentPhvageMultImputeMaleMean.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="Total area", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,70))
ggsave(file.path(paperPath, "m1ttarCentPhvageMultImputeFemaleMean.png"))

# get quantiles for derivative -------------------------------------------------
derivMale <- attach.big.matrix("derivMale.desc")
derivFemale <- attach.big.matrix("derivFemale.desc")

Hmale <- matrix(nrow=2, ncol=ncol(derivMale))
Hfemale <- matrix(nrow=2, ncol=ncol(derivFemale))
meanMale <- rep(NA, ncol(derivMale))
meanFemale <- rep(NA, ncol(derivFemale))

for (i in 1:ncol(derivMale)){
  print(i)
  Hmale[,i] <- quantile(derivMale[,i], probs=c(alpha/2, 1-alpha/2))
  Hfemale[,i] <- quantile(derivFemale[,i], probs=c(alpha/2, 1-alpha/2))
  meanMale[i] <- mean(derivMale[,i])
  meanFemale[i] <- mean(derivFemale[,i])
}

z <- list()
  
  z[[1]] <- list(l=matrix(Hmale[1,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          u=matrix(Hmale[2,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          mean=matrix(meanMale, nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)))
          
  z[[2]] <- list(l=matrix(Hfemale[1,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          u=matrix(Hfemale[2,], nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)),
          mean=matrix(meanFemale, nrow=nSkelage,
            dimnames=list(skelage=skelagePred, birthday=birthdayPred)))
            
  names(z) <- c("Male","Female")

outd <- list(z=z, skelagePred=skelagePred, birthdayPred=birthdayPred,
  alpha=alpha)

# 2d plots: derivatives
slice <- sliceFun(outd, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    geom_hline(yintercept=0, color="red")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y=expression(paste(partialdiff," Total area / ", 
          partialdiff," birthdate", sep="")), 
        x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.4, 0.2))
ggsave(file.path(paperPath, "m1ttarCentPhvageMultImputeMaleDer.png"))

slice <- sliceFun(outd, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(16)+
    geom_hline(yintercept=0, color="red")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y=expression(paste(partialdiff," Total area / ", 
          partialdiff," birthdate", sep="")), 
        x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.4, 0.2))
ggsave(file.path(paperPath, "m1ttarCentPhvageMultImputeFemaleDer.png"))

# get quantities for percent change --------------------------------------------
pcMale <- attach.big.matrix("pcMale.desc")
pcFemale <- attach.big.matrix("pcFemale.desc")

Hmale <- matrix(nrow=2, ncol=ncol(pcMale))
Hfemale <- matrix(nrow=2, ncol=ncol(pcFemale))
meanMale <- rep(NA, ncol(pcMale))
meanFemale <- rep(NA, ncol(pcFemale))

for (i in 1:ncol(pcMale)){
  print(i)
  Hmale[,i] <- quantile(pcMale[,i], probs=c(alpha/2, 1-alpha/2))
  Hfemale[,i] <- quantile(pcFemale[,i], probs=c(alpha/2, 1-alpha/2))
  meanMale[i] <- mean(pcMale[,i])
  meanFemale[i] <- mean(pcFemale[,i])
}

z <- list()
  
z[[1]] <- data.frame(skelage=skelagePredPc,
        l=Hmale[1,],
        u=Hmale[2,],
        mean=meanMale,
        comparison="Male")
          
z[[2]] <- data.frame(skelage=skelagePredPc,
        l=Hfemale[1,],
        u=Hfemale[2,],
        mean=meanFemale,
        comparison="Female")
        
names(z) <- c("Male","Female")

zAll <- do.call(rbind,z)

ggAll <- ggplot(aes(x=skelage, y=mean*100), data=zAll)+
  geom_line()+
  geom_line(aes(y=l*100),linetype="dashed")+
  geom_line(aes(y=u*100),linetype="dashed")+
  theme_bw(16)+
  facet_wrap(~comparison)+
  labs(
    # title= paste("Percent change ", initialBirthday, " to ", endBirthday, ": ",(1-alpha)*100, "% credible intervals", sep=""),
    y="Percent change in\ntotal area", x="Skeletal age")+
  geom_hline(yintercept=0,color="red")+
  scale_x_continuous(breaks=seq(min(skelagePredPc), max(skelagePredPc), 2))+
  scale_y_continuous(breaks=seq(-30,30,5))

dev.new(height=3.5, width=7.25)
ggAll
ggsave(file.path(paperPath, "m1ttarCentPhvageMultImputepc.png"))

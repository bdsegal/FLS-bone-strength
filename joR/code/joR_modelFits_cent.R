# Fit several models for joR with centered bodysize,
# with 1) different random effect structures without phvage and 
# 2) with phvage and various sample restrictions

library(mgcv)
library(dplyr)

# set computer-specific paths (note: 'path' points to the folder
# containing the data, and is used in 'data_prep.R')
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

# path for paper
paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/joR/paper")
# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/joR/code"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# load functions for obtaining predictions and integrals
source(file.path(dataPrepPath,"predict_functions_2prod_centered.R"))

# prep/load data
source(file.path(dataPrepPath,"data_prep.R"))

# save analytic dataset for main analyses
main_analysis_data <- dataSub %>%
  select(pedno, ptno, joR, skelage, birthday, male, bodysizeCentered)

write.csv(main_analysis_data,
          row.names = FALSE, 
          file = file.path(path, "main_analysis_data.cvs"))

# 1) Different random effect structures ---------------------------------------

# random intercept for subject and pedigree
# k=5 basis functions for skelage and birthday (for quicker fitting)

m1cent <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub)

save(m1cent, file="m1cent.Rdata")

m2cent <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1 + skelage),
				data=dataSub)
# Maximum number of PQL iterations:  20 
# iteration 1
# iteration 2
# Error in lme.formula(fixed = fixed, random = random, data = data, correlation = correlation,  : 
  # nlminb problem, convergence error code = 1
  # message = iteration limit reached without convergence (10)
       
# save(m2cent, file="m2cent.Rdata")  

m3cent <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(ptno=~1),
				data=dataSub)
        
save(m3cent, file="m3cent.Rdata")  

m4cent <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(ptno=~1 + skelage),
				data=dataSub)
# Maximum number of PQL iterations:  20 
# iteration 1
# iteration 2
# Error in lme.formula(fixed = fixed, random = random, data = data, correlation = correlation,  : 
  # nlminb problem, convergence error code = 1
  # message = false convergence (8)
        
# save(m4cent, file="m4cent.Rdata")  

m5cent <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1 + skelage, ptno=~1),
				data=dataSub)
# Maximum number of PQL iterations:  20 
# iteration 1
# Error in lme.formula(fixed = fixed, random = random, data = data, correlation = correlation,  : 
  # nlminb problem, convergence error code = 1
  # message = iteration limit reached without convergence (10)
        
# save(m5cent, file="m5cent.Rdata")  
  
m6cent <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1 + skelage, ptno=~1 + skelage),
				data=dataSub)
 # Maximum number of PQL iterations:  20 
# iteration 1
# Error in lme.formula(fixed = fixed, random = random, data = data, correlation = correlation,  : 
  # nlminb problem, convergence error code = 1
  # message = iteration limit reached without convergence (10)

m7cent <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1 + skelage + I(skelage^2)),
				data=dataSub)
# Maximum number of PQL iterations:  20 
# iteration 1
# Error in lme.formula(fixed = fixed, random = random, data = data, correlation = correlation,  : 
  # nlminb problem, convergence error code = 1
  # message = iteration limit reached without convergence (10)

# Same as above, but fitting random effects as though they were
# penalized components
# DID NOT FINISH RUNNING - CRASHED R
# m0Gamcent <- gam(joR ~ 
				# te(skelage, birthday, k=c(5,5), bs="cr")+
				# te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				# te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				# te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr")+
        # s(ptno, bs="re")+
        # s(pedno, bs="re"),
				# family=Gamma(link=log),
				# # random=list(ptno=~1, pedno=~1),
				# gamma=1.4,
        # # method="REML"),
				# data=dataSub)

# m3 with gamm4
# took LONG time to run
library(gamm4)
mgamm4cent <- gamm4(joR ~ 
				t2(skelage, birthday, k=c(5,5), bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random = ~ (1|ptno),
				data=dataSub)
        
# 2) With phvage and various sample restrictions ------------------------------

# adding phvage to model m1cent
m1centPhvage <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=phvageCentered, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=male*phvageCentered, bs="cr"),
				family=Gamma(link=log),
        control=list(opt='optim'),
				random=list(pedno=~1, ptno=~1),
				data=dataSub)
        
save(m1centPhvage, file="m1centPhvage.Rdata")
     
m1centNoNAPhvage <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub[which(!is.na(dataSub$phvageCentered)),])
save(m1centNoNAPhvage, file="m1centNoNAPhvage.Rdata")

m1centPre1995 <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub[which(dataSub$birthday <=1995),])
save(m1centPre1995, file="m1centPre1995.Rdata")

m1centPre1990 <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub[which(dataSub$birthday <=1990),])
save(m1centPre1990, file="m1centPre1990.Rdata")

  dataSub$post1990andNotMiss <- with(dataSub, birthday>=1990 & !is.na(phvage))
  sub <- dataSub[dataSub$post1990andNotMiss,]

  nrow(sub)
  # [1] 226
  length(unique(sub$ptno))
  # [1] 29

  dataSub$post1990 <- with(dataSub, birthday>=1990)
  sub <- dataSub[dataSub$post1990,]

  nrow(sub)
  length(unique(sub$ptno))

dataSub$remove <- with(dataSub, birthday >=1990 & birthday <=1995)
  sub <- dataSub[dataSub$remove,]
  nrow(sub)
  length(unique(sub$ptno))

m1centRemove1990_1995 <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub[which(!dataSub$remove),])
save(m1centRemove1990_1995, file="m1centRemove1990_1995.Rdata")

# pedigrees entering the study on or before 1950
dataSub <- as.data.frame(mutate(group_by(dataSub,pedno),
	firstPedBY = min(birthday))
)

max(dataSub[which(dataSub$firstPedBY <=1950),"birthday"])
sub <- dataSub[which(dataSub$firstPedBY <=1950),]
nrow(sub)
length(unique(sub$ptno))

m1centPedPre1950 <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub[which(dataSub$firstPedBY <=1950),])
save(m1centPedPre1950, file="m1centPedPre1950.Rdata")

# removing outliers
outliers <- unique(dataSub[which(dataSub$joR >=230),"ptno"])
dataSub[which(dataSub$joR >=230),]
m1centNoOutliers <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub[which(! dataSub$ptno %in% outliers),])
save(m1centNoOutliers, file="m1centNoOutliers.Rdata")

# restricting BMI
bmiGrowth <- read.csv("http://www.cdc.gov/growthcharts/data/zscore/bmiagerev.csv",
  stringsAsFactors=FALSE)
# there is one row of text in the middle of the file
# which the following code deals with
bmiGrowth <- as.data.frame(sapply(bmiGrowth, as.numeric))
bmiGrowth <- bmiGrowth[-which(is.na(bmiGrowth[,1])),]

bmiGrowth$age <- bmiGrowth$Agemos/12

# 1=male; 2=female
dataSub$bmi5 <- NA
dataSub$bmi95 <- NA
for (i in 1:nrow(dataSub)){
  sex <- ifelse(dataSub$sex[i]=="Female", 2, 1)
  bmiSub <- bmiGrowth[which(bmiGrowth$Sex==sex),]
  row <- which.min(abs(dataSub$age[i] - bmiSub$age))
  dataSub$bmi5[i] <- bmiSub[row, "P5"]
  dataSub$bmi95[i] <- bmiSub[row, "P95"]
}

# note: 6 missing BMI values
dataSub$bmiObsFlag <- with(dataSub, (bmi <= bmi5 | bmi >= bmi95)*1)
dataSub[which(is.na(dataSub$bmiObsFlag)),]
max(dataSub[which(!dataSub$bmiObsFlag),"birthday"])

sub <- dataSub[which(!dataSub$bmiObsFlag),]
nrow(sub)
length(unique(sub$ptno))

m1centRestrictBMI <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub[which(dataSub$bmiObsFlag == 0),])
save(m1centRestrictBMI, file="m1centRestrictBMI.Rdata")

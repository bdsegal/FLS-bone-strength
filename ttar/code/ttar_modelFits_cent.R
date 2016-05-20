# Fit several models for ttar with centered bodysize,
# with different random effect structures without phvage

library(mgcv)

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
paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/ttar/paper")
# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/ttar/code"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))

# random intercept for subject and pedigree
# k=5 basis functions for skelage and birthday (for quicker fitting)

m1ttarCent <- gamm(ttar ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub)

save(m1ttarCent, file="m1ttarCent.Rdata")

m2ttarCent <- gamm(ttar ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1 + skelage),
				data=dataSub)
     
save(m2ttarCent, file="m2ttarCent.Rdata")  

m3ttarCent <- gamm(ttar ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(ptno=~1),
				data=dataSub)
        
save(m3ttarCent, file="m3ttarCent.Rdata")  

m4ttarCent <- gamm(ttar ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(ptno=~1 + skelage),
				data=dataSub)

save(m4ttarCent, file="m4ttarCent.Rdata")  


m5ttarCent <- gamm(ttar ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1 + skelage, ptno=~1),
				data=dataSub)
  
save(m5ttarCent, file="m5ttarCent.Rdata")  
  
m6ttarCent <- gamm(ttar ~ 
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
# save(m6ttarCent, file="m6ttarCent.Rdata")  

m7ttarCent <- gamm(ttar ~ 
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
# save(m7ttarCent, file="m7ttarCent.Rdata")  

# m3 with gamm4
# scaled covariates, took LONG time to run
library(gamm4)
mgamm4cent <- gamm4(ttar ~ 
				t2(skelage, birthday, k=c(5,5), bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random = ~ (1|ptno),
				data=dataSub)

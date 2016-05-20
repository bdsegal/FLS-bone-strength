# fit different random effect structures with scaled body size
# (had less convergence problems than with centered body size 
# when trying to fit random effect structures)

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
paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/joR/paper")
# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/joR/code"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# prep data
source(file.path(dataPrepPath,"data_prep_joR.R"))

# random intercept for subject and pedigree
# k=5 basis functions for skelage and birthday (for quicker fitting)

m1scale <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeScaled, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeScaled, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1),
				data=dataSub)

save(m1scale, file="m1scale.Rdata")

m2scale <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeScaled, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeScaled, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1 + skelage),
				data=dataSub)

save(m2scale, file="m2scale.Rdata")  

m3scale <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeScaled, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeScaled, bs="cr"),
				family=Gamma(link=log),
				random=list(ptno=~1),
        control=list(opt='optim'),
				data=dataSub)

save(m3scale, file="m3scale.Rdata")  

m4scale <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeScaled, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeScaled, bs="cr"),
				family=Gamma(link=log),
				random=list(ptno=~1 + skelage),
				data=dataSub)

save(m4scale, file="m4scale.Rdata")  

m5scale <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeScaled, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeScaled, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1 + skelage, ptno=~1),
				data=dataSub)

save(m5scale, file="m5scale.Rdata") 

m6scale <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeScaled, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeScaled, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1 + skelage, ptno=~1 + skelage),
				data=dataSub)
 # Maximum number of PQL iterations:  20 
# iteration 1
# iteration 2
# iteration 3
# Error in chol.default((value + t(value))/2) : 
  # the leading minor of order 5 is not positive definite

save(m6scale, file="m6scale.Rdata") 

m7scale <- gamm(joR ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeScaled, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeScaled, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1 + skelage + I(skelage^2)),
				data=dataSub)
 # Maximum number of PQL iterations:  20 
# iteration 1
# Error in lme.formula(fixed = fixed, random = random, data = data, correlation = correlation,  : 
  # nlminb problem, convergence error code = 1
  # message = iteration limit reached without convergence (10)
        
save(m7scale, file="m7scale.Rdata")

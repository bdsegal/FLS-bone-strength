# Fit several models for cti with centered bodysize,
# with different random effect structures

# note: 'tci' and 'cti' refer to the same variable
# it is labeled as tci in the dataset

library(mgcv)

# set computer-specific paths (note: 'path' used in data_prep.R)
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
paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/cti/paper")
# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/cti/code"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))

# random intercept for subject and pedigree
# k=5 basis functions for skelage and birthday (for quicker fitting)

m1tciCent <- gamm(mc2tci ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				random=list(pedno=~1, ptno=~1),
				data=dataSub)

save(m1tciCent, file="m1tciCent.Rdata")

m2tciCent <- gamm(mc2tci ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				random=list(pedno=~1, ptno=~1 + skelage),
				data=dataSub)
     
save(m2tciCent, file="m2tciCent.Rdata")  

m3tciCent <- gamm(mc2tci ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				random=list(ptno=~1),
				data=dataSub)
        
save(m3tciCent, file="m3tciCent.Rdata")  

m4tciCent <- gamm(mc2tci ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				random=list(ptno=~1 + skelage),
				data=dataSub)

save(m4tciCent, file="m4tciCent.Rdata")  

m5tciCent <- gamm(mc2tci ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				random=list(pedno=~1 + skelage, ptno=~1),
				data=dataSub)
  
save(m5tciCent, file="m5tciCent.Rdata")  
  
m6tciCent <- gamm(mc2tci ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				random=list(pedno=~1 + skelage, ptno=~1 + skelage),
				data=dataSub)
        
save(m6tciCent, file="m6tciCent.Rdata")  

m7tciCent <- gamm(mc2tci ~ 
				te(skelage, birthday, k=c(5,5), bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				family=Gamma(link=log),
				random=list(pedno=~1, ptno=~1 + skelage + I(skelage^2)),
				data=dataSub)
 # Maximum number of PQL iterations:  20 
# iteration 1
# Error in chol.default((value + t(value))/2) : 
  # the leading minor of order 5 is not positive definite       
# save(m7tciCent, file="m7tciCent.Rdata")

# m3 with gamm4
library(gamm4)
mgamm4tciCent <- gamm4(tci ~ 
				t2(skelage, birthday, k=c(5,5), bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=male, bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
				t2(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr"),
				random = ~ (1|ptno),
				data=dataSub)
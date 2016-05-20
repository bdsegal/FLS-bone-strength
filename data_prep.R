# Prepare data for regressions 

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

remove <- which(!complete.cases(dataSub[,c("skelage","birthday","bodysize","male", "ptno",
  "pedno")]))

# view subjects with missing data
# dataSub[remove,]

dataSub <- dataSub[keep,]

# get percent male to make population average predictions
# using percent of observations, not subjects
sexTable <- table(dataSub[,"sex"])
percMale <- as.numeric(sexTable[which(names(sexTable) == "Male")]/sum(sexTable))

# data for predictions
dataPred <- as.data.frame(cbind(skelage = dataSub$skelage, male=percMale))

# center bodysize on mean trajectory ----------------------------
# Note: this script also includes code for standardizing and scaling bodysize,
# and for standardizing phv age. In the final analysis, we did not used the
# standardized or scaled versions.

# Run model b1 once, then comment out
# b1 <- gamm(bodysize ~ s(skelage, k=25)+
#   s(skelage, k=25,  by=male), 
#   family=Gamma(link=log),
#   data=dataSub,
#   random=list(pedno=~1, ptno=~1))
# save(b1, file=file.path(computer, "Dropbox/Research/Bones/FinalAnalysis", "b1.Rdata"))
load(file=file.path(computer,"Dropbox/Research/Bones/final_analysis", "b1.Rdata"))

# get marginal predictions (not including random effects)
dataSub$bodysizePred <- exp(predict(b1$gam, newdata=dataPred))

# center and standardize body size
dataSub$bodysizeCentered <- with(dataSub, bodysize - bodysizePred)

# run b2 once, then comment out
# regress squared residuals on skelage, and get expected variance
# b2 <- gamm(bodysizeCentered^2 ~ s(skelage, k=25), 
    # random=list(pedno=~1, ptno=~1),
    # data=dataSub, family=Gamma(link=log))
# save(b2, file=file.path(computer,"Dropbox/Research/Bones/FinalAnalysis","b2.Rdata"))
load(file.path(computer,"Dropbox/Research/Bones/final_analysis", "b2.Rdata"))

dataSub$bodysizeVar <- predict(b2$gam, type="response")
dataSub$bodysizeStandard <- with(dataSub, bodysizeCentered/sqrt(bodysizeVar))

# scaled body size
dataSub$bodysizeScaled <- with(dataSub, bodysizeCentered/500)

# center and standardized phvage
dataSub$phvageCentered <- with(dataSub, phvage - mean(phvage, na.rm=TRUE))
dataSub$phvageStandard <- with(dataSub, phvageCentered/sd(phvage, na.rm=TRUE))

save(dataSub, file=file.path(path,"dataSub.Rdata"))

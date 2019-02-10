# Impute single values for phvage
# Get model summaries and plots

library(mgcv)
library(lme4)

library(dplyr)
library(ggplot2)
library(reshape2)

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

# path for paper
paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/plots_report/all")

# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/joR/code"))
dataPrepPath <- dirname(dirname(getwd()))

# load functions for obtaining predictions
source(file.path(dataPrepPath,"predict_functions_2prod_centered.R"))

# prep/load data
source(file.path(dataPrepPath,"data_prep.R"))

# indicator variable for missing, and indices of missingness
dataSub$phvageMiss <- is.na(dataSub$phvage)
phvageNA <- which(dataSub$phvageMiss)

length(phvageNA)

# histogram of subjects with and without phvage
dataSubSubj <- dataSub %>% group_by(ptno) %>%
  summarize(phvageMiss = first(phvageMiss),
            birthday = first(birthday),
            sex = first(sex)
           )

dataSubSubj$phvageObs <- sapply(dataSubSubj$phvageMiss,
  function(x){ifelse(x,"Missing phv age", "Observed phv age")})

ggplot(aes(birthday), data=dataSubSubj)+
  geom_histogram()+
  theme_bw(20)+
  facet_grid(phvageObs~sex)+
  labs(x="Birthdate", y="Number of subjects")
ggsave(file.path(paperPath,"bday_phvageNA.png"))

dev.new(width=6, height=5)
ggplot(aes(x=birthday, y=phvage), data=dataSub)+
  geom_point()+
  geom_smooth(method="gam")+
  facet_wrap(~sex)+
  theme_bw(20)+
  labs(x="Birthdate", y="Peak height velocity age")
ggsave(file.path(paperPath,"phvage_bday.png"))

# no random effects
fitPhv <- lm(phvage ~ birthday*male, data=dataSub) 

# try with random intercept for pedigree
# note: did not converge with random slope for pedigree
fitPhvRandPed <- lmer(phvage ~ birthday*male + (1|pedno), data=dataSub)

cbind(lm=AIC(fitPhv), lmer=AIC(fitPhvRandPed))
cbind(fixedEffect=fixef(fitPhvRandPed)[1], randomVariance=as.numeric(VarCorr(fitPhvRandPed)$pedno^2))

summary(fitPhv)

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
plot(y=varPred, x=bd)

# mean of predicted phvage
mu <- predict(fitPhv, newdata=dataSub[phvageNA,])

# impute phvage as draws from posterior
dataSub$phvage[phvageNA] <- rnorm(n=length(phvageNA), mean=mu, sd=sqrt(varPred))

# results seem reasonable
dataSubj <- dataSub %>% group_by(ptno) %>%
  summarize(birthday = first(birthday),
            phvage = first(phvage),
            sex = first(sex),
            phvageMiss = first(phvageMiss)
           )

dev.new(width=7.25, height=5)           
qplot(x=birthday, y=phvage, data=dataSubj)+
  geom_point(aes(color=phvageMiss))+
  geom_smooth(method="gam")+
  facet_wrap(~sex)+
  theme_bw(20)+
  labs(y="Peak height velocity age", x="Birthdate")+
  scale_color_manual("", labels=c("Observed","Imputed"), values=c("black","red"))
ggsave(file.path(paperPath, "phvage_postImpute.png"))

dataSub$phvageCentered <- with(dataSub, phvage - mean(phvage))

# adding phvage to model m1cent
m1centPhvageImpute <- gamm(joR ~ 
        te(skelage, birthday, k=c(5,5), bs="cr")+
        te(skelage, birthday, k=c(5,5), by=male, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=bodysizeCentered, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=male*bodysizeCentered, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=phvageCentered, bs="cr")+
        te(skelage, birthday, k=c(5,5), by=male*phvageCentered, bs="cr"),
        family=Gamma(link=log),
        random=list(pedno=~1, ptno=~1),
        data=dataSub)
        
save(m1centPhvageImpute, file="m1centPhvageImpute.Rdata")

load("m1centPhvageImpute.Rdata")

summary(m1centPhvageImpute$gam)
AIC(m1centPhvageImpute$lme)
plot(m1centPhvageImpute$gam, scheme=1)
concurvity(m1centPhvageImpute$gam)
gam.check(m1centPhvageImpute$gam)
plot(m1centPhvageImpute$lme)

# m1cent
dataSub$joRhat <- exp(predict(m1centPhvageImpute$lme))

dataMelt <- melt(dataSub[,c("skelage","ptno","birthday","joR","joRhat","sex")],
  measure.vars=c("joR","joRhat"),
  value.name="BSI"
  )
levels(dataMelt$variable) <- c("Observed", "Predicted")

ggplot(aes(x=skelage, y=BSI, group=ptno, color=birthday),data=dataMelt)+
  geom_line()+
  theme_bw(18)+
  facet_grid(sex ~ variable)+
  labs(y="BSI", x="Skeletal age")+
  scale_color_continuous("Birthdate")
ggsave(file.path(paperPath, "m1centPhvageImpute.png"))

png(file.path(paperPath,'m1centPhvageImputeResid.png'))
plot(m1centPhvageImpute$lme)
dev.off()

# plots of change over birthday ---------------------------------------------------   

alpha=0.05

# percent change
pc <- plotPerChange(m1centPhvageImpute, initialBirthday=1930, endBirthday=2000, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1centPhvageImputepc.png"))

# m1centPhvageImpute -------------------------------------------------
vc <- VarCorr(m1centPhvageImpute$lme)
sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])
sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)

out <- plotTrend(fit0 = m1centPhvageImpute,
                 alpha = alpha,
                 birthdayPred = seq(1930, 2000, 0.5),
                 sigmaOutcome = sigmaOutcome,

                 integrateRandomEffects = "MC",
                 Bmc = 1000)

outd <- plotDerivative(fit0 = m1centPhvageImpute,
                       alpha = alpha,
                       birthdayPred = seq(1930, 2000, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = 1000)

save(out, outd, file = "out_outd_BSI_single_impute.Rdata")

# mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
    scale_x_continuous(breaks=seq(1930,2000,10),
                       labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centPhvageImputeMaleMean.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
    scale_x_continuous(breaks=seq(1930,2000,10),
                       labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centPhvageImputeFemaleMean.png"))

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
    labs(y=expression(paste(partialdiff," BSI / ", 
                            partialdiff," birthdate", sep="")), 
        x="Birthdate")+
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centPhvageImputeMaleDer.png"))
  
slice <- sliceFun(outd, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    geom_hline(yintercept=0, color="red")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y=expression(paste(partialdiff," BSI / ", 
                            partialdiff," birthdate", sep="")), 
        x="Birthdate")+
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centPhvageImputeFemaleDer.png"))

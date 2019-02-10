# Evaluate the results of the ttar models,
# which are fit in 'ttar_modelFits_cent.R'.
# Produces summary plots of trends, derivatives, and percent change.

library(mgcv)
library(reshape2)
library(ggplot2)

# set computer-specific paths (note: 'path' points to the folder
# containing the data, and is used in 'dataPred.R')
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
paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/plots_report/all")
# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/ttar/code"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))

# load functions for obtaining predictions
source(file.path(dataPrepPath,"predict_functions_2prod_centered.R"))

# 1-alpha (conditional) credible intervals
alpha=0.05

# transparency for 3d plots
alphaPlot=0.7

# evaluate models -------------------------------------------------------------
load("m1ttarCent.Rdata")  
load("m2ttarCent.Rdata")  
load("m3ttarCent.Rdata")  
load("m4ttarCent.Rdata")  
load("m5ttarCent.Rdata")  

aic <- data.frame(model=c("m1ttarCent",
                          "m2ttarCent",
                          "m3ttarCent",
                          "m4ttarCent",
                          "m5ttarCent"),
                  aic=c(AIC(m1ttarCent$lme),
                        AIC(m2ttarCent$lme),
                        AIC(m3ttarCent$lme),
                        AIC(m4ttarCent$lme),
                        AIC(m5ttarCent$lme))
  )
aic

Rsq <- sapply(list(m1ttarCent$gam,
                   m2ttarCent$gam,
                   m3ttarCent$gam,
                   m4ttarCent$gam,
                   m5ttarCent$gam),
              function(x){summary(x)$r.sq})
round(Rsq, 3)

summary(m1ttarCent$gam)
plot(m1ttarCent$gam, scheme=1)
concurvity(m1ttarCent$gam)
gam.check(m1ttarCent$gam) 

# variance components
vc <- VarCorr(m1ttarCent$lme)
nvc <- nrow(vc)
vc[(nvc-4):nvc,]

vc <- VarCorr(m2ttarCent$lme)
nvc <- nrow(vc)
vc[(nvc-5):nvc,]

vc <- VarCorr(m3ttarCent$lme)
nvc <- nrow(vc)
vc[(nvc-2):nvc,]

vc <- VarCorr(m4ttarCent$lme)
nvc <- nrow(vc)
vc[(nvc-3):nvc,]

vc <- VarCorr(m5ttarCent$lme)
nvc <- nrow(vc)
vc[(nvc-5):nvc,]

# plot predicted quantities ------------------------------------

# m1cent
dataSub$ttarHat <- exp(predict(m1ttarCent$lme))

dataSubphv <- dataSub[which(!is.na(dataSub$phvage)),]
max(dataSubphv$birthday)

dataMelt <- melt(dataSub[,c("skelage","ptno","birthday","ttar","ttarHat","sex")],
                 measure.vars=c("ttar","ttarHat"),
                 value.name="ttar")
levels(dataMelt$variable) <- c("Observed", "Predicted")

dev.new(height = 7, width = 8)
ggplot(aes(x=skelage, y=ttar, group=ptno, color=birthday),data=dataMelt)+
  geom_line()+
  theme_bw(20)+
  facet_grid(sex ~ variable)+
  labs(y="Total area", x="Skeletal area")+
  scale_color_continuous("Birthdate")
ggsave(file.path(paperPath, "m1ttarCent.png"))
 
png(file.path(paperPath, "m1ttarCentResid.png"))
plot(m1ttarCent$lme)
dev.off()

qqnorm(residuals(m1ttarCent$gam, type="deviance"))

# plots of change over birthday ---------------------------------------------------   

# percent change
pc <- plotPerChange(m1ttarCent,
                    initialBirthday=1930,
                    endBirthday=1995,
                    ttar=TRUE)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1ttarPc.png"))


# averages over all ages
minAge <- 8
mean(pc$simFemale[which(pc$skelagePred >= minAge),])*100
quantile(pc$simFemale[which(pc$skelagePred >= minAge),], c(alpha/2, 1-alpha/2))*100

mean(pc$simMale[which(pc$skelagePred >= minAge),])*100
quantile(pc$simMale[which(pc$skelagePred >= minAge),], c(alpha/2, 1-alpha/2))*100

# separate plots for males and females
dev.new(height=3.5, width=3.9)
ggplot(aes(x=skelage, y=mean*100), data=subset(pc$zAll, comparison=="Male"))+
    geom_line()+
    geom_line(aes(y=l*100),linetype="dashed")+
    geom_line(aes(y=u*100),linetype="dashed")+
    theme_bw(16)+
    labs(y=expression(paste("Percent change in TA",sep="")), 
         x="Skeletal age")+
    geom_hline(yintercept=0, color="red")+
    scale_x_continuous(breaks=seq(min(pc$skelage), max(pc$skelage), 2))+
    scale_y_continuous(breaks=seq(-40,40,5), lim=c(-27,7))
ggsave(file.path(paperPath, "m1ttarPc_male.png"))

dev.new(height=3.5, width=3.9)
ggplot(aes(x=skelage, y=mean*100), data=subset(pc$zAll, comparison=="Female"))+
    geom_line()+
    geom_line(aes(y=l*100),linetype="dashed")+
    geom_line(aes(y=u*100),linetype="dashed")+
    theme_bw(16)+
    labs(y=expression(paste("Percent change in TA",sep="")),
         x="Skeletal age")+
    geom_hline(yintercept=0, color="red")+
    scale_x_continuous(breaks=seq(min(pc$skelage), max(pc$skelage), 2))+
    scale_y_continuous(breaks=seq(-40,40,5), lim=c(-27,7))
ggsave(file.path(paperPath, "m1ttarPc_female.png"))

# 2d plots --------------------------------------------------------------------

# m1cent -------------------------------------------------
vc <- VarCorr(m1ttarCent$lme)
sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])
sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)

out <- plotTrend(fit0 = m1ttarCent,
                 alpha = alpha,
                 birthdayPred = seq(1930, 2000, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = 1000)

outd <- plotDerivative(fit0 = m1ttarCent,
                       alpha = alpha,
                       birthdayPred = seq(1930, 2000, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = 1000)

save(out, outd, file = "out_outd_TA.Rdata")

# mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="TA", x="Birthdate")+
    scale_x_continuous(breaks=seq(1930,2000,10),
                       labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,70))
ggsave(file.path(paperPath, "m1ttarMaleMean.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
  geom_line()+
  geom_line(aes(y=l),linetype="dashed")+
  geom_line(aes(y=u),linetype="dashed")+
  theme_bw(17)+
  facet_grid(~skelage)+
  labs(y="TA", x="Birthdate")+
  scale_x_continuous(breaks=seq(1930,2000,10),
                     labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
  scale_y_continuous(lim=c(0,70))
ggsave(file.path(paperPath, "m1ttarFemaleMean.png"))

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
  labs(y=expression(paste(partialdiff," TA / ", 
                          partialdiff," birthdate", sep="")), 
       x="Birthdate")+
  scale_x_continuous(breaks=seq(1930,2000,10),
                     labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
  scale_y_continuous(lim=c(-0.4, 0.2))
ggsave(file.path(paperPath, "m1ttarMaleDer.png"))
  
slice <- sliceFun(outd, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
  geom_line()+
  geom_line(aes(y=l),linetype="dashed")+
  geom_line(aes(y=u),linetype="dashed")+
  geom_hline(yintercept=0, color="red")+
  theme_bw(17)+
  facet_grid(~skelage)+
  labs(y=expression(paste(partialdiff," TA / ", 
                          partialdiff," birthdate", sep="")), 
       x="Birthdate")+
  scale_x_continuous(breaks=seq(1930,2000,10),
                     labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
  scale_y_continuous(lim=c(-0.4, 0.2))
ggsave(file.path(paperPath, "m1ttarFemaleDer.png"))

# Evaluate the results of the joR models,
# which are fit in 'joR_modelFits_cent.R'.
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
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/joR/code"))

dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# load functions for obtaining predictions
source(file.path(dataPrepPath,"predict_functions_2prod_centered.R"))

# prep/load data
source(file.path(dataPrepPath,"data_prep.R"))


# 1-alpha (conditional) credible intervals
alpha <- 0.05

Bmc <- 1000

# evaluate models -------------------------------------------------------------
load("m1cent.Rdata")
load("m1centPhvage.Rdata")
load("m1centNoNAPhvage.Rdata")
load("m1centPre1995.Rdata")
load("m1centPre1990.Rdata")
load("m1centPedPre1950.Rdata")
load("m1centNoOutliers.Rdata")
load("m1centRestrictBMI.Rdata")

aic <- data.frame(model=c("m1cent","m1centPhvage","m1centNoNAPhvage",
  "m1centPre1995", "m1centPre1990", "m1centPedPre1950", "m1centNoOutliers"),
  aic=c(AIC(m1cent$lme), AIC(m1centPhvage$lme), AIC(m1centNoNAPhvage$lme),
    AIC(m1centPre1995$lme), AIC(m1centPre1990$lme), AIC(m1centPedPre1950$lme),
    AIC(m1centNoOutliers$lme))
  )
aic

summary(m1cent$gam)

plot(m1cent$gam, scheme=1)
plot(m1centPhvage$gam, scheme=1)
plot(m1centNoNAPhvage$gam, scheme=1)
plot(m1centPre1995$gam, scheme=1)

concurvity(m1cent$gam)
concurvity(m1centPhvage$gam)
concurvity(m1centNoNAPhvage$gam)
concurvity(m1centPre1995$gam)

gam.check(m1cent$gam) 
gam.check(m1centPhvage$gam) 
gam.check(m1centNoNAPhvage$gam) 
gam.check(m1centPre1995$gam) 

plot(m1cent$lme)
plot(m1centPhvage$lme)
plot(m1centNoNAPhvage$lme)
plot(m1centPre1995$lme)

# variance components
vc <- VarCorr(m1cent$lme)
nvc <- nrow(vc)
vc[(nvc-4):nvc,]

dataSub$expectedVals <- exp(predict(m1cent$lme))
dataSub$errors <- with(dataSub, joR - expectedVals)

residFit <- lm(log(errors^2) ~ log(expectedVals), data = dataSub)
coef(residFit)

# plot predicted quantities ------------------------------------

# m1cent
dataSub$joRhat <- exp(predict(m1cent$lme))

dataSubphv <- dataSub[which(!is.na(dataSub$phvage)),]
max(dataSubphv$birthday)
# 1995.562

dataMelt <- melt(dataSub[,c("skelage","ptno","birthday","joR","joRhat","sex")],
                 measure.vars=c("joR","joRhat"),
                 value.name="joR")
  
levels(dataMelt$variable) <- c("Observed", "Predicted")

dev.new(height = 7, width = 8)
ggplot(aes(x=skelage, y=joR, group=ptno, color=birthday),data=dataMelt)+
  geom_line()+
  theme_bw(20)+
  facet_grid(sex ~ variable)+
  labs(y="BSI", x="Skeletal age")+
  scale_x_continuous(breaks = seq(6, 18, 2)) +
  scale_color_continuous("Birthdate")
ggsave(file.path(paperPath, "m1cent_joR.png"))
 
png(file.path(paperPath,'m1centResid_joR.png'))
plot(m1cent$lme)
dev.off()

qqnorm(residuals(m1cent$gam, type="deviance"))

# plots of change over birthday ---------------------------------------------------   

# percent change
pc <- plotPerChange(fit0 = m1cent, initialBirthday=1930, endBirthday=2000, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1pc_joR.png"))

# averages over all ages
minAge <- 8
mean(pc$simFemale[which(pc$skelagePred >= minAge),])*100
# [1] -13.80297
quantile(pc$simFemale[which(pc$skelagePred >= minAge),], c(alpha/2, 1-alpha/2))*100
#       2.5%      97.5% 
# -19.824719  -7.254802 

mean(pc$simMale[which(pc$skelagePred >= minAge),])*100
# [1] -18.09804
quantile(pc$simMale[which(pc$skelagePred >= minAge),], c(alpha/2, 1-alpha/2))*100
#      2.5%     97.5% 
# -23.81684 -11.87419 

# separate plots for males and females
dev.new(height=3.5, width=3.9)
ggplot(aes(x=skelage, y=mean*100), data=subset(pc$zAll, comparison=="Male"))+
    geom_line()+
    geom_line(aes(y=l*100),linetype="dashed")+
    geom_line(aes(y=u*100),linetype="dashed")+
    theme_bw(16)+
    labs(
      # title= paste("Percent change ", initialBirthday, " to ", endBirthday, ": ",(1-alpha)*100, "% credible intervals", sep=""),
      y=expression(paste("Percent change in BSI",sep="")), x="Skeletal age")+
    geom_hline(yintercept=0, color="red")+
    scale_x_continuous(breaks=seq(min(pc$skelage), max(pc$skelage), 2))+
    scale_y_continuous(breaks=seq(-40,40,5), lim=c(-27,7))
ggsave(file.path(paperPath, "m1pc_joR_male.png"))

dev.new(height=3.5, width=3.9)
ggplot(aes(x=skelage, y=mean*100), data=subset(pc$zAll, comparison=="Female"))+
    geom_line()+
    geom_line(aes(y=l*100),linetype="dashed")+
    geom_line(aes(y=u*100),linetype="dashed")+
    theme_bw(16)+
    labs(
      # title= paste("Percent change ", initialBirthday, " to ", endBirthday, ": ",(1-alpha)*100, "% credible intervals", sep=""),
      y=expression(paste("Percent change in BSI",sep="")), x="Skeletal age")+
    geom_hline(yintercept=0, color="red")+
    scale_x_continuous(breaks=seq(min(pc$skelage), max(pc$skelage), 2))+
    scale_y_continuous(breaks=seq(-40,40,5), lim=c(-27,7))
ggsave(file.path(paperPath, "m1pc_joR_female.png"))

# birthdays only go through 1995 for subjects with phv age observations
pc <- plotPerChange(m1centPhvage, initialBirthday=1930, endBirthday=1995, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1phvpc_joR.png"))

# birthdays only go through 1995 for subjects with phv age observations
pc <- plotPerChange(m1centNoNAPhvage, initialBirthday=1930, endBirthday=1995, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1NoNAphvpc_joR.png"))

pc <- plotPerChange(m1centPre1995, initialBirthday=1930, endBirthday=1995, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1centPre1995pc_joR.png"))

pc <- plotPerChange(m1centPre1990, initialBirthday=1930, endBirthday=1990, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1centPre1990pc_joR.png"))

pc <- plotPerChange(m1centPedPre1950, initialBirthday=1930, endBirthday=2000, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1centPedPre1950pc_joR.png"))

pc <- plotPerChange(m1centNoOutliers, initialBirthday=1930, endBirthday=2000, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1centNoOutlierspc_joR.png"))

pc <- plotPerChange(m1centRestrictBMI, initialBirthday=1930, endBirthday=2000, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1centRestrictBMIpc_joR.png"))

outList <- list()
outdList <- list()
# 2d plots --------------------------------------------------------------------

# m1cent -------------------------------------------------
vc <- VarCorr(m1cent$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

limOutcome <- c(sigmaPedOutcome, sigmaPtnoOutcome) * 10

out <- plotTrend(fit0 = m1cent,
                 alpha = alpha,
                 birthdayPred = seq(1930, 2000, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = 1000)

outd <- plotDerivative(fit0 = m1cent,
                       alpha = alpha,
                       birthdayPred = seq(1930, 2000, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = 1000)

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
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,150,20))

ggsave(file.path(paperPath, "m1maleMean_joR.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))

ggsave(file.path(paperPath, "m1femaleMean_joR.png"))

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
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))

ggsave(file.path(paperPath, "m1maleDer_joR.png"))
  
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
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))

ggsave(file.path(paperPath, "m1femaleDer_joR.png"))

# m1centPhvage -- birthdays only go through 1995
out <- plotTrend(m1centPhvage, alpha=alpha, birthdayPred=seq(1930,1995, 0.5))
outd <- plotDerivative(m1centPhvage, alpha=alpha, birthdayPred=seq(1930,1995, 0.5))

# 2d plots: mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1phvMaleMean_joR.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1phvFemaleMean_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1phvMaleDer_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1phvFemaleDer_joR.png"))


# m1centNoNAPhvage -- birthdays only go through 1995
out <- plotTrend(m1centNoNAPhvage, alpha=alpha, birthdayPred=seq(1930, 1995, 0.5))
outd <- plotDerivative(m1centNoNAPhvage, alpha=alpha, birthdayPred=seq(1930, 1995, 0.5))

# 2d plots: mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1NoNAphvMaleMean_joR.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1NoNAphvFemaleMean_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1NoNAphvMaleDer_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1NoNAphvFemaleDer_joR.png"))


# m1centPre1995
out <- plotTrend(m1centPre1995, alpha=alpha, birthdayPred=seq(1930, 1995, 0.5))
outd <- plotDerivative(m1centPre1995, birthdayPred=seq(1930, 1995, 0.5), alpha=alpha)

# 2d plots: mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centPre1995MaleMean_joR.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centPre1995FemaleMean_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centPre1995MaleDer_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1995,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centPre1995FemaleDer_joR.png"))


# m1centPre1990
out <- plotTrend(m1centPre1990, alpha=alpha, birthdayPred=seq(1930, 1990, 0.5))
outd <- plotDerivative(m1centPre1990, alpha=alpha, birthdayPred=seq(1930, 1990, 0.5))

# 2d plots: mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1990,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centPre1990MaleMean_joR.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1990,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centPre1990FemaleMean_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1990,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centPre1990MaleDer_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,1990,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990"))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centPre1990FemaleDer_joR.png"))

# m1centPedPre1950
out <- plotTrend(m1centPedPre1950, alpha=alpha, birthdayPred=seq(1930, 2000, 0.5))
outd <- plotDerivative(m1centPedPre1950, alpha=alpha, birthdayPred=seq(1930, 2000, 0.5))

# 2d plots: mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centPedPre1950MaleMean_joR.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centPedPre1950FemaleMean_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centPedPre1950MaleDer_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centPedPre1950FemaleDer_joR.png"))


# m1centNoOutliers
out <- plotTrend(m1centNoOutliers, alpha=alpha, birthdayPred=seq(1930, 2000, 0.5))
outd <- plotDerivative(m1centNoOutliers, alpha=alpha, birthdayPred=seq(1930, 2000, 0.5))

# 2d plots: mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))

# m1centRestrictBMI
out <- plotTrend(m1centRestrictBMI, alpha=alpha, birthdayPred=seq(1930, 2000, 0.5))
outd <- plotDerivative(m1centRestrictBMI, alpha=alpha, birthdayPred=seq(1930, 2000, 0.5))

# 2d plots: mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centRestrictBMIMaleMean_joR.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="BSI", x="Birthdate")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,145), breaks=seq(0,140,20))
ggsave(file.path(paperPath, "m1centRestrictBMIFemaleMean_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centRestrictBMIMaleDer_joR.png"))

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
        x="Birthdate")+      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.9, 0.9), breaks=seq(-0.8,0.8,0.4))
ggsave(file.path(paperPath, "m1centRestrictBMIFemaleDer_joR.png"))


# 3d plots
# m1cent -------------------------------------------------
c <- VarCorr(m1cent$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

limOutcome <- c(sigmaPedOutcome, sigmaPtnoOutcome) * 10

out <- plotTrend(m1cent,
                 alpha=alpha,
                 skelagePred=seq(8,18, 0.25), 
                 birthdayPred=seq(1930, 2000, 1),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = 1000)

# outd <- plotDerivative(m1cent, alpha=alpha, skelagePred=seq(8,18, 0.5), birthdayPred=seq(1930, 2000, 0.5))

# males
i=1
# colors <- heat.colors(1000)[cut(out$z[[i]]$mean, quantile(out$z[[i]]$mean, seq(0,1,.001)))]
cuts <- cut(t(out$z[[i]]$mean), quantile(t(out$z[[i]]$mean), seq(0,1,.001)))

# manually fix one cell, which is showing as NA -- I'm not sure why
breakPoints <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", levels(cuts)) ),
      upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(cuts)) ))
index <- which(t(out$z[[i]]$mean)[2,1] > breakPoints[,1] & t(out$z[[i]]$mean)[2,1] <= breakPoints[,2])
cuts[71] <- levels(cuts)[index]

colors <- heat.colors(1000)[cuts]
# zMinMax <- c(min(out$z[[i]]$l), max(out$z[[i]]$u))
zMinMax <- c(30, 140)

png(file.path(paperPath,"m1MaleMean_joR_3d.png"))
persp(y=out$skelagePred, x=out$birthdayPred, z=t(out$z[[i]]$mean), col=colors, 
  ylab="Skeletal age", xlab="Birthdate", zlab="BSI", zlim=zMinMax,
  shade=NA, ticktype="detailed", theta=25,
  cex.lab=1.25, cex.axis=1.25)
dev.off()

# females
i=2
# colors <- heat.colors(1000)[cut(out$z[[i]]$mean, quantile(out$z[[i]]$mean, seq(0,1,.001)))]
cuts <- cut(t(out$z[[i]]$mean), quantile(t(out$z[[i]]$mean), seq(0,1,.001)))

# manually fix one cell, which is showing as NA -- I'm not sure why
breakPoints <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", levels(cuts)) ),
      upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(cuts)) ))
index <- which(t(out$z[[i]]$mean)[2,1] > breakPoints[,1] & t(out$z[[i]]$mean)[2,1] <= breakPoints[,2])
cuts[71] <- levels(cuts)[index]
colors <- heat.colors(1000)[cuts]
# zMinMax <- c(min(out$z[[i]]$l), max(out$z[[i]]$u))
zMinMax <- c(30, 140)

png(file.path(paperPath,"m1FemaleMean_joR_3d.png"))
persp(y=out$skelagePred, x=out$birthdayPred, z=t(out$z[[i]]$mean), col=colors, 
  ylab="Skeletal age", xlab="Birthdate", zlab="BSI", zlim=zMinMax,
    ticktype="detailed", theta=25,
    cex.lab=1.25, cex.axis=1.25)
dev.off()

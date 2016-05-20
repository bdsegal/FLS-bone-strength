# Evaluate the results of the cti models,
# which are fit in 'cti_modelFits_cent.R'.
# Produces summary plots of trends, derivatives, and percent change.

# note: 'tci' and 'cti' refer to the same variable
# it is labeled as tci in the dataset

library(mgcv)

library(reshape2)
library(ggplot2)
library(rgl)

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
paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/cti/paper")
# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/cti/code"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))

# load functions for obtaining predictions
source(file.path(dataPrepPath,"predict_functions_2prod_centered_normal.R"))

# 1-alpha (conditional) credible intervals
alpha=0.01

# transparency for 3d plots
alphaPlot=0.7

# evaluate models -------------------------------------------------------------
load("m1tciCent.Rdata")  
load("m2tciCent.Rdata")  
load("m3tciCent.Rdata")  
load("m4tciCent.Rdata")  
load("m5tciCent.Rdata")  

aic <- data.frame(model=c("m1tciCent","m2tciCent","m3tciCent",
  "m4tciCent", "m5tciCent"),
  aic=c(AIC(m1tciCent$lme), AIC(m2tciCent$lme), AIC(m3tciCent$lme),
    AIC(m4tciCent$lme), AIC(m5tciCent$lme))
  )
aic

Rsq <- sapply(list(m1tciCent$gam, m2tciCent$gam, m3tciCent$gam, m4tciCent$gam, m5tciCent$gam), function(x){summary(x)$r.sq})
signif(Rsq, 3)

summary(m1tciCent$gam)

plot(m1tciCent$gam, scheme=1)

concurvity(m1tciCent$gam)

gam.check(m1tciCent$gam) 

plot(m1tciCent$lme)

# variance components
vc <- VarCorr(m1tciCent$lme)
nvc <- nrow(vc)
vc[(nvc-4):nvc,]

vc <- VarCorr(m2tciCent$lme)
nvc <- nrow(vc)
vc[(nvc-5):nvc,]

vc <- VarCorr(m3tciCent$lme)
nvc <- nrow(vc)
vc[(nvc-2):nvc,]

vc <- VarCorr(m4tciCent$lme)
nvc <- nrow(vc)
vc[(nvc-3):nvc,]

vc <- VarCorr(m5tciCent$lme)
nvc <- nrow(vc)
vc[(nvc-5):nvc,]


# average slope of skeleton age
outSkel <- plotSkelageDeriv(m1tciCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,] 0.008141 0.01164 0.01514

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,] 0.008606 0.01211 0.0156

outSkel <- plotSkelageDeriv(m2tciCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,] 0.008289 0.01111 0.01393

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,] 0.009357 0.01194 0.01453

outSkel <- plotSkelageDeriv(m3tciCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,] 0.008164 0.01166 0.01515

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,]  0.008582 0.01209 0.01558

outSkel <- plotSkelageDeriv(m4tciCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,]  0.008333 0.01113 0.01391


i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,]  0.009321 0.01193 0.01455

outSkel <- plotSkelageDeriv(m5tciCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,] 0.007798 0.01122 0.01463

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,] 0.008242 0.01168 0.01512


# plot predicted quantities ------------------------------------

# m1cent
dataSub$tciHat <- predict(m1tciCent$lme)

dataMelt <- melt(dataSub[,c("skelage","ptno","birthday","mc2tci","tciHat","sex")],
  measure.vars=c("mc2tci","tciHat"),
  # variable.name="joR"
  value.name="cti"
  )
levels(dataMelt$variable) <- c("Observed", "Predicted")

ggplot(aes(x=skelage, y=cti, group=ptno, color=birthday),data=dataMelt)+
  geom_line()+
  theme_bw(20)+
  facet_grid(sex ~ variable)+
  labs(y="Cortical thickness index", x="Skeletal age")+
  scale_color_continuous("Birthday")
ggsave(file.path(paperPath, "m1tciCent.png"))
 
png(file.path(paperPath, "m1tciCentResid.png"))
plot(m1tciCent$lme)
dev.off()

qqnorm(residuals(m1ttarCent$gam, type="deviance"))

# plots of change over birthday ---------------------------------------------------   

# percent change
pc <- plotPerChange(m1tciCent, initialBirthday=1930, endBirthday=2000, alpha=alpha)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1tciPc.png"))


# 2d plots --------------------------------------------------------------------

# m1cent -------------------------------------------------
out <- plotTrend(m1tciCent, alpha=alpha, birthdayPred=seq(1930,2000, 0.5))
outd <- plotDerivative(m1tciCent, alpha=alpha, birthdayPred=seq(1930,2000, 0.5))

# mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="Cortical thickness index", x="Birthday")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990", ""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,0.65))
ggsave(file.path(paperPath, "m1tciMaleMean.png"))

slice <- sliceFun(out, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="Cortical thickness index", x="Birthday")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(0,0.65))
ggsave(file.path(paperPath, "m1tciFemaleMean.png"))

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
    labs(y=expression(paste(partialdiff," CTI / ", 
          partialdiff," birthday", sep="")), 
        x="Birthday")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.0025, 0.0075))
ggsave(file.path(paperPath, "m1tciMaleDer.png"))
  
slice <- sliceFun(outd, i=2)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    geom_hline(yintercept=0, color="red")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y=expression(paste(partialdiff," CTI / ", 
          partialdiff," birthday", sep="")), 
        x="Birthday")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.0025, 0.0075))
ggsave(file.path(paperPath, "m1tciFemaleDer.png"))


# 3d plots

out <- plotTrend(m1tciCent, alpha=alpha, birthdayPred=seq(1930,1995, 0.5))
outd <- plotDerivative(m1tciCent, alpha=alpha, birthdayPred=seq(1930,1995, 0.5))

i=1
colors <- heat.colors(1000)[cut(out$z[[i]]$mean, quantile(out$z[[i]]$mean, seq(0,1,.001)))]

zMinMax <- c(min(out$z[[i]]$l), max(out$z[[i]]$u))
zMinMax <- c(0, 0.65)

open3d()
persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$mean, col=colors, 
main=paste(names(out$z)[i],": ", (1-out$alpha)*100,"% credible Intervals", sep=""), xlab="skeleton age", ylab="birthday", zlab="jo/R",
zlim=zMinMax, add=FALSE)
persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$l, add=TRUE, col="grey", alpha=alphaPlot)
persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$u, add=TRUE, col="grey", alpha=alphaPlot)

# partial derivative with respect to birthday

i=1
colors <- heat.colors(1000)[cut(outd$z[[i]]$mean, quantile(outd$z[[i]]$mean, seq(0,1,.001)))]
zMinMax <- c(min(outd$z[[i]]$l), max(outd$z[[i]]$u))
zMinMax <- c(-0.0025, 0.0075)

open3d()
persp3d(x=outd$skelagePred, y=outd$birthdayPred, z=outd$z[[i]]$mean, col=colors, 
main=paste(names(outd$z)[i],": ", (1-out$alpha)*100,"% credible Intervals", sep=""), xlab="skeleton age", ylab="birthday", zlab="jo/R",
zlim=zMinMax, add=FALSE)
persp3d(x=outd$skelagePred, y=outd$birthdayPred, z=outd$z[[i]]$l, add=TRUE, col="grey", alpha=alphaPlot)
persp3d(x=outd$skelagePred, y=outd$birthdayPred, z=outd$z[[i]]$u, add=TRUE, col="grey", alpha=alphaPlot)
persp3d(x=outd$skelagePred, y=outd$birthdayPred, z=array(0,dim=dim(outd$z[[i]]$mean)), add=TRUE, col="grey", alpha=1)
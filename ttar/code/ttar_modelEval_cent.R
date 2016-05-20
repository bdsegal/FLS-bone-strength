# Evaluate the results of the ttar models,
# which are fit in 'ttar_modelFits_cent.R'.
# Produces summary plots of trends, derivatives, and percent change.

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
paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/ttar/paper")
# directory with code
setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis/ttar/code"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))

# load functions for obtaining predictions
source(file.path(dataPrepPath,"predict_functions_2prod_centered.R"))

# 1-alpha (conditional) credible intervals
alpha=0.01

# transparency for 3d plots
alphaPlot=0.7

# evaluate models -------------------------------------------------------------
load("m1ttarCent.Rdata")  
load("m2ttarCent.Rdata")  
load("m3ttarCent.Rdata")  
load("m4ttarCent.Rdata")  
load("m5ttarCent.Rdata")  

aic <- data.frame(model=c("m1ttarCent","m2ttarCent","m3ttarCent",
  "m4ttarCent", "m5ttarCent"),
  aic=c(AIC(m1ttarCent$lme), AIC(m2ttarCent$lme), AIC(m3ttarCent$lme),
    AIC(m4ttarCent$lme), AIC(m5ttarCent$lme))
  )
aic

Rsq <- sapply(list(m1ttarCent$gam, m2ttarCent$gam, m3ttarCent$gam, m4ttarCent$gam, m5ttarCent$gam), function(x){summary(x)$r.sq})
signif(Rsq, 3)
# [1] 0.702 0.702 0.702 0.702 0.703


summary(m1ttarCent$gam)

plot(m1ttarCent$gam, scheme=1)

concurvity(m1ttarCent$gam)

gam.check(m1ttarCent$gam) 

plot(m1ttarCent$lme)

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


# average slope of skeleton age
outSkel <- plotSkelageDeriv(m1ttarCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,] 0.06151 0.06651 0.07152

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,] 0.04209 0.04621 0.05035


outSkel <- plotSkelageDeriv(m2ttarCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,] 0.06225 0.06705 0.07184

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,] 0.04178 0.04591 0.05007

outSkel <- plotSkelageDeriv(m3ttarCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,] 0.06153 0.06651 0.07152

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,]  0.04206 0.04621 0.05035

outSkel <- plotSkelageDeriv(m4ttarCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,]  0.06225 0.06702 0.07178

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,]  0.04174 0.0459 0.05005

outSkel <- plotSkelageDeriv(m5ttarCent)
i=1
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
      # lower   mean  upper
# [1,] 0.06163 0.06643 0.07123

i=2
skelageSlope <- cbind(mean(outSkel$z[[i]]$l),
  mean(outSkel$z[[i]]$mean),
  mean(outSkel$z[[i]]$u))
colnames(skelageSlope) <- c("lower","mean","upper")
signif(skelageSlope,4)
       # lower    mean   upper
# [1,]  0.04088 0.04586 0.05083


# plot predicted quantities ------------------------------------

# m1cent
dataSub$ttarHat <- exp(predict(m1ttarCent$lme))

dataSubphv <- dataSub[which(!is.na(dataSub$phvage)),]
max(dataSubphv$birthday)
# 1995.562

dataMelt <- melt(dataSub[,c("skelage","ptno","birthday","ttar","ttarHat","sex")],
  measure.vars=c("ttar","ttarHat"),
  # variable.name="joR"
  value.name="ttar"
  )
levels(dataMelt$variable) <- c("Observed", "Predicted")

ggplot(aes(x=skelage, y=ttar, group=ptno, color=birthday),data=dataMelt)+
  geom_line()+
  theme_bw(20)+
  facet_grid(sex ~ variable)+
  labs(y="Total area", x="Skeletal area")+
  scale_color_continuous("Birthday")
ggsave(file.path(paperPath, "m1ttarCent.png"))
 
png(file.path(paperPath, "m1ttarCentResid.png"))
plot(m1ttarCent$lme)
dev.off()

qqnorm(residuals(m1ttarCent$gam, type="deviance"))

# plots of change over birthday ---------------------------------------------------   

# percent change
pc <- plotPerChange(m1ttarCent, initialBirthday=1930, endBirthday=1995, alpha=alpha, ttar=TRUE)
dev.new(height=3.5, width=7)
pc$ggAll
ggsave(file.path(paperPath, "m1ttarPc.png"))


# 2d plots --------------------------------------------------------------------


# m1cent -------------------------------------------------
out <- plotTrend(m1ttarCent, alpha=alpha, birthdayPred=seq(1930,2000, 0.5))
outd <- plotDerivative(m1ttarCent, alpha=alpha, birthdayPred=seq(1930,2000, 0.5))

# mean trend
slice <- sliceFun(out, i=1)
dev.new(width=10, height=3.5)
ggplot(aes(x=birthday, y=mean), data=slice)+
    geom_line()+
    geom_line(aes(y=l),linetype="dashed")+
    geom_line(aes(y=u),linetype="dashed")+
    theme_bw(17)+
    facet_grid(~skelage)+
    labs(y="Total area", x="Birthday")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
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
    labs(y="Total area", x="Birthday")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
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
    labs(y=expression(paste(partialdiff," Total area / ", 
          partialdiff," birthday", sep="")), 
        x="Birthday")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
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
    labs(y=expression(paste(partialdiff," Total area / ", 
          partialdiff," birthday", sep="")), 
        x="Birthday")+
      # title=paste("Male ", (1-out$alpha)*100, "% credible intervals", sep="")
    scale_x_continuous(breaks=seq(1930,2000,10),
      labels=c("1930", "", "1950", "", "1970", "", "1990",""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))+
    scale_y_continuous(lim=c(-0.4, 0.2))
ggsave(file.path(paperPath, "m1ttarFemaleDer.png"))


# 3d plots

out <- plotTrend(m1ttarCent, alpha=alpha)
outd <- plotDerivative(m1ttarCent, alpha=alpha)

i=1
colors <- heat.colors(1000)[cut(out$z[[i]]$mean, quantile(out$z[[i]]$mean, seq(0,1,.001)))]

zMinMax <- c(min(out$z[[i]]$l), max(out$z[[i]]$u))

open3d()
persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$mean, col=colors, 
main=paste(names(out$z)[i],": ", (1-out$alpha)*100,"% credible Intervals", sep=""), xlab="skeleton age", ylab="birthday", zlab="jo/R",
zlim=zMinMax, add=FALSE)
persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$l, add=TRUE, col="grey", alpha=alphaPlot)
persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$u, add=TRUE, col="grey", alpha=alphaPlot)

# partial derivative with respect to birthday

i=2
colors <- heat.colors(1000)[cut(outd$z[[i]]$mean, quantile(outd$z[[i]]$mean, seq(0,1,.001)))]
zMinMax <- c(min(outd$z[[i]]$l), max(outd$z[[i]]$u))

open3d()
persp3d(x=outd$agePred, y=outd$yearsPred, z=outd$z[[i]]$mean, col=colors, 
main=paste(names(outd$z)[i],": ", (1-out$alpha)*100,"% credible Intervals", sep=""), xlab="skeleton age", ylab="birthday", zlab="jo/R",
zlim=zMinMax, add=FALSE)
persp3d(x=outd$agePred, y=outd$yearsPred, z=outd$z[[i]]$l, add=TRUE, col="grey", alpha=alphaPlot)
persp3d(x=outd$agePred, y=outd$yearsPred, z=outd$z[[i]]$u, add=TRUE, col="grey", alpha=alphaPlot)
persp3d(x=outd$agePred, y=outd$yearsPred, z=array(0,dim=dim(outd$z[[i]]$mean)), add=TRUE, col="grey", alpha=1)
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

outList <- list()
outdList <- list()
# 2d plots --------------------------------------------------------------------

# m1cent --------------------------------------------------
vc <- VarCorr(m1cent$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

limOutcome <- c(sigmaPedOutcome, sigmaPtnoOutcome) * 10

system.time(
out <- plotTrend(fit0 = m1cent,
                 alpha = alpha,
                 birthdayPred = seq(1930, 2000, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc))

outList[[1]] <- out
save(outList, file = "outList.Rdata")

system.time(
outd <- plotDerivative(fit0 = m1cent,
                       alpha = alpha,
                       birthdayPred = seq(1930, 2000, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = Bmc))

outdList[[1]] <- outd
save(outdList, file = "outdList.Rdata")

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

# 3d plots ------------------------------------------------

system.time(
out <- plotTrend(m1cent,
                 alpha=alpha,
                 skelagePred=seq(8,18, 0.25), 
                 birthdayPred=seq(1930, 2000, 1),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc))

save(out, file = "out_m1cent_3d.Rdata")
load("out_m1cent_3d.Rdata")

# males
i=1
cuts <- cut(t(out$z[[i]]$mean), quantile(t(out$z[[i]]$mean), seq(0,1,.001)))

# manually fix one cell, which is showing as NA -- I'm not sure why
breakPoints <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", levels(cuts)) ),
                     upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(cuts)) ))
index <- which(t(out$z[[i]]$mean)[2,1] > breakPoints[,1] & t(out$z[[i]]$mean)[2,1] <= breakPoints[,2])
cuts[71] <- levels(cuts)[index]

colors <- heat.colors(1000)[cuts]
zMinMax <- c(30, 145)

png(file.path(paperPath,"m1MaleMean_joR_3d.png"))
persp(y=out$skelagePred, x=out$birthdayPred, z=t(out$z[[i]]$mean), col=colors, 
  ylab="Skeletal age", xlab="Birthdate", zlab="BSI", zlim=zMinMax,
  shade=NA, ticktype="detailed", theta=25,
  cex.lab=1.25, cex.axis=1.25)
dev.off()

# females
i=2
cuts <- cut(t(out$z[[i]]$mean), quantile(t(out$z[[i]]$mean), seq(0,1,.001)))

breakPoints <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", levels(cuts)) ),
                     upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(cuts)) ))
index <- which(t(out$z[[i]]$mean)[2,1] > breakPoints[,1] & t(out$z[[i]]$mean)[2,1] <= breakPoints[,2])
cuts[71] <- levels(cuts)[index]
colors <- heat.colors(1000)[cuts]
zMinMax <- c(30, 145)

png(file.path(paperPath,"m1FemaleMean_joR_3d.png"))
persp(y=out$skelagePred, x=out$birthdayPred, z=t(out$z[[i]]$mean), col=colors, 
  ylab="Skeletal age", xlab="Birthdate", zlab="BSI", zlim=zMinMax,
    ticktype="detailed", theta=25,
    cex.lab=1.25, cex.axis=1.25)
dev.off()

# m1centPhvage -- birthdays only go through 1995 ----------
vc <- VarCorr(m1centPhvage$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

limOutcome <- c(sigmaPedOutcome, sigmaPtnoOutcome) * 10

out <- plotTrend(fit0 = m1centPhvage,
                 alpha = alpha,
                 birthdayPred = seq(1930, 1995, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc)

outList[[2]] <- out
save(outList, file = "outList.Rdata")

outd <- plotDerivative(fit0 = m1centPhvage,
                       alpha = alpha,
                       birthdayPred = seq(1930, 1995, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = Bmc)

outdList[[2]] <- outd
save(outdList, file = "outdList.Rdata")

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

# m1centNoNAPhvage -- birthdays only go through 1995 ------
vc <- VarCorr(m1centNoNAPhvage$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

limOutcome <- c(sigmaPedOutcome, sigmaPtnoOutcome) * 10

out <- plotTrend(fit0 = m1centNoNAPhvage,
                 alpha = alpha,
                 birthdayPred = seq(1930, 1995, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc)

outList[[3]] <- out
save(outList, file = "outList.Rdata")

outd <- plotDerivative(fit0 = m1centNoNAPhvage,
                       alpha = alpha,
                       birthdayPred = seq(1930, 1995, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = Bmc)

outdList[[3]] <- outd
save(outdList, file = "outdList.Rdata")

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


# m1centPre1995 -------------------------------------------
vc <- VarCorr(m1centPre1995$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

out <- plotTrend(fit0 = m1centPre1995,
                 alpha = alpha,
                 birthdayPred = seq(1930, 1995, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc)

outList[[4]] <- out
save(outList, file = "outList.Rdata")

outd <- plotDerivative(fit0 = m1centPre1995,
                       alpha = alpha,
                       birthdayPred = seq(1930, 1995, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = Bmc)

outdList[[4]] <- outd
save(outdList, file = "outdList.Rdata")

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


# m1centPre1990 -------------------------------------------
vc <- VarCorr(m1centPre1990$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

out <- plotTrend(fit0 = m1centPre1990,
                 alpha = alpha,
                 birthdayPred = seq(1930, 1990, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc)

outList[[5]] <- out
save(outList, file = "outList.Rdata")

outd <- plotDerivative(fit0 = m1centPre1990,
                       alpha = alpha,
                       birthdayPred = seq(1930, 1990, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = Bmc)

outdList[[5]] <- outd
save(outdList, file = "outdList.Rdata")

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

# m1centPedPre1950 ----------------------------------------
vc <- VarCorr(m1centPedPre1950$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

out <- plotTrend(fit0 = m1centPedPre1950,
                 alpha = alpha,
                 birthdayPred = seq(1930, 2000, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc)

outList[[6]] <- out
save(outList, file = "outList.Rdata")

outd <- plotDerivative(fit0 = m1centPedPre1950,
                       alpha = alpha,
                       birthdayPred = seq(1930, 2000, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = Bmc)

outdList[[6]] <- outd
save(outdList, file = "outdList.Rdata")

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


# m1centNoOutliers ----------------------------------------
vc <- VarCorr(m1centNoOutliers$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

out <- plotTrend(fit0 = m1centNoOutliers,
                 alpha = alpha,
                 birthdayPred = seq(1930, 2000, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc)

outList[[7]] <- out
save(outList, file = "outList.Rdata")

outd <- plotDerivative(fit0 = m1centNoOutliers,
                       alpha = alpha,
                       birthdayPred = seq(1930, 2000, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = Bmc)

outdList[[7]] <- outd
save(outdList, file = "outdList.Rdata")

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

# m1centRestrictBMI ---------------------------------------
vc <- VarCorr(m1centRestrictBMI$lme)

sigmaPedOutcome <- as.numeric(vc[which(grepl("pedno", rownames(vc))) + 1, 2])
sigmaPtnoOutcome <- as.numeric(vc[which(grepl("ptno", rownames(vc))) + 1, 2])

sigmaOutcome <- diag(c(sigmaPedOutcome, sigmaPtnoOutcome)^2)
logdetOutcome <- sum(log(eigen(sigmaOutcome, symmetric = TRUE, only.values = TRUE)$values))

out <- plotTrend(fit0 = m1centRestrictBMI,
                 alpha = alpha,
                 birthdayPred = seq(1930, 2000, 0.5),
                 sigmaOutcome = sigmaOutcome,
                 logdetOutcome = logdetOutcome,
                 limOutcome = limOutcome,
                 integrateRandomEffects = "MC",
                 Bmc = Bmc)

outList[[8]] <- out
save(outList, file = "outList.Rdata")

outd <- plotDerivative(fit0 = m1centRestrictBMI,
                       alpha = alpha,
                       birthdayPred = seq(1930, 2000, 0.5),
                       sigmaOutcome = sigmaOutcome,
                       logdetOutcome = logdetOutcome,
                       limOutcome = limOutcome,
                       integrateRandomEffects = "MC",
                       Bmc = Bmc)

outdList[[8]] <- outd
save(outdList, file = "outdList.Rdata")

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

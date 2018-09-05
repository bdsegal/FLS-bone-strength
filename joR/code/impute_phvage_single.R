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
#             fixedEffect randomVariance
# (Intercept)    28.02805       0.930322

summary(fitPhv)
# Call:
# lm(formula = phvage ~ birthday * male, data = dataSub)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -4.4463 -0.5950  0.0838  0.7043  2.8583 

# Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   21.099251   2.002271  10.538  < 2e-16 ***
# birthday      -0.004888   0.001023  -4.779 1.79e-06 ***
# male           9.652101   2.846108   3.391 0.000699 ***
# birthday:male -0.003851   0.001454  -2.649 0.008091 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.063 on 8891 degrees of freedom
  # (1385 observations deleted due to missingness)
# Multiple R-squared:  0.5005,    Adjusted R-squared:  0.5003 
# F-statistic:  2969 on 3 and 8891 DF,  p-value: < 2.2e-16

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

# fit0 <- m1centPhvageImpute

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
  # variable.name="joR"
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

# qqnorm(residuals(m1cent$gam, type="deviance"))


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


# # 3d plots
# alphaPlot=0.7

# i=2
# colors <- heat.colors(1000)[cut(out$z[[i]]$mean, quantile(out$z[[i]]$mean, seq(0,1,.001)))]
# zMinMax <- c(min(out$z[[i]]$l), max(out$z[[i]]$u))

# open3d()
# persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$mean, col=colors, 
# main=paste(names(out$z)[i],": ", (1-out$alpha)*100,"% credible Intervals", sep=""), xlab="skeleton age", ylab="birthday", zlab="jo/R",
# zlim=zMinMax, add=FALSE)
# persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$l, add=TRUE, col="grey", alpha=alphaPlot)
# persp3d(x=out$skelagePred, y=out$birthdayPred, z=out$z[[i]]$u, add=TRUE, col="grey", alpha=alphaPlot)

# # partial derivative with respect to birthday

# i=2
# colors <- heat.colors(1000)[cut(outd$z[[i]]$mean, quantile(outd$z[[i]]$mean, seq(0,1,.001)))]
# zMinMax <- c(min(outd$z[[i]]$l), max(outd$z[[i]]$u))

# open3d()
# persp3d(x=outd$skelagePred, y=outd$birthdayPred, z=outd$z[[i]]$mean, col=colors, 
# main=paste(names(outd$z)[i],": ", (1-out$alpha)*100,"% credible Intervals", sep=""), xlab="skeleton age", ylab="birthday", zlab="jo/R",
# zlim=zMinMax, add=FALSE)
# persp3d(x=outd$skelagePred, y=outd$birthdayPred, z=outd$z[[i]]$l, add=TRUE, col="grey", alpha=alphaPlot)
# persp3d(x=outd$skelagePred, y=outd$birthdayPred, z=outd$z[[i]]$u, add=TRUE, col="grey", alpha=alphaPlot)
# persp3d(x=outd$skelagePred, y=outd$birthdayPred, z=array(0,dim=dim(outd$z[[i]]$mean)), add=TRUE, col="grey", alpha=1)

# Get basic statistics on number of subjects and observations in sample

library(mgcv)
library(dplyr)
library(ggplot2)
library(pheatmap)

# set computer-specific paths (note: 'path' points to the folder
# containing the data, and is used in 'data_pred.R')
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

setwd(file.path(computer,"Dropbox/Research/Bones/final_analysis"))
dataPrepPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis")

paperPath <- file.path(computer,"Dropbox/Research/Bones/final_analysis/plots_report/all")

# prep data
source(file.path(dataPrepPath,"data_prep.R"))

sampleStats <- rbind(
    # subjects
    c(length(unique(data$ptno)),
    length(unique(dataSub$ptno))),
    # male
    c(length(unique(data$ptno[which(data$sex == "Male")])), 
      length(unique(dataSub$ptno[which(dataSub$sex == "Male")]))),
    #female
    c(length(unique(data$ptno[which(data$sex == "Female")])),
      length(unique(dataSub$ptno[which(dataSub$sex == "Female")]))),
    #families
    c(length(unique(data$pedno)),
      length(unique(dataSub$pedno))),
    # total observations
    c(nrow(data),
     nrow(dataSub))
  )

colnames(sampleStats) <- c("all", "analyzed")
rownames(sampleStats) <- c("subjects", "male", "female", 
  "families", "total observations")
sampleStats
#                      all analyzed
# subjects            1204     1093
# male                 619      559
# female               585      534
# families             209      177
# total observations 16937    10280


obsPerSubj <- summarise(group_by(dataSub, ptno),
    count = n(),
    sex = sex[1]
)

ggplot(aes(x = count), data=obsPerSubj) +
    geom_histogram(breaks=0:max(obsPerSubj$count))+
    theme_bw(22)+
    facet_wrap(~sex)+
    labs(y = "Count", x = "Number of observations")

# number of subjects by pedigree
dataSub$birthyear <- trunc(dataSub$birthday)

uniq <- unique(dataSub[,c("pedno","ptno")])
uniq <- uniq[ order(uniq[,1],uniq[,2]),]

# add birthyear to unique subjects
uniq$birthyear <- NA
for (i in 1:nrow(uniq)){
    ptnoTemp <- uniq$ptno[i]
    uniq$birthyear[i] <- dataSub$birthyear[which(dataSub$ptno==ptnoTemp)][1]
}

# change NA to 999 to make sure all subjects are counted
# uniq[which(is.na(uniq$pedno)),"pedno"] <- 999

# make matrix for plotting
yrRange <- with(uniq,c(min(birthyear),max(birthyear)))
yrSeq <- seq(yrRange[1], yrRange[2], 1)

pedByBirth <- matrix(nrow=length(unique(uniq$pedno)),
                    ncol=length(yrSeq))
dimnames(pedByBirth) <- list(Pedigree=unique(uniq$pedno),
    Birthyear=yrSeq)
rownames(pedByBirth) <- unique(uniq$pedno)
colnames(pedByBirth) <- yrSeq

# fill matrix
for (i in 1:nrow(pedByBirth)){
    for (j in 1:ncol(pedByBirth)){
            
        temp <- uniq[which(uniq$pedno == unique(uniq$pedno)[i] &
                    uniq$birthyear == yrSeq[j]),]

        pedByBirth[i,j] <- nrow(temp)
    }
}

# plotting labels and colors
colLabels <- rep("",ncol(pedByBirth))
colLabels[seq(1,71,5)] <- as.character(seq(1930,2000,5))
colors <- colorRampPalette(c("white","navy blue"))

dimnames(pedByBirth) <- list

png(file.path(paperPath,"pedByBirth.png"))
pheatmap(pedByBirth,
    cluster_cols=FALSE, cluster_rows=FALSE,
    show_rownames=FALSE,labels_col=colLabels,
    color=colors(6),fontsize_col=12,
    annotation_names_row=TRUE,
    cex=1.2)
    # main="Counts: pedigree (row) by birth year (column)")
dev.off()

range(dataSub$birthday)

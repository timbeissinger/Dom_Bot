################################################################################
### Use this script to map teosintes from Maize HapMap 2                     ###
################################################################################

### Timothy M. Beissinger
### 6-18-2014

### Load libraries
library(maps)
library(mapdata)
library(mapplots)


### Load the data
data <- read.csv("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot_Data/LR_teosinte_coord.csv",header=T,stringsAsFactors=F)


############################### All Accessions ################################
### Create pdf
pdf("../Dom_Bot_Figs/maizeAmericas.pdf")

### Make the map
map("worldHires",col="lightgrey",fill=T,ylim=c(-60,48),xlim=c(-130,-35),main="Location of Hapmap 2 accessions")
box()
map.scale(cex=0.5)
title("Location of Hapmap 2 accessions")

### Subset genotypes
parv <- data[which(data$subspecies=="parviglumis"),]
mex <- data[which(data$subspecies=="mexicana"),]
maize <- data[which(data$subspecies=="mays"),]

### Plot genotypes
points(mex$Long_Dec,mex$Latit_Dec,pch=19,col="blue",cex=0.75)
points(maize$Long_Dec,maize$Latit_Dec,pch=19,col="darkgreen",cex=0.75)
points(parv$Long_Dec,parv$Latit_Dec,pch=19,col="red",cex=0.75)

### Make legend
legend("topright","(x,y)",c("parviglumis", "mexicana","landrace"),col=c("red","blue","darkgreen"),pch=19,cex=0.75)

### End pdf
dev.off()



############################### Parviglumis Accessions##########################
### Create pdf
pdf("../Dom_Bot_Figs/parvMexico.pdf")

### Make the map
map("worldHires",col="lightgrey",fill=T,ylim=c(14,25),xlim=c(-110,-95),main="Location of Hapmap 2 accessions")
box()
map.scale(cex=0.5)
title("Location of Hapmap 2 parviglumis accessions")

### Plot parviglumis accessions
points(parv$Long_Dec,parv$Latit_Dec,pch=19,col="red",cex=2)
points(parv$Long_Dec[which(parv$Altitude <= 500)],parv$Latit_Dec[which(parv$Altitude <= 500)],col="black",pch=19,cex=2)

### Legend
legend("top","(x,y)",c("Altitude < 500", "Altitude > 500"),col=c("black","red"),pch=19,cex=1)

### End pdf
dev.off()

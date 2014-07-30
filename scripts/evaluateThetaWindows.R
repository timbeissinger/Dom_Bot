################################################################################
### Use this script to evaluate thetaStats computed over windows            ####
################################################################################

### Timothy M. Beissinger
### 7/29/2014


### Load pestPG files
BKNgenic <- read.table("../THETAS/BKN_genic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
BKNintergenic <- read.table("../THETAS/BKN_intergenic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
TILgenic <- read.table("../THETAS/TIL_genic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)

### Include only windows with many sites included
BKNgenicCovered <- BKNgenic[which(BKNgenic$nSites>0),]
BKNintergenicCovered <- BKNintergenic[which(BKNintergenic$nSites>0),]
TILgenicCovered <- TILgenic[which(TILgenic$nSites>0),]

### Simple plot
plot(BKNgenicCovered$WinCenter,BKNgenicCovered$Tajima)
plot(BKNintergenicCovered$WinCenter,BKNintergenicCovered$Tajima)
plot(TILgenicCovered$WinCenter,TILgenicCovered$Tajima)

### Mean
mean(BKNintergenicCovered$Tajima)
mean(BKNgenicCovered$Tajima)
mean(TILgenicCovered$Tajima)

### Hist
hist(BKNintergenicCovered$Tajima)
hist(BKNgenicCovered$Tajima)
hist(TILgenicCovered$Tajima)


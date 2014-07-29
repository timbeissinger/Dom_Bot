################################################################################
### Use this script to evaluate thetaStats computed over windows            ####
################################################################################

### Timothy M. Beissinger
### 7/29/2014


### Load pestPG files
BKNgenic <- read.table("../THETAS/BKN_genic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
TILgenic <- read.table("../THETAS/TIL_genic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)

### Include only windows with many sites included
BKNgenicCovered <- BKNgenic[which(BKNgenic$nSites>500),]
TILgenicCovered <- TILgenic[which(TILgenic$nSites>500),]

### Simple plot
plot(BKNgenicCovered$WinCenter,BKNgenicCovered$Tajima)
plot(TILgenicCovered$WinCenter,TILgenicCovered$Tajima)

### Mean
mean(BKNgenicCovered$Tajima)
mean(TILgenicCovered$Tajima)

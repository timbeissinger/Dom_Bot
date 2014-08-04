################################################################################
### Use this script to evaluate thetaStats computed over windows            ####
################################################################################

### Timothy M. Beissinger
### 7/29/2014


### Load pestPG files
BKNgenic <- read.table("../THETAS/BKN_genic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
BKNintergenic <- read.table("../THETAS/BKN_intergenic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
TILgenic <- read.table("../THETAS/TIL_genic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
TILintergenic <- read.table("../THETAS/TIL_intergenic_10_windows.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)

### Include only windows with many sites included
BKNgenicCovered <- BKNgenic[which(BKNgenic$nSites>100),]
BKNintergenicCovered <- BKNintergenic[which(BKNintergenic$nSites>100),]
TILgenicCovered <- TILgenic[which(TILgenic$nSites>100),]
TILintergenicCovered <- TILintergenic[which(TILintergenic$nSites>100),]

### Simple plot
par(mfrow=c(2,2))
plot(BKNgenicCovered$WinCenter,BKNgenicCovered$Tajima)
plot(BKNintergenicCovered$WinCenter,BKNintergenicCovered$Tajima)
plot(TILgenicCovered$WinCenter,TILgenicCovered$Tajima)
plot(TILintergenicCovered$WinCenter,TILintergenicCovered$Tajima)

### Mean
mean(BKNintergenicCovered$Tajima)
mean(BKNgenicCovered$Tajima)
mean(TILgenicCovered$Tajima)
mean(TILintergenicCovered$Tajima,na.rm=T)

### Hist
breaks=seq(-2.5,4.5,0.2)
hist(BKNintergenicCovered$Tajima,breaks=breaks)
hist(BKNgenicCovered$Tajima,breaks=breaks)
hist(TILgenicCovered$Tajima,breaks=breaks)
hist(TILintergenicCovered$Tajima,breaks=breaks)

### Overlapping hists

pdf("../THETAS/TajimaHists.pdf")
breaks=seq(-2.5,4.5,0.2)
### BKN
hist(BKNgenicCovered$Tajima,breaks=breaks,freq=F,col=rgb(1,0,1,0.5),ylim=c(0,0.7),main="Maize: Tajima's D on Chromosome 10",xlab="Tajima's D",ylab="Density")
par(new=T)
hist(BKNintergenicCovered$Tajima,breaks=breaks,freq=F,col=rgb(1,0.5,0,0.5),ylim=c(0,0.7),ylab="",xlab="",main="")
legend("topright","(x,y)",pch=15,col=c(rgb(1,0,1,0.5),rgb(1,0.5,0,0.5)),c(paste("Genic (mean = ", signif(mean(BKNgenicCovered$Tajima),2),")",sep=""),paste("Non-genic (mean = ", signif(mean(BKNintergenicCovered$Tajima),2),")",sep="")))



hist(TILgenicCovered$Tajima,breaks=breaks,freq=F,col=rgb(1,0,1,0.5),ylim=c(0,0.7),main="TIL: Tajima's D on Chromosome 10",xlab="Tajima's D",ylab="Density")
par(new=T)
hist(TILintergenicCovered$Tajima,breaks=breaks,freq=F,col=rgb(1,0.5,0,0.5),ylim=c(0,0.7),ylab="",xlab="",main="")
legend("topright","(x,y)",pch=15,col=c(rgb(1,0,1,0.5),rgb(1,0.5,0,0.5)),c(paste("Genic (mean = ", signif(mean(TILgenicCovered$Tajima),2),")",sep=""),paste("Non-genic (mean = ", signif(mean(TILintergenicCovered$Tajima),2),")",sep="")))

dev.off()





###############################################################################################
##### THETAS ON ONLY OVERLAPPING SITES (BETWEEN BKN AND TIL)                              #####
###############################################################################################

### Load pestPG files
BKNgenicOverlapping <- read.table("../THETAS/BKN_genic_10_windows_overlapping.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
BKNintergenicOverlapping <- read.table("../THETAS/BKN_intergenic_10_windows_overlapping.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
TILgenicOverlapping <- read.table("../THETAS/TIL_genic_10_windows_overlapping.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)
TILintergenicOverlapping <- read.table("../THETAS/TIL_intergenic_10_windows_overlapping.thetas.gz.pestPG",header=T,stringsAsFactors=F,comment.char="",skip=1)

### Include only windows with many sites included
BKNgenicOverlappingCovered <- BKNgenicOverlapping[which(BKNgenicOverlapping$nSites>100),]
BKNintergenicOverlappingCovered <- BKNintergenicOverlapping[which(BKNintergenicOverlapping$nSites>100),]
TILgenicOverlappingCovered <- TILgenicOverlapping[which(TILgenicOverlapping$nSites>100),]
TILintergenicOverlappingCovered <- TILintergenicOverlapping[which(TILintergenicOverlapping$nSites>100),]

### Simple plot
par(mfrow=c(2,2))
plot(BKNgenicOverlappingCovered$WinCenter,BKNgenicOverlappingCovered$Tajima)
plot(BKNintergenicOverlappingCovered$WinCenter,BKNintergenicOverlappingCovered$Tajima)
plot(TILgenicOverlappingCovered$WinCenter,TILgenicOverlappingCovered$Tajima)
plot(TILintergenicOverlappingCovered$WinCenter,TILintergenicOverlappingCovered$Tajima)

### Mean
mean(BKNintergenicOverlappingCovered$Tajima)
mean(BKNgenicOverlappingCovered$Tajima)
mean(TILgenicOverlappingCovered$Tajima)
mean(TILintergenicOverlappingCovered$Tajima,na.rm=T)

### Hist
breaks=seq(-2.5,4.5,0.2)
hist(BKNintergenicOverlappingCovered$Tajima,breaks=breaks)
hist(BKNgenicOverlappingCovered$Tajima,breaks=breaks)
hist(TILgenicOverlappingCovered$Tajima,breaks=breaks)
hist(TILintergenicOverlappingCovered$Tajima,breaks=breaks)

### Overlapping hists

pdf("../THETAS/TajimaHists_sitesInBoth.pdf")
breaks=seq(-2.5,4.5,0.2)
### BKN
hist(BKNgenicOverlappingCovered$Tajima,breaks=breaks,freq=F,col=rgb(1,0,1,0.5),ylim=c(0,0.7),main="Maize: Tajima's D on Chromosome 10",xlab="Tajima's D",ylab="Density")
par(new=T)
hist(BKNintergenicOverlappingCovered$Tajima,breaks=breaks,freq=F,col=rgb(1,0.5,0,0.5),ylim=c(0,0.7),ylab="",xlab="",main="")
legend("topright","(x,y)",pch=15,col=c(rgb(1,0,1,0.5),rgb(1,0.5,0,0.5)),c(paste("Genic (mean = ", signif(mean(BKNgenicOverlappingCovered$Tajima),2),")",sep=""),paste("Non-genic (mean = ", signif(mean(BKNintergenicOverlappingCovered$Tajima),2),")",sep="")))



hist(TILgenicOverlappingCovered$Tajima,breaks=breaks,freq=F,col=rgb(1,0,1,0.5),ylim=c(0,0.7),main="TIL: Tajima's D on Chromosome 10",xlab="Tajima's D",ylab="Density")
par(new=T)
hist(TILintergenicOverlappingCovered$Tajima,breaks=breaks,freq=F,col=rgb(1,0.5,0,0.5),ylim=c(0,0.7),ylab="",xlab="",main="")
legend("topright","(x,y)",pch=15,col=c(rgb(1,0,1,0.5),rgb(1,0.5,0,0.5)),c(paste("Genic (mean = ", signif(mean(TILgenicOverlappingCovered$Tajima),2),")",sep=""),paste("Non-genic (mean = ", signif(mean(TILintergenicOverlappingCovered$Tajima),2),")",sep="")))

dev.off()

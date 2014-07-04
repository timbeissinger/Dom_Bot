##############################################################################
### This script can be used to plot the sfs output from Dadi for the BKN   ###
### population, the TIL population, and the 2dsfs output from both.        ###
### Designed to run locally, on my imac.                                   ###
##############################################################################

### Timothy M. Beissinger
### 7-3-2014

### Set working directory
setwd("~/Documents/DomesticationBottleneck/Dom_Bot_Git/SFS/")

### Load package
library(gplots)

### Read in the spectra
sfsBKN <- exp(scan("BKN.sfs"))
sfsTIL <- exp(scan("TIL.sfs"))
sfs2d <- read.table( "2dsfs.BKN.TIL.sfs",header=F,stringsAsFactors=F)

### Plot the sfs by population
pdf("SFS_Plots_individual_pops.pdf")
barplot(sfsBKN[2:{length(sfsBKN)-1}],main="Landrace SFS")
barplot(sfsTIL[2:{length(sfsTIL)-1}],main="Teosinte SFS")
dev.off()

### Make a 2d sfs matrix
sfs2d.mat <- exp(as.matrix(sfs2d))
rownames(sfs2d.mat) <- 0:{nrow(sfs2d.mat)-1}
colnames(sfs2d.mat) <- 0:{ncol(sfs2d.mat)-1}

### 2d sfs heatmap
# define color breaks
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue", "red", "yellow"))(n = 299)

col_breaks = c(seq(0,0.00005,length=100),  # for blue
  seq(0.00005,0.0001,length=100),              # for red
  seq(0.0001,0.001,length=100))              # for yellow

heatmap.2(sfs2d.mat[2:{nrow(sfs2d.mat)-1},2:{ncol(sfs2d.mat)-1}],Rowv=F,Colv="NA",dendrogram="none",density.info="none",trace="none",col=my_palette, breaks=col_breaks,revC=T,keysize=1)

pdf("BKN-TIL-2d-sfs.pdf")
heatmap.2(sfs2d.mat[2:{nrow(sfs2d.mat)-1},2:{ncol(sfs2d.mat)-1}],Rowv=F,Colv="NA",dendrogram="none",density.info="none",trace="none",col=my_palette, breaks=col_breaks,revC=T,keysize=0.5,lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(2, 8, 2 ), main="2d sfs for teo & landraces (0-10MB on chr 10)",xlab="No. Teosinte chromosomes with derived allele",ylab="No. Landrace chromosomes with derived allele")
dev.off()

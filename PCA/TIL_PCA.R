################################################################################
### Compute PCA between BKN lines using v2 data downloaded from panzea       ###
################################################################################

### Timothy M. Beissinger
### 7/28/2014

library(adegenet)

### Load the data (trimmed to every thousandth line with the shell commands below)
# head -n 1 HapMapv2_data_TIL.txt > HapMapv2_data_TIL_trimmed.txt
# awk '!(NR%1000)' HapMapv2_data_TIL.txt >> HapMapv2_data_TIL_trimmed.txt
markers <- read.table("~/Documents/DomesticationBottleneck/Dom_Bot_Data/HapMapv2_data_TIL_trimmed.txt", header=T, stringsAsFactors=F,comment.char="")


### Remove indels (only want SNPs)
SNPs <- which(markers[,15]!="+" & markers[,15]!="-")
markers <- markers[SNPs,]

### Make a matrix of just markers
data <- as.matrix(markers[,12:ncol(markers)])

data[which(data=="N")] <- NA
data[which(data=="X")] <- NA
data[which(data=="A")] <- "A/A"
data[which(data=="C")] <- "C/C"
data[which(data=="G")] <- "G/G"
data[which(data=="T")] <- "T/T"
data[which(data=="K")] <- "G/T"
data[which(data=="M")] <- "A/C"
data[which(data=="R")] <- "A/G"
data[which(data=="S")] <- "C/G"
data[which(data=="W")] <- "A/T"
data[which(data=="Y")] <- "C/T"


### Determine which markers are NOT polymporphic
isPolymorphic <- function(data){
    if(length(table(data)) > 1) return(TRUE)
    else return(FALSE)
}

poly <- apply(data,1,isPolymorphic)


### Pick 10,000 SNPs to use for structure ### DO NOT USE SNPS THAT ARE COMPLETELY MISSING!!!
index <- sample(which(poly==T),replace=F,size=10000)
index <- index[order(index)]
data <- data[index,]

### Format ames matrix into a data frame with rownames equal to taxa,
### no additional data
tTIL <- data.frame(t(data),stringsAsFactors=F)

### Make genind object
tilGenind <- df2genind(tTIL,ploidy=2,sep="/",ind.names=rownames(tTIL))

### Make data matrix
tilMat <- scaleGen(tilGenind,missing="mean")

#### Perform PCA
tilPCA <- dudi.pca(tilMat,cent=F,scale=F,scannf=F,nf=34)

### Plot PCA PCA
pdf("tilPCA.pdf",height=8,width=12)
par(mfrow=c(1,2))
barplot(tilPCA$eig[1:34]/sum(tilPCA$eig),main="TIL PCA eigenvalues",col=heat.colors(34))
plot(tilPCA$li[,1],tilPCA$li[,2],pch=16,cex=0.5,main="TIL PC1 vs PC2",xlab="PC 1",ylab="PC 2")
text(labels=substr(rownames(tTIL),1,6),tilPCA$li[,1],tilPCA$li[,2],cex=0.5,pos=1)
dev.off()

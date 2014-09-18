################################################################################
### Use this script to make a plot of diversity surrounding Syn, non-syn     ###
### substitutions between tripsicum and maize.                               ###
################################################################################

### 9/16/2014

### Load effects estimates
effects <- read.table("../SNPs/TvMeffects.txt",header=T,stringsAsFactors=F,sep="\t",comment.char="",na.strings="-",skip=8)
levels(as.factor(effects$Consequence))

effects <- effects[,c(1,2,3,4,5,6,7)]

### Make a set of syn, non variants
syn <- effects[which(effects$Consequence == "synonymous_variant"),]
mis <- effects[which(effects$Consequence == "missense_variant"),]

### Remove ambiguous subs (positions in both syn, mis)
amb <- intersect(syn$Location,mis$Location)
syn <- syn[-which(syn$Location %in% amb),]
mis <- mis[-which(mis$Location %in% amb),]

### Remove duplicate positions syn (multiple transcripts)
syn0 <- syn[NULL,]
levels <- levels(as.factor(syn$Location))
nlevels <- length(levels(as.factor(syn$Location)))

for(i in 1:nlevels){
print(i)
uniqueRow <- which(syn$Location==levels[i])[1]
syn0[nrow(syn0)+1,] <- syn[uniqueRow,]
}

### Remove duplicate positions mis (multiple transcripts)
mis0 <- mis[NULL,]
levels <- levels(as.factor(mis$Location))
nlevels <- length(levels(as.factor(mis$Location)))

for(i in 1:nlevels){
print(i)
uniqueRow <- which(mis$Location==levels[i])[1]
mis0[nrow(mis0)+1,] <- mis[uniqueRow,]
}

### Put syn0 and mis0 in order
options(scipen=10)
syn0$chr <- as.numeric(unlist(strsplit(syn0$Location,split=":"))[seq(1,2*nrow(syn0),2)])
syn0$pos <- as.numeric(unlist(strsplit(syn0$Location,split=":"))[seq(2,2*nrow(syn0),2)])

mis0$chr <- as.numeric(unlist(strsplit(mis0$Location,split=":"))[seq(1,2*nrow(mis0),2)])
mis0$pos <- as.numeric(unlist(strsplit(mis0$Location,split=":"))[seq(2,2*nrow(mis0),2)])

syn0 <- syn0[order(syn0$chr,syn0$pos),]
mis0 <- mis0[order(mis0$chr,mis0$pos),]


### Load genetic map
map <- read.table("../SNPs/NAM_phasedImputed_1cM_AllZeaGBSv2.3_allChrs/NAM_phasedImputed_1cM_AllZeaGBSv2.3_allChrs.hmp.txt",header=T,stringsAsFactors=F,sep="\t",comment.char="")
map <- map[,1:5]
ensemblUp <- map[,c(3,4,4)]
#write.table(file="../SNPs/ensemblUp.txt",ensemblUp,quote=F,col.names=F,row.names=F) # upload this file to ensembl to convert to maize v3
ensemblDown <- read.table("../SNPs/ensemblDown.gff",header=F,stringsAsFactors=F,sep="\t")

rem <- which(abs(as.numeric(ensemblUp[,2])-as.numeric(ensemblDown[,4])) > 2000000) #identify positions with massive shifts

map$posV3 <- ensemblDown[,4]

map <- map[-rem,] #remove positions with massive shifts

### CHECKPOINT ###
save.image("plotDiversity_TvM.RData")

### Interpolate genetic position for every syn0 SNP
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
syn0$cm <- NA

for(i in 1:nrow(syn0)){
print(i)
lowerIndex <- which(map$chrom == syn0$chr[i] & map$posV3 <= syn0$pos[i]) #find index of map anchors smaller than observed position
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==syn0$chr[i])][1]-1,na.rm=T) #take corresponding genetic position

higherIndex <- which(map$chrom == syn0$chr[i] & map$posV3 >= syn0$pos[i]) #find index of map anchors larger than observed position
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[syn0$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==syn0$chr[i])][length(which(map$chrom==syn0$chr[i]))]+1,na.rm=T) #take corresponding genetic position

scale <- {syn0$pos[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

syn0$cm[i] <- newGen
}


### Interpolate genetic position for every mis0 SNP
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
mis0$cm <- NA

for(i in 1:nrow(mis0)){
print(i)
lowerIndex <- which(map$chrom == mis0$chr[i] & map$posV3 <= mis0$pos[i]) #find index of map anchors smaller than observed position
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==mis0$chr[i])][1]-1,na.rm=T) #take corresponding genetic position

higherIndex <- which(map$chrom == mis0$chr[i] & map$posV3 >= mis0$pos[i]) #find index of map anchors larger than observed position
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[mis0$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==mis0$chr[i])][length(which(map$chrom==mis0$chr[i]))]+1,na.rm=T) #take corresponding genetic position

scale <- {mis0$pos[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

mis0$cm[i] <- newGen
}


### CHECKPOINT ###
save.image("plotDiversity_TvM.RData")

### Load diversity data (from angsd)
diversity <- read.table("../../WholeGenome/OUTS/BKN_WholeGenome_windows.thetas.gz.pestPG",comment.char="",skip=1,stringsAsFactors=F,header=T)
div0 <- diversity[,c(2,3,5,14)]

### Interpolate genetic position for every angsd position with diversity info.
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
divCm <- rep(NA,nrow(div0))
rows <- nrow(div0)

for(i in 1:nrow(div0)){
cat( 100*i/rows, "% done", "\r")
lowerIndex <- which(map$chrom == div0$Chr[i] & map$posV3 <= div0$WinCenter[i]) #find index of map anchors smaller than observed position
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==div0$Chr[i])][1]-1,na.rm=T) #take corresponding genetic position

higherIndex <- which(map$chrom == div0$Chr[i] & map$posV3 > div0$WinCenter[i]) #find index of map anchors larger than observed position
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[div0$Chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==div0$Chr[i])][length(which(map$chrom==div0$Chr[i]))]+1,na.rm=T) #take corresponding genetic position

scale <- {div0$WinCenter[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

divCm[i] <- newGen
}

div0$cm <- divCm

### CHECKPOINT ###
save.image("plotDiversity_TvM.RData")

### Determine p.theta up- and down- stream of each sub in mis0 
divMisRight <- matrix(NA,nrow=nrow(mis0),ncol=100)
colnames(divMisRight) <- seq(0.01,1,0.01)
divMisLeft <-  matrix(NA,nrow=nrow(mis0),ncol=100)

for(i in 1:nrow(mis0)){
print(i)
pos <- mis0$cm[i]
chr <- mis0$chr[i]
divTemp <- div0[which(div0$Chr==chr & abs(div0$cm-pos)<= 1 ),]

for(j in 1:100){
maxRight <- j/100
minRight <- j/100-0.01
divMisRight[i,j] <- mean(divTemp$tP[which({divTemp$cm-pos}<=maxRight & {divTemp$cm-pos}>=minRight)])

maxLeft <- -j/100
minLeft <- -j/100+0.01
divMisLeft[i,j] <- mean(divTemp$tP[which({divTemp$cm-pos}>=maxLeft & {divTemp$cm-pos}<=minLeft)])
}
}

### Make plot (right only, for now)
x <- rep(seq(0.01,1,0.01),nrow(mis0))
y <- as.vector(t(divMisRight))
plot(x,y)
lines(spline(x,y),col="red",lwd=3)


### Determine p.theta up- and down- stream of each sub in syn0 
divSynRight <- matrix(NA,nrow=nrow(syn0),ncol=100)
colnames(divSynRight) <- seq(0.01,1,0.01)
divSynLeft <-  matrix(NA,nrow=nrow(syn0),ncol=100)

for(i in 1:nrow(syn0)){
print(i)
pos <- syn0$cm[i]
chr <- syn0$chr[i]
divTemp <- div0[which(div0$Chr==chr & abs(div0$cm-pos)<= 1 ),]

for(j in 1:100){
maxRight <- j/100
minRight <- j/100-0.01
divSynRight[i,j] <- mean(divTemp$tP[which({divTemp$cm-pos}<=maxRight & {divTemp$cm-pos}>=minRight)])

maxLeft <- -j/100
minLeft <- -j/100+0.01
divSynLeft[i,j] <- mean(divTemp$tP[which({divTemp$cm-pos}>=maxLeft & {divTemp$cm-pos}<=minLeft)])
}
}

### Make plot (right only, for now)
x <- rep(seq(0.01,1,0.01),nrow(syn0))
y <- as.vector(t(divSynRight))
plot(x,y)
lines(spline(x,y),col="red",lwd=3)


### CHECKPOINT ###
save.image("plotDiversity_TvM.RData")

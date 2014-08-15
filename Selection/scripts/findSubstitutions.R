##############################################################################
### Use this script to identify substitutions between trip, teo, and maize ###
### from the SNPs called by angsd.                                         ###
##############################################################################

### Timothy M. Beissinger
### 8-15-2014

### Set working directory
getwd()
setwd("SNPs")

### Load SNP calls
rawSNPs <- read.table("angsd_snps.geno",header=F,stringsAsFactors=F,sep="\t",na.strings="NN")
SNPs <- rawSNPs[,1:41] #remove extra column

### Load and format header
head <- read.table("../INS/snpCallList.txt",stringsAsFactors=F) #load
head <- unlist(strsplit(head[,1],split="mergedBams/"))[seq(2,74,2)] #seperate directory
head <- unlist(strsplit(head,split="_merged.bam")) #remove tail
head <- c("chr","pos","Major","Minor",head) #add early columns before taxa

### Add header
names(SNPs) <- head

### Function to identify substitutions
whichAreSubs <- function(data,outCols=5,popCols=6:41,missing.limit=0.5){
    if(length(which(is.na(data[popCols]))) > (missing.limit*length(popCols)) ) return ("population NA Issue")
    if(length(which(is.na(data[outCols]))) > (missing.limit*length(outCols)) ) return ("outgroup NA Issue")

    out <- unique(unlist(strsplit(as.character(data[outCols]),split="")))
    out <- out[is.na(out)==F]
    pop <- unique(unlist(strsplit(as.character(data[popCols]),split="")))
    pop <- pop[is.na(pop)==F]
    if(length(out)>0 & length(pop)>0 & length(intersect(out,pop))==0) return ("Sub")
}

# Trip vs. Maize & Teo
TvMT_result<- apply(SNPs,1,whichAreSubs,outCols=5,popCols=6:41,missing.limit=0.2) # trip vs mays
TvMT <- SNPs[which(as.character(TvMT_result)=="Sub"),c(1:4,5,6)]
write.table(TvMT,file="TvMT.txt",row.names=F,col.names=T,quote=F,sep="\t")

# Trip vs. Maize
TvM_result<- apply(SNPs,1,whichAreSubs,outCols=5,popCols=19:41,missing.limit=0.2) # trip vs mays
TvM <- SNPs[which(as.character(TvM_result)=="Sub"),c(1:4,5,19)]
write.table(TvM,file="TvM.txt",row.names=F,col.names=T,quote=F,sep="\t")

# Trip vs. Teo
TvT_result<- apply(SNPs,1,whichAreSubs,outCols=5,popCols=6:18,missing.limit=0.2) # trip vs mays
TvT <- SNPs[which(as.character(TvT_result)=="Sub"),c(1:4,5,6)]
write.table(TvT,file="TvT.txt",row.names=F,col.names=T,quote=F,sep="\t")

# Teo vs. Maize
TevM_result <- apply(SNPs,1,whichAreSubs,outCols=6:18,popCols=19:41,missing.limit=0.2) # teo vs maize
TevM <- SNPs[which(as.character(TevM_result)=="Sub"),c(1:4,6,19)]
write.table(TevM,file="TevM.txt",row.names=F,col.names=T,quote=F,sep="\t")

# Print some summary stats
cat(c("Trip-Maize/Teo Subs = ",nrow(TvMT), "\n"))
cat(c("Trip-Maize Subs = ",nrow(TvM), "\n"))
cat(c("Trip-Teo Subs = ",nrow(TvT), "\n"))
cat(c("Teo-Maize Subs = ",nrow(TevM), "\n"))

# Save image
save.image("SubstitutionSpace.RData")

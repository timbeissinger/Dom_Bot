###############################################################################
### This script converts SFS output from angsd to the format readable by dadi #
###############################################################################

### Timothy M. Beissinger
### 7-14-2014


## ### Convert BKN.sfs
## sfsBKN <- exp(scan("../SFS/BKN.sfs"))*10000000 #we calculated 10 MB of DNA
## sfsBKN <- sfsBKN[seq(1,length(sfsBKN),2)] # treat as sample from haploid
## sfsBKNdadi <- matrix(nrow=2,ncol=1)
## sfsBKNdadi[1,1] <- paste(length(sfsBKN), "unfolded",sep=" ")
## sfsBKNdadi[2,1] <- paste(sfsBKN,collapse=" ")
## #sfsBKNdadi[3,1] <- paste(c(0,rep(1,{length(sfsBKN)-2}),0),collapse=" ")
## write.table(sfsBKNdadi,file="../SFS/sfsBKN.dadi",quote=F,row.names=F,col.names=F)


## ### Convert TIL.sfs
## sfsTIL <- exp(scan("../SFS/TIL.sfs"))*10000000 #we calculated 10 MB of DNA
## sfsTIL <- sfsTIL[seq(1,length(sfsTIL),2)] # treat as sample from haploid
## sfsTILdadi <- matrix(nrow=2,ncol=1)
## sfsTILdadi[1,1] <- paste(length(sfsTIL), "unfolded",sep=" ")
## sfsTILdadi[2,1] <- paste(sfsTIL,collapse=" ")
## #sfsTILdadi[3,1] <- paste(c(0,rep(1,{length(sfsTIL)-2}),0),collapse=" ")
## write.table(sfsTILdadi,file="../SFS/sfsTIL.dadi",quote=F,row.names=F,col.names=F)



### Convert BKN_genic10.sfs
sfsBKNgenic <- exp(scan("../SFS/BKN_genic_10.sfs"))
sfsBKNgenic <- sfsBKNgenic[seq(1,length(sfsBKNgenic),2)] # treat as sample from haploid
sfsBKNgenicDadi <- matrix(nrow=2,ncol=1)
sfsBKNgenicDadi[1,1] <- paste(length(sfsBKNgenic), "unfolded",sep=" ")
sfsBKNgenicDadi[2,1] <- paste(sfsBKNgenic,collapse=" ")
write.table(sfsBKNgenicDadi,file="../SFS/sfsBKN_genic_10.dadi",quote=F,row.names=F,col.names=F)

### Convert BKN_intergenic10.sfs
sfsBKNintergenic <- exp(scan("../SFS/BKN_intergenic_10.sfs"))
sfsBKNintergenic <- sfsBKNintergenic[seq(1,length(sfsBKNintergenic),2)] # treat as sample from haploid
sfsBKNintergenicDadi <- matrix(nrow=2,ncol=1)
sfsBKNintergenicDadi[1,1] <- paste(length(sfsBKNintergenic), "unfolded",sep=" ")
sfsBKNintergenicDadi[2,1] <- paste(sfsBKNintergenic,collapse=" ")
write.table(sfsBKNintergenicDadi,file="../SFS/sfsBKN_intergenic_10.dadi",quote=F,row.names=F,col.names=F)



### Convert TIL_genic10.sfs
sfsTILgenic <- exp(scan("../SFS/TIL_genic_10.sfs"))
sfsTILgenic <- sfsTILgenic[seq(1,length(sfsTILgenic),2)] # treat as sample from haploid
sfsTILgenicDadi <- matrix(nrow=2,ncol=1)
sfsTILgenicDadi[1,1] <- paste(length(sfsTILgenic), "unfolded",sep=" ")
sfsTILgenicDadi[2,1] <- paste(sfsTILgenic,collapse=" ")
write.table(sfsTILgenicDadi,file="../SFS/sfsTIL_genic_10.dadi",quote=F,row.names=F,col.names=F)

### Convert TIL_intergenic10.sfs
sfsTILintergenic <- exp(scan("../SFS/TIL_intergenic_10.sfs"))
sfsTILintergenic <- sfsTILintergenic[seq(1,length(sfsTILintergenic),2)] # treat as sample from haploid
sfsTILintergenicDadi <- matrix(nrow=2,ncol=1)
sfsTILintergenicDadi[1,1] <- paste(length(sfsTILintergenic), "unfolded",sep=" ")
sfsTILintergenicDadi[2,1] <- paste(sfsTILintergenic,collapse=" ")
write.table(sfsTILintergenicDadi,file="../SFS/sfsTIL_intergenic_10.dadi",quote=F,row.names=F,col.names=F)


### Convert 2d sfs. Script adapted from Jacob Crawford, https://github.com/mfumagalli/ngsToolsDev/blob/master/convert.2Dsfs.to.dadi.R

############################ Genic first ############################
# Read in 2D sfs genic
sfsgenic <- exp(read.table("../SFS/2dsfs_genic.TIL.BKN.sfs"))
sfsgenic <- sfsgenic[seq(1,nrow(sfsgenic),2),seq(1,ncol(sfsgenic),2)]

# Get sample sizes and make header
n1=nrow(sfsgenic);
n2=ncol(sfsgenic);
ns=paste(n1,n2,collapse=' ');
ns=paste(ns,"unfolded",collapse=' ');

# Convert 2D sfs to dadi array format
dadi=NULL;
for(i in 1:n1){
	dadi=c(dadi,as.numeric(sfsgenic[i,]))
}

# Write out dadi format to file with same name as 2D sfs file with .fs appeded
write.table(ns,file="../SFS/2d_genic_TIL_BKN.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE);
write.table(paste(dadi,collapse=' '),file="../SFS/2d_genic_TIL_BKN.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE);

##


############################ Intergenic next ############################
# Read in 2D sfs genic
sfsintergenic <- exp(read.table("../SFS/2dsfs_intergenic.TIL.BKN.sfs"))
sfsintergenic <- sfsintergenic[seq(1,nrow(sfsintergenic),2),seq(1,ncol(sfsintergenic),2)]

# Get sample sizes and make header
n1=nrow(sfsintergenic);
n2=ncol(sfsintergenic);
ns=paste(n1,n2,collapse=' ');
ns=paste(ns,"unfolded",collapse=' ');

# Convert 2D sfs to dadi array format
dadi=NULL;
for(i in 1:n1){
	dadi=c(dadi,as.numeric(sfsintergenic[i,]))
}

# Write out dadi format to file with same name as 2D sfs file with .fs appeded
write.table(ns,file="../SFS/2d_intergenic_TIL_BKN.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE);
write.table(paste(dadi,collapse=' '),file="../SFS/2d_intergenic_TIL_BKN.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE);

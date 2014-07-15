###############################################################################
### This script converts SFS output from angsd to the format readable by dadi #
###############################################################################

### Timothy M. Beissinger
### 7-14-2014


### Convert BKN.sfs
sfsBKN <- exp(scan("BKN.sfs"))
sfsBKNdadi <- matrix(nrow=2,ncol=1)
sfsBKNdadi[1,1] <- paste(length(sfsBKN), "unfolded",sep=" ")
sfsBKNdadi[2,1] <- paste(sfsBKN,collapse=" ")
#sfsBKNdadi[3,1] <- paste(c(0,rep(1,{length(sfsBKN)-2}),0),collapse=" ")
write.table(sfsBKNdadi,file="sfsBKN.dadi",quote=F,row.names=F,col.names=F)

### Convert TIL.sfs
sfsTIL <- exp(scan("TIL.sfs"))
sfsTILdadi <- matrix(nrow=2,ncol=1)
sfsTILdadi[1,1] <- paste(length(sfsTIL), "unfolded",sep=" ")
sfsTILdadi[2,1] <- paste(sfsTIL,collapse=" ")
#sfsTILdadi[3,1] <- paste(c(0,rep(1,{length(sfsTIL)-2}),0),collapse=" ")
write.table(sfsTILdadi,file="sfsTIL.dadi",quote=F,row.names=F,col.names=F)

### Convert 2d sfs. Script adapted from Jacob Crawford, https://github.com/mfumagalli/ngsToolsDev/blob/master/convert.2Dsfs.to.dadi.R

# Read in 2D sfs
sfs <- exp(read.table("2dsfs.BKN.TIL.sfs"))*10000000 #we calculated 10 MB of DNA

# Get sample sizes and make header
n1=nrow(sfs);
n2=ncol(sfs);
ns=paste(n1,n2,collapse=' ');
ns=paste(ns,"unfolded",collapse=' ');

# Convert 2D sfs to dadi array format
dadi=NULL;
for(i in 1:n1){
	dadi=c(dadi,as.numeric(sfs[i,]))
}

# Write out dadi format to file with same name as 2D sfs file with .fs appeded
write.table(ns,file="2d_BKN_TIL.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE);
write.table(paste(dadi,collapse=' '),file="2d_BKN_TIL.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE);

##

################################################################################
### Use this script to identify and make a list of maize version 3           ###
### genic regions.                                                           ###
################################################################################

### Timothy M. Beissinger
### 7-22-2014

### Install biomart
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

### Load package
library(biomaRt)

### Select mart
zea <- useMart("plants_mart_22")

### Select dataset: "Zea mays genes (AGPv3 (5b))"
zea <- useDataset("zmays_eg_gene", mart=zea)

### Check out filters
filters <- listFilters(zea)

### Pull maize genes
genes <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"),mart=zea)


###########################
### Modify genes matrix ###
###########################
genesMod <- genes

### Throw out genes that are not on chromosomes 1-10
genesMod$chromosome_name <- as.numeric(genesMod$chromosome_name)
genesMod <- genesMod[-which(is.na(genesMod$chromosome_name)),]

### Sort genes
genesMod <- genesMod[order(genesMod$chromosome_name, genesMod$start_position, genesMod$end_position),]

### DON"T DO THE FOLLOWING SECTION BECAUSE WE WANT GENES ONLY ###

## ### Add 5kb to gene start/stop positions
## genesMod$start_position <- genesMod$start_position-5000
## genesMod$start_position[which(genesMod$start_position<0)] <- 1 # set negative starts to 1
## genesMod$end_position <- genesMod$end_position+5000

##################################################################

### Merge genes that overlap
genesMod$merged <- "No"
i <- 1
while(i < nrow(genesMod)){
    print(i)
    if(genesMod$chromosome_name[i] == genesMod$chromosome_name[i+1] &
       genesMod$end_position[i] < genesMod$start_position[i+1]){
        i <- i+1
    }
    else if(genesMod$chromosome_name[i] == genesMod$chromosome_name[i+1] &
            genesMod$end_position[i] >= genesMod$start_position[i+1]
            ){
                genesMod$ensembl_gene_id[i] <- paste(genesMod$ensembl_gene_id[i],
                                                     genesMod$ensembl_gene_id[i+1],sep="_")
                genesMod$end_position[i] <- genesMod$end_position[i+1]
                genesMod$merged[i] <- "Yes"
                genesMod <- genesMod[-c(i+1),]
              }
    else if(genesMod$chromosome_name[i] != genesMod$chromosome_name[i+1]){
              print("Chromosome Up")
              i <- i+1
            }
}


### Matrix check
gap <- genesMod$start_position-c(0,genesMod$end_position[-nrow(genesMod)]) # test passed! Only 9 < 0.


### Create regionfile
regions <- paste(genesMod$chromosome_name, ":", genesMod$start_position, "-", genesMod$end_position,sep="")

### Create chromosome 10 regionfile
regions10 <- paste(genesMod$chromosome_name[which(genesMod$chromosome_name==10)], ":",
                   genesMod$start_position[which(genesMod$chromosome_name==10)], "-",
                   genesMod$end_position[which(genesMod$chromosome_name==10)],sep="")

### Write full table
write.table(regions,file="genicRegionFile.txt",row.names=F,col.names=F,quote=F)
write.table(regions10,file="genicRegionFile_chr10.txt",row.names=F,col.names=F,quote=F)

################################################################################
### This script compares SFS estimated with Inb = 1 to that estimated with   ###
### Inb = 0                                                                  ###
################################################################################

### Timothy M. Beissinger
### 7-16-2014

### Load each sfs
sfsInb <- exp(scan("../SFS/BKN.sfs"))
sfsNoInb <- exp(scan("../SFS/BKN_no_inb.sfs"))

### Compare
sfsInb - sfsNoInb # identical...

### Plot
pdf("../SFS/BKN.Inbreeding.Correct.pdf")
barplot(sfsInb[-c(1,length(sfsInb))],names.arg=1:39)
dev.off()

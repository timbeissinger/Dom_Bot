##############################################################################
### This script runs some basic dadi data analysis on ouput generated from ###
### angsd.                                                                 ###
##############################################################################

### Timothy M. Beissinger
### 7-14-2014

### Import dadi
import numpy
import dadi

### Load the teosinte spectrum object
fsSingle = dadi.Spectrum([0,100,20,10,1,0])
fsBKN = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/sfsBKN.dadi" )
fsTIL = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/sfsTIL.dadi" )

### Some summary statistics for single-population spectrum
thetaW_BKN = fsBKN.Watterson_theta()
pi_BKN = fsBKN.pi()
D_BKN = fsBKN.Tajima_D()
S_BKN = fsBKN.S()

thetaW_TIL = fsTIL.Watterson_theta()
pi_TIL = fsTIL.pi()
D_TIL = fsTIL.Tajima_D()
S_TIL = fsTIL.S()


### Plotting a single-population spectra
import pylab
dadi.Plotting.plot_1d_fs(fsBKN)
dadi.Plotting.plot_1d_fs(fsTIL)

#########################################################################
#########################################################################

### Load a 2d fs
fs2d = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/2d_BKN_TIL.dadi" )

### Time to plot the fs
import pylab
dadi.Plotting.plot_single_2d_sfs(fs2d,vmin=1)
pylab.show()

### Compute some summary statistics
S = fs2d.S() # Number of segregating sites -- percentage?
Fst = fs2d.Fst() # Fst




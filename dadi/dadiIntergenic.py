##############################################################################
### This script runs a reasonable dadi bottleneck model using n ouput      ###
### generated from angsd. Intergenic DNA is used.                          ###
##############################################################################


### Timothy M. Beissinger
### 7-30-2014

import numpy
from numpy import array

import dadi

### Load the single-pop spectra
fsBKNintergenic = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/sfsBKN_intergenic_10.dadi" )
fsTILintergenic = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/sfsTIL_intergenic_10.dadi" )

### Plotting a single-population spectra
import pylab
dadi.Plotting.plot_1d_fs(fsBKNintergenic)
dadi.Plotting.plot_1d_fs(fsTILintergenic)

### Load a 2d fs (BKN and TIL)
fs2dIntergenic = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/2d_intergenic_TIL_BKN.dadi" )

### Time to plot the fs
import pylab
dadi.Plotting.plot_single_2d_sfs(fs2dIntergenic,vmin=1,pop_ids=("TIL","BKN"))
pylab.show()

### Compute some summary statistics
S = fs2dIntergenic.S() # Number of segregating sites
Fst = fs2dIntergenic.Fst() # Fst

######################################
######### Better Bottleneck  #########
######################################

import realisticBottleneck

# determine sample sizes
ns = fs2dIntergenic.sample_sizes
#ns= fsBKN.sample_sizes


# set grid sizes
pts_l = [40,50,60]

# Let's work with the "custom" model in the bottleneck.py file
#func = realisticBottleneck.maizeBottleneck

#func = realisticBottleneck.maizeBottleneck
func = realisticBottleneck.maizeBottleneck_teoFixed

#params = array([.01, 3, 1, 0.6 ]) # junk values for testing
params = array([.01, 3, 0.6 ]) # junk values for testing

#upper_bound = [100, 100, 100, 100] # junk values for testing
upper_bound = [100, 100, 100] # junk values for testing

#lower_bound = [ 0.0001, 0.0001, 0.0001, 0.0000001]
lower_bound = [ 0.0001, 0.0001, 0.0000001]

# Make the extrapolating version of the demographic model function
func_ex = dadi.Numerics.make_extrap_log_func(func)
# Calculate the model AFS
model = func_ex(params, ns, pts_l)

# Likelihood of the data given the model
ll_model = dadi.Inference.ll_multinom(model, fs2dIntergenic)

# The optimal value of theta given the model
theta = dadi.Inference.optimal_sfs_scaling(model,fs2dIntergenic)

# Perturb parameter array before optimization. In other words, take each parameter up or
# down by a factor of two or less
p0 = dadi.Misc.perturb_params(params,fold=1,upper_bound=upper_bound)

# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters.
# The maxiter argument restricts how long the optimizer will run. For production
# runs, you may want to set this value higher, to encourage better convergence.
popt = dadi.Inference.optimize_log(p0, fs2dIntergenic, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=500)
print 'Optimized parameters', repr(popt)
# update model and ll
model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, fs2dIntergenic)
print 'Optimized log-likelihood:', ll_opt
# The optimal value of theta given the model
theta = dadi.Inference.optimal_sfs_scaling(model,fs2dIntergenic)


# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure()
dadi.Plotting.plot_single_2d_sfs(theta*model,vmin=1,pop_ids=("TIL","BKN"))
dadi.Plotting.plot_2d_comp_multinom(model, fs2dIntergenic, pop_ids=("TIL","BKN"),vmin=1,residual='linear')
pylab.savefig('Bottleneck_teoFixed.png',dpi=300)







"""
Optimized parameters array([  5.54853085e-02,   1.79185626e+01,   4.38606255e-03,
         3.48659372e-02]) # decent set...

Optimized parameters array([  5.91906178e-02,   1.81637622e+01,   5.99466187e-03,
         2.44746920e-02]) #another decent set
"""

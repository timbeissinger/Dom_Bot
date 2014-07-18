##############################################################################
### This script runs a basic dadi bottleneck model using n ouput generated ###
### from angsd.                                                            ###
##############################################################################

### Timothy M. Beissinger
### 7-14-2014

import numpy
from numpy import array

import dadi

### Load the single-pop spectra
fsBKN = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/sfsBKN.dadi" )
fsTIL = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/sfsTIL.dadi" )

### Load a 2d fs (BKN and TIL)
fs2d = dadi.Spectrum.from_file("/Users/beissinger/Documents/DomesticationBottleneck/Dom_Bot/SFS/2d_TIL_BKN.dadi" )

### Time to plot the fs
import pylab
dadi.Plotting.plot_single_2d_sfs(fs2d,vmin=1,pop_ids=("TIL","BKN"))
pylab.show()

### Compute some summary statistics
S = fs2d.S() # Number of segregating sites
Fst = fs2d.Fst() # Fst

######################################
######### Better Bottleneck  #########
######################################

import realisticBottleneck_teoFixed

# determine sample sizes
ns = fs2d.sample_sizes
#ns= fsBKN.sample_sizes


# set grid sizes
pts_l = [80,90,100]

# Let's work with the "custom" model in the bottleneck.py file
#func = realisticBottleneck.maizeBottleneck
func = realisticBottleneck_teoFixed.maizeBottleneck
# ll for this model:
#params = array([1,1,1,1,1 ]) # junk values for testing
params = array([1,1,1,1 ]) # junk values for testing
#upper_bound = [10000, 100000, 3000, 2000, 2000] # junk values for testing
upper_bound = [100, 10000, 5, 5] # junk values for testing
#lower_bound = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001]
lower_bound = [ 0.00001, 0.00001, 0.00001, 0.00001]
# Make the extrapolating version of the demographic model function
func_ex = dadi.Numerics.make_extrap_log_func(func)
# Calculate the model AFS
model = func_ex(params, ns, pts_l)

# Likelihood of the data given the model
ll_model = dadi.Inference.ll_multinom(model, fs2d)

# The optimal value of theta given the model
theta = dadi.Inference.optimal_sfs_scaling(model,fs2d)

# Perturb parameter array before optimization. In other words, take each parameter up or
# down by a factor of two or less
p0 = dadi.Misc.perturb_params(params,fold=1,upper_bound=upper_bound)

# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters.
# The maxiter argument restricts how long the optimizer will run. For production
# runs, you may want to set this value higher, to encourage better convergence.
popt = dadi.Inference.optimize_log(p0, fs2d, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params),
                                   maxiter=5000)
print 'Optimized parameters', repr(popt)
# update model and ll
model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, fs2d)
print 'Optimized log-likelihood:', ll_opt
# The optimal value of theta given the model
theta = dadi.Inference.optimal_sfs_scaling(model,fs2d)


# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure()
dadi.Plotting.plot_single_2d_sfs(theta*model,vmin=1,pop_ids=("TIL","BKN"))
dadi.Plotting.plot_2d_comp_multinom(model, fs2d, pop_ids=("TIL","BKN"),vmin=1,residual='linear')
pylab.savefig('Bottleneck_teoFixed.png',dpi=300)







"""
Optimized parameters array([  5.54853085e-02,   1.79185626e+01,   4.38606255e-03,
         3.48659372e-02]) # decent set...

Optimized parameters array([  5.91906178e-02,   1.81637622e+01,   5.99466187e-03,
         2.44746920e-02]) #another decent set
"""

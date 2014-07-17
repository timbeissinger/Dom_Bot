"""
Realistic bottleneck model, no migration.
"""

"""
In words:
At time TF + TB in the past, an equilibrium teosinte population will split into
a maize population and a teosinte population. The maize population will go
through a bottleneck of size nuMB. Then, it will grow exponentially  at
time TF. The Teosinte population will not go through a bottleneck. It will end
at size nuTF.
"""
import numpy
import dadi

def maizeBottleneck ((nuMB , nuMF, nuTF , TB , TF ), ns, pts):
    # define the grid
    xx = yy = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # now do the population split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Send maize population through a bottleneck
    phi = dadi.Integration.two_pops(phi, xx, TB, nuMB, nuTF, m12=0, m21=0)

    # define the exponential growth function
    nu_func = lambda t: numpy.exp(numpy.log(nuMF) * t/TF)

    # Recover maize population exponentially
    phi = dadi.Integration.two_pops(phi, xx, TF, nu_func, nuTF, m12=0, m21=0)

    #phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    #phi = dadi.Integration.one_pop(phi,xx,TF,nuF)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,yy))
    return fs

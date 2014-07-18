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

def maizeBottleneck ((nuMB , nuMF, TF ), ns, pts):
    # define the grid
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # now do the population split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # define the exponential growth function
    nu_func = lambda t: nuMB*(nuMF/nuMB)**(t/TF)

    # Recover maize population exponentially, Teo constant
    phi = dadi.Integration.two_pops(phi, xx, TF, nu1=1, nu2=nu_func)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

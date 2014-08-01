"""
Realistic bottleneck model, no migration.
"""

"""
In words:
At time TF in the past, an equilibrium teosinte population will split into
a maize population and a teosinte population. The maize population will immediately go
through a bottleneck of size nuMB. Then, it will grow exponentially  at
time TF. The Teosinte population will not go through a bottleneck. It will remain at a
fixed size.
"""
import numpy
import dadi

def maizeBottleneck ((nuMB , nuMF, nuTeoF , TF ), ns, pts):
    # define the grid
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # now do the population split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Send maize population through a bottleneck
    phi = dadi.Integration.two_pops(phi, xx, nu1 = 1, nu2 = nuMB)

    # define the maize exponential growth function
    nu_funcM = lambda t: nuMB*(nuMF/nuMB) ** (t/TF)
    
    # define the teosinte linear growth function
    nu_funcTeo = lambda t: (nuTeoF-1) * (t/TF) + 1

    # Recover maize population exponentially
    phi = dadi.Integration.two_pops(phi, xx, T = TF, nu1 = nu_funcTeo, nu2 = nu_funcM)


    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

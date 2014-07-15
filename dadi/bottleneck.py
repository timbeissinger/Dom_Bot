"""
Simple bottleneck model
"""
import numpy
import dadi

def bottleneck ((nuB , nuF , TB , TF ), (ns), pts):
#    nuB , nuF , TB , TF = params
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi,xx,TF,nuF)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,)) 
    return fs

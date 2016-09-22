import numpy as np
import matplotlib.pyplot as plt

import SiteModel as sm

#--------------------------------------------------------------
# Manually creating a site model

site = sm.Site1D(Id=1,
                 Name='Test Site',
                 Latitude=10.,
                 Longitude=10.,
                 Elevation=[])

site.AddLayer({'Hl': 10.0,'Qp': 50.0,'Qs': 20.0,'Vp': 300.0,'Vs': 200.0,'Dn': 1900.0})
site.AddLayer({'Hl': 20.0,'Qp': 50.0,'Qs': 20.0,'Vp': 500.0,'Vs': 300.0,'Dn': 1900.0})
site.AddLayer({'Hl': 0.0,'Qp': 100.0,'Qs': 50.0,'Vp': 1000.0,'Vs': 800.0,'Dn': 2100.0})

hl = np.array(site.GetProfile('Hl'))
vs = np.array(site.GetProfile('Vs'))
dn = np.array(site.GetProfile('Dn'))

# Frequency axis
site.Freq = np.logspace(np.log10(0.1),
                        np.log10(10),
                        100)

print 'Vs30:', site.ComputeTTAV(), 'm/s'

import quartwlen as qwl

qwhl, qwvs, qwdn, ampf = qwl.QwlSolver(hl, vs, dn, site.Freq)
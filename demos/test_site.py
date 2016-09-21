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
site.AddLayer({'Hl': 10.0,'Qp': 50.0,'Qs': 20.0,'Vp': 500.0,'Vs': 300.0,'Dn': 1900.0})
site.AddLayer({'Hl': 0.0,'Qp': 100.0,'Qs': 50.0,'Vp': 1000.0,'Vs': 800.0,'Dn': 2100.0})

#--------------------------------------------------------------
# Alternative method using arbitrary parameter keys

site = sm.Site1D(Keys=['Depth','Vs'])

site.AddLayer([10.0, 200.0])
site.AddLayer([10.0, 300.0])
site.AddLayer([10.0, 800.0])

#--------------------------------------------------------------
# Parsing standar CSV data format

site = sm.Site1D()

site.ImportModel('data/site01.csv')

#--------------------------------------------------------------
# Parsing Geopsy data format

site = sm.Site1D(Keys=['Hl','Vp','Vs','Dn','Qp','Qs'])

site.ImportModel('data/site01.mod',
                  read_header='no',
                  delimiter=' ',
                  skipline=1)

#--------------------------------------------------------------
# Extracting profile information

print 'Vs Profile:', site.GetProfile('Vs')

#--------------------------------------------------------------
# Add sites to a regional container

zone = sm.SiteBox(Id=1,
                  Name='Site Collection')

zone.AddSite(site)
zone.AddSite(site)
zone.AddSite(site)

#--------------------------------------------------------------
# Computing Travel-Time average velocity

# Direct call (using default values)
print 'Vs30:', site.ComputeTTAV(), 'm/s'

# Using stored values for arbitrary depth
site.ComputeTTAV('Vs',50.)
print 'Vs50:', site.EngPar['Vz']['50.0'], 'm/s'

#--------------------------------------------------------------
# Compute Geotechnical soil class

site.ComputeGTClass('EC8')

print 'EC8 Site Class:', site.EngPar['EC8']

#--------------------------------------------------------------
# Computing SH transfer-function

# Frequency axis
site.Freq = np.logspace(np.log10(0.1),
                        np.log10(100),
                        1000)

site.ComputeSHTF()

plt.plot(site.Freq,np.abs(site.AmpFun['ShTF']))
plt.show(block=False)

#--------------------------------------------------------------
# Identify resonance frequencies and their amplitude

site.ComputeFnRes()

print 'F0 is:', site.AmpFun['Fn'][0][0], 'Hz'
print 'A0 is:', site.AmpFun['Fn'][1][0]
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

print 'Hl array:', site.GetProfile('Hl')
print 'Vs array:', site.GetProfile('Vs')

#--------------------------------------------------------------
# Add sites to a regional container

zone = sm.SiteBox(Id=1,
                  Name='Site Collection')

zone.AddSite(site)
zone.AddSite(site)
zone.AddSite(site)

#--------------------------------------------------------------
# Computing Travel-Time average velocity

# Direct call
print 'Vs30:', site.ComputeTTAV('Vs',30.)

# Using database
site.ComputeTTAV('Vs',50.)
print 'Vs50:', site.EngPar['Vz']['50.0']

#--------------------------------------------------------------
# Computing SH transfer-function

# Frequency axis
freq = np.logspace(np.log10(0.1),np.log10(100),1000)

shtf = site.ComputeSHTF(freq,0.)

plt.plot(freq,np.abs(shtf))
plt.show(block=False)
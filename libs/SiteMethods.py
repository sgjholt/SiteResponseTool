# -*- coding: utf-8 -*-

# Copyright (C) 2016 GEM Foundation
#
# SeismicSiteTool is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SeismicSiteTool is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.
#
# Author: Poggi Valerio

'''
Collection of standard functions for site response analysis
'''

import numpy as np
import scipy.optimize as spo


def TTAverageVelocity(hl, vs, z):
  '''
  The function calculates the travel-time
  average seismic velocity down to specified
  depth (e.g. the Vs30).

  Input parameters:
    hl = array of n layer thickness (m)
    vs = attay of n seismic velocities (m/s)
    z = averaging depth (m)

  Output:
    vsz = average velocity over depth 'z'
  '''

  # Initialisation
  lnum = len(hl)

  hl = np.array(hl)
  vs = np.array(vs)

  # Depth averaging is done on slowness
  vsz = 1./DepthAverage(lnum, hl, 1./vs, z)

  return vsz


def QwlApproxSolver(hl, vs, dn, fr):
  '''
  This function solves the quarter-wavelength problem
  (Boore 2003) and return the frequency-dependent
  depth, velocity, density and amplification factor.

  Input parameters:

    hl = vector of n thickness (m)
    vs = vector of n S-wave velocites (m/s)
    dn = vector of n densities (gr/m3)
    fr = vector of discrete frequencies (Hz)

  Output:

    qwhl = vector of quarter-wavelength depths
    qwvs = vector of quarter-wavelength velocities
    qwdn = vector of quarter-wavelength densities
    qwaf = vector of quarter-wavelength amp. factors
  '''

  # Initialisation
  fnum = len(fr)
  lnum = len(hl)

  hl = np.array(hl)
  vs = np.array(vs)
  dn = np.array(dn)

  qwhl = np.zeros(fnum)
  qwvs = np.zeros(fnum)
  qwdn = np.zeros(fnum)
  qwaf = np.zeros(fnum)

  # Rock reference (last layer)
  refv = vs[-1]
  refd = dn[-1]

  for nf in range(fnum):

    # Upper depth bound for the search
    ubnd = np.max(vs)/(4.*fr[nf])

    # Search for quarter-wavelength depth
    qwhl[nf] = spo.fminbound(QwlFitFunc, 0., ubnd,
                             args=(lnum,hl,1./vs,fr[nf]))

    # Computing average velocity (note: slowness is used)
    qwvs[nf] = 1./DepthAverage(lnum, hl, 1./vs, qwhl[nf])

    # Computing average density (for amplification function)
    qwdn[nf] = DepthAverage(lnum, hl, dn, qwhl[nf])

    # Computing amplification function
    qwaf[nf] = np.sqrt((refd*refv)/(qwdn[nf]*qwvs[nf]))

  return qwhl, qwvs, qwdn, qwaf


def QwlFitFunc(z, lnum, hl, sl, fr):
  '''
  Misfit function (simple L1 norm)
  '''

  qwsl = DepthAverage(lnum, hl, sl, z)
  obj = np.abs(z-(1/(4.*fr*qwsl)))
  
  return obj


def DepthAverage(lnum, hl, par, z):
  '''
  Search function to compute depth-weighted
  average of a generic parameter (slowness, density...)
  '''

  ztot = 0.
  sum = 0.
  islast = False

  for nl in range(lnum):

    if ((ztot + hl[nl]) < z) and (nl != lnum-1):
      sum += (hl[nl]*par[nl])

    else:
      if islast == False:
        sum += (z - ztot)*par[nl]
        islast = True

    ztot += hl[nl]

  return sum/z


def ShTransferFunction(Hl, Vs, Dn, Qs, Freq, Iang=0., Elastic=False):

  """
  SH wave transfer function using Knopoff formalism.
  Authors: Poggi Valerio, Marwan Irnaka
  """

  # Variable recasting
  nlayer = len(Hl)
  hl = np.array(Hl,dtype='complex128')
  vs = np.array(Vs,dtype='complex128')
  dn = np.array(Dn,dtype='complex128')
  qs = np.array(Qs,dtype='complex128')
  freq = np.array(Freq)
  iang = np.array(Iang)

  # Angular frequency conversion
  angf = 2.*np.pi*freq

  # Attenuation using complex velocities
  if not Elastic:
    vs = vs*((2.*qs*1j)/(2.*qs*1j-1.))

  # Angle of propagation within layers
  iD = np.zeros((nlayer,1))
  iCORE = np.zeros((nlayer,nlayer),dtype='complex128')

  iD[0] = np.sin(iang)
  iCORE[0,-1] = 1.

  for nl in range(nlayer-1):
      iCORE[nl+1,nl] = 1./vs[nl]
      iCORE[nl+1,nl+1] = -1./vs[nl+1]

  iA = np.linalg.solve(iCORE,iD)

  iS = np.arcsin(iA)

  # Lame Parameter(s)
  mu = np.zeros((nlayer,1),dtype='complex128')

  for nl in range(nlayer):
      mu[nl]=dn[nl]*(vs[nl]**2)

  # Horizontal and vertical slowness
  ns = np.zeros((nlayer,1),dtype='complex128')

  for nl in range(nlayer):
      ns[nl]=np.cos(iS[nl])/vs[nl]

  # Building data vector
  A = np.zeros((nlayer*2,1))
  D = np.zeros((nlayer*2,1))
  D[-1] = 1.

  # Dispacement and transfer function initialisation
  fnum = len(Freq)
  htu = np.zeros((fnum,1),dtype='complex128')
  hbu = np.zeros((fnum,1),dtype='complex128')
  htf = np.zeros((fnum,1),dtype='complex128')

  for nf, af in enumerate(angf):

    # Building core matrix
    CORE = np.zeros((nlayer*2,nlayer*2),dtype='complex128')

    # Free surface constraints
    CORE[0,0] = 1.
    CORE[0,1] = -1.

    # Interfaces constraints
    for nl in range(nlayer-1):

      row = (nl*2)+1
      col = nl*2

      expDSA = np.exp(1j*af*ns[nl]*hl[nl])
      expUSA = np.exp(-1j*af*ns[nl]*hl[nl])

      CORE[row,col+0] = expDSA[0]
      CORE[row,col+1] = expUSA[0]
      CORE[row,col+2] = -1.
      CORE[row,col+3] = -1.

      CORE[row+1,col+0] =  mu[nl][0]*ns[nl][0]*expDSA[0]
      CORE[row+1,col+1] = -mu[nl][0]*ns[nl][0]*expUSA[0]
      CORE[row+1,col+2] = -mu[nl+1][0]*ns[nl+1][0]
      CORE[row+1,col+3] =  mu[nl+1][0]*ns[nl+1][0]

    # input constraints
    CORE[-1,-1] = 1.

    # solving linear system
    try:
        A = np.linalg.solve(CORE,D)
    except:
        A[:] = np.nan

    # Computing displacements
    htu[nf] = A[0]+A[1]
    hbu[nf] = 2*A[-1]
    htf[nf] = htu[nf]/hbu[nf]

  return htf


def GetResFreq(Freq, AmpF):
  '''
  Identify resonance frequencies of an amplification function.
  Output are two arrays: resonce frequencies and amplitude maxima.
  '''

  Fn = np.array([])
  An = np.array([])

  for nf, fr in enumerate(Freq[:-2]):

    # Three-points search for local maxima
    a0 = AmpF[nf]
    a1 = AmpF[nf+1]
    a2 = AmpF[nf+2]

    if (a1-a0) > 0 and (a2-a1) < 0:

      # Storing value
      Fn = np.append(Fn,fr)
      An = np.append(An,a1)

  return Fn, An


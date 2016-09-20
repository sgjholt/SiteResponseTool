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

def TTAverageVelocity(Hl, Vl, Z):
  '''
  The function calculates the travel-time
  average seismic velocity down to specified
  depth (e.g. the Vl30).

  Input parameters:
    Hl = array of n layer thickness (m)
    Vl = attay of n seismic velocities (m/s)
    Z = averaging depth (m)

  Output:
    VlZ = average over depth 'Z'
  '''

  # Initialisation
  tt = 0.
  Htot = 0.
  Hl = np.array(Hl)
  Vl = np.array(Vl)
  lnum = Hl.size;

  # Check if the model is an homogenous half-space
  if lnum == 1:
    VlZ = Vl
  else:
    # Loop over layers to compute the travel-times(tt)
    for nl in range(0,lnum):
      # Check for special cases (last and mid layers)
      if nl != (lnum-1) and (Htot+Hl[nl]) <= Z:
        tt += Hl[nl]/Vl[nl]
      else:
        tt += (Z-Htot)/Vl[nl]
        break
      Htot += Hl[nl]
    VlZ = Z/tt

  return VlZ


def ShTransferFunction(Hl, Vs, Dn, Qs, Freq, Iang):

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
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
Collection of functions for site analysis
'''

#--------------------------------------------------

import numpy as np

#--------------------------------------------------

def TTAverageVelocity(Hl, Vs, Z):

  '''
  The function calculates the travel-time
  average seismic velocity down to specified
  depth (e.g. the Vs30).

  Input parameters:
    Hl = array of n layer thickness (m)
    Vs = attay of n seismic velocities (m/s)
    Z = averaging depth (m)

  Output:
    VsZ = average over depth 'Z'
  '''

  # Initialisation
  tt = 0
  Htot = 0
  Hl = np.array(Hl)
  Vs = np.array(Vs)
  lnum = Hl.size;

  # Check if the model is an homogenous half-space
  if lnum == 1:
    VsZ = Vs
  else:
    # Loop over layers to compute the travel-times(tt)
    for nl in range(0,lnum):
      # Check for special cases (last and mid layers)
      if nl != (lnum-1) and (Htot+Hl[nl]) <= Z:
        tt += Hl[nl]/Vs[nl]
      else:
        tt += (Z-Htot)/Vs[nl]
        break
      Htot += Hl[nl]
    VsZ = Z/tt

  return VsZ
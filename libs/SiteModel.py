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
Collection of base classes to create site models.
'''

#--------------------------------------------------

import numpy as np
import SiteMethods as SM

#--------------------------------------------------

class Site1D(object):

  def __init__(self, Id=[],
                     Name=[],
                     Longitude=[],
                     Latitude=[]):

    self.meta = {'Id': Id,
                 'Name': Name,
                 'Longitude': Longitude,
                 'Latitude': Latitude}

    self.layer = []
    self.layer_number = 0

  #--------------------------------------------------

  def AddLayer(self, data=[]):

    # Create a basic empty layer structure
    L = {'Hl': [],
         'Vp': [],
         'Vs': [],
         'Qp': [],
         'Qs': []}

    # Inflate structure with data (if any)
    if data:
      for k in data.keys():
        L[k] = float(data[k])

    # Add the layer to the model
    self.layer.append(L)
    self.layer_number += 1

  #--------------------------------------------------

  def GetProfile(self, key):

    return np.array([i[key] for i in self.layer])

  #--------------------------------------------------

  def TTAverageVelocity(self, key, Z):

    Hl = self.GetProfile('Hl')
    Vs = self.GetProfile(key)
    VsZ = SM.TTAverageVelocity(Hl, Vs, Z)

    return VsZ

  #--------------------------------------------------

  def ImportAsciiModel(self, model_file,
                             header=[],
                             delimeter=',',
                             skipline=0,
                             comment='#'):

    # Opening input model
    with open(model_file, 'r') as f:

      # Reinitialise layer parameters
      self.layer = []
      self.layer_number = 0

      # Read and ignore initial lines
      for i in range(0, skipline):
        f.readline()

      # Import header line if not specified
      if not header:
        line = f.readline()
        header = line.strip().split(delimeter)

      # Loop over data
      for line in f:

         # Skip comment lines
        if line[0] != comment:

          line = line.strip().split(delimeter)
          value = np.array(line, dtype=float)

          # Loop over header keys
          data = {}
          for i, k in enumerate(header):
            data[k] = value[i]

          self.AddLayer(data)

      f.close()
      return

    # If file does not exist
    print 'File not found.'

#--------------------------------------------------

class SiteBox(object):

  def __init__(self, Id=[],
                     Name=[]):

    self.meta = {'Id': Id,
                 'Name': Name}

    self.site = []
    self.site_number = 0

  #--------------------------------------------------

  def AddSite(self, Site1D):

    self.site.append(Site1D)
    self.site_number += 1
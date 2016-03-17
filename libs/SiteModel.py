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
# along with OpenQuake.  If not, see <http://www.gnu.org/licenses/>.

'''
Collection of base classes to create site models.
'''

#--------------------------------------------------

import numpy as np

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

  def AddLayer(self, Hl=[],
                     Vp=[],
                     Vs=[],
                     Qp=[],
                     Qs=[]):

    L = {'Hl': Hl,
         'Vp': Vp,
         'Vs': Vs,
         'Qp': Qp,
         'Qs': Qs}

    self.layer.append(L)
    self.layer_number += 1

  #--------------------------------------------------

  def GetProfile(self, key):

    return np.array([i[key] for i in self.layer])

  #--------------------------------------------------

  def ImportAsciiModel(self, model_file, delimeter=' ', skipline=0):

    with open(model_file, 'r') as f:

      self.layer = []
      self.layer_number = 0

      # Read and ignore header lines
      for i in range(0, skipline):
        f.readline()

      for i, line in enumerate(f):

        line = line.strip().split(delimeter)
        value = np.array(line, dtype=float)

        self.AddLayer(value[0],
                      value[1],
                      value[2],
                      value[3],
                      value[4])

      f.close()

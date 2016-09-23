# -*- coding: utf-8 -*-
#
# Copyright (C) 2010-2016 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.
#
# Author: Poggi Valerio

'''
Collection of base classes to create and manipulate site models.
'''

import numpy as np
import SiteMethods as SM
import AsciiTools as AT

class Site1D(object):
  '''
  SINGLE SITE ELEMENT (ONE-DIMENSIONAL)
  The object is structured to act as minimal database for
  soil properties, site metadata and methods to compute
  site parameters (Vs30, amplification...)
  '''


  def __init__(self, Id=[],
                     Name=[],
                     Longitude=[],
                     Latitude=[],
                     Elevation=[],
                     Keys=['Hl','Vp','Vs','Dn','Qp','Qs']):

    self.Meta = {'Id': Id,
                 'Name': Name,
                 'Longitude': Longitude,
                 'Latitude': Latitude,
                 'Elevation': Elevation}

    self.Keys = Keys
    self.Layer = []
    self.Freq = []
    self.EngPar = {}
    self.AmpFun = {}


  def AddLayer(self, Data=[]):
    '''
    Method to add a single layer (and its properties)
    to the site structure, at the bottom of an existing stack.
    Data can be a list of values, sorted according to Site1D.Keys
    or a dictionary element with corresponding format.
    '''

    LS = {}

    # Inflate the layer structure with data (if any)

    if type(Data) is list:
      for i, k in enumerate(self.Keys):
        if Data:
          LS[k] = Data[i]
        else:
          LS[k] = []

    if type(Data) is dict:
      for k in Data.keys():
        LS[k] = Data[k]

    # Add the new layer structure to the layer stack
    self.Layer.append(LS)


  def Size(self):
    '''
    Method to return size of the data matrix.
    '''

    lnum = len(self.Layer)
    knum = len(self.Keys)

    return [lnum, knum]


  def ImportModel(self, ascii_file,
                        read_header='yes',
                        dtype='float',
                        delimiter=',',
                        skipline=0,
                        comment='#'):
    '''
    Method to parse soil properties from an ascii file.
    '''

    if read_header == 'yes':
      header = []

    if read_header == 'no':
      header = self.Keys

    at = AT.AsciiTable()
    at.Import(ascii_file, header=header,
                          dtype=dtype,
                          delimiter=delimiter,
                          skipline=skipline,
                          comment=comment)

    self.Keys = at.header
    self.Layer = at.data


  def GetProfile(self, key):
    '''
    Utility method to extract a column of soil properties
    from the layer stack. It returns a numpy array.
    '''

    return np.array([i[key] for i in self.Layer])


  def ComputeTTAV(self, key='Vs', Z=30.):
    '''
    Compute and store travel-time average velocity at
    a given depth (Z) and for a specific key.
    Default is Vs30
    '''

    # Formatting model parameters
    Hl = self.GetProfile('Hl')
    Vl = self.GetProfile(key)

    # Compute average velocity
    Vz = SM.TTAverageVelocity(Hl, Vl, Z)

    # Check data structure
    if 'Vz' not in self.EngPar:
      self.EngPar['Vz'] = {}

    # Store results into database
    self.EngPar['Vz'][str(Z)] = Vz

    return Vz


  def ComputeQWL(self, key='Vs'):
    '''
    Compute and store the quarter-wavelength representation
    of the soil profile. Qwl-amplification is also stored.
    '''

    # Formatting model parameters
    Hl = self.GetProfile('Hl')
    Vl = self.GetProfile(key)
    Dn = self.GetProfile('Dn')

    # Compute average velocity
    QwHl, QwVs, QwDn, QwAf = SM.QwlApproxSolver(Hl, Vl, Dn, self.Freq)

    # Check data structure
    if 'Qwl' not in self.EngPar:
      self.EngPar['Qwl'] = {}

    # Store results into database
    self.EngPar['Qwl']['Hl'] = QwHl
    self.EngPar['Qwl']['Vs'] = QwVs
    self.EngPar['Qwl']['Dn'] = QwDn

    # Check data structure
    if 'ShTF' not in self.AmpFun:
      self.AmpFun['Qwl'] = {}

    # Store results into database
    self.AmpFun['Qwl'] = QwAf

    return QwHl, QwVs, QwDn, QwAf


  def ComputeKappa0(self, key='Vs', Z=[]):
    '''
    Compute the Kappa parameter from the Qs profile
    of the site, down to a given depth (default is
    the whole profile)
    '''

    # Formatting model parameters
    Hl = self.GetProfile('Hl')
    Vl = self.GetProfile(key)
    Qs = self.GetProfile('Qs')

    lnum = len(Hl)

    # If Z not given, using the whole profile
    if not Z:
      Z = np.sum(Hl)

    Par = Z/(Vl*Qs)
    Kappa0 = SM.DepthAverage(lnum, Hl, Par, Z)

    # Check data structure
    if 'Kappa0' not in self.EngPar:
      self.EngPar['Kappa0'] = {}

    # Store results into database
    self.EngPar['Kappa0'] = Kappa0

    return Kappa0

  def ComputeGTClass(self,BCode='EC8'):
    '''
    Compute geotechnical classification according
    to specified building code. Default is EC8.
    '''

    Vs30 = self.EngPar['Vz']['30.0']

    # Check data structure
    if BCode not in self.EngPar:
      self.EngPar[BCode] = {}

    # Check the values
    if Vs30 >= 800.:
      GClass = 'A'
    if Vs30 >= 360. and Vs30 < 800.:
      GClass = 'B'
    if Vs30 >= 180. and Vs30 < 360.:
      GClass = 'C'
    if Vs30 < 180.:
      GClass = 'D'

    self.EngPar['EC8'] = GClass

    return GClass


  def ComputeSHTF(self, Iang=0., Elastic=False):
    '''
    Compute the SH transfer function for an arbitrary
    incidence angle. Default angle is 0.
    '''

    # Formatting model parameters
    Hl = self.GetProfile('Hl')
    Vs = self.GetProfile('Vs')
    Dn = self.GetProfile('Dn')
    Qs = self.GetProfile('Qs')

    # TF calculation
    ShTF = SM.ShTransferFunction(Hl, Vs, Dn, Qs, self.Freq, Iang, Elastic)

    # Check data structure
    if 'ShTF' not in self.AmpFun:
      self.AmpFun['ShTF'] = {}

    # Store results into database
    self.AmpFun['ShTF'] = ShTF

    return ShTF


  def ComputeFnRes(self):
    '''
    Identify resonance frequencies of the SH-wave
    transfer function.
    '''

    Fn, An = SM.GetResFreq(self.Freq, np.abs(self.AmpFun['ShTF']))

    # Check data structure
    if 'Fn' not in self.AmpFun:
      self.AmpFun['Fn'] = {}

    # Store results into database
    self.AmpFun['Fn'] = Fn

    # Check data structure
    if 'An' not in self.AmpFun:
      self.AmpFun['An'] = {}

    # Store results into database
    self.AmpFun['An'] = An

    return Fn, An


  def ComputeAttFun(self):
    '''
    Compute the frequency-dependent attenuation
    functio from a give Kappa0.
    '''

    # Compute exponential decay function
    AttF = np.exp(-np.pi*self.EngPar['Kappa0']*self.Freq)

    # Check data structure
    if 'AttF' not in self.AmpFun:
      self.AmpFun['AttF'] = {}

    # Store results into database
    self.AmpFun['AttF'] = AttF

    return AttF


  def GenFreqAx(self, Fmin=0.1, Fmax=10., Fnum=100, Log=True):
    '''
    Method to generate a lin/log spaced frequency axis.
    '''

    if Log:
      self.Freq = np.logspace(np.log10(Fmin),
                              np.log10(Fmax),
                              Fnum)
    else:
      self.Freq = np.linspace(Fmin, Fmax, Fnum)


class SiteBox(object):
  '''
  A simple container class to group sites of a region
  '''

  def __init__(self, Id=[],
                     Name=[]):

    self.Meta = {'Id': Id,
                 'Name': Name}

    self.Site = []
    self.Size = 0


  def AddSite(self, Site1D):
    '''
    Method to just add a single site object to the container.
    '''

    self.Site.append(Site1D)
    self.Size += 1


def AddDict(target, key=[], value=[]):
  '''
  Test function.
  '''

  if not target:
    target = {}

  if key:
    if key not in target:
      target[key] = {}
      if value:
        target[key] = value

  return target

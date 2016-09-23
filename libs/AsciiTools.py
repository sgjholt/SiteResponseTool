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

"""
A simple tool to manipulate data from/to ascii files.
"""

import numpy as np
import fnmatch as fnm

class AsciiTable():

  def __init__(self, header=[]):

    if header:
      self.header = header
    else:
      self.header = []

    self.data = []


  def AddElement(self, data=[]):
    """
    Add an element (with header's format) to the data structure.
    Element can be empty or filled with data.
    """

    newitem = {}

    for i, key in enumerate(self.header):

      if not data:
        newitem[key] = []
      else:
        newitem[key] = data[i]

    self.data.append(newitem)


  def Size(self):
    """
    Method to return size of the data matrix.
    """

    enum = len(self.data)
    hnum = len(self.header)

    return [enum, hnum]


  def Import(self, ascii_file,
                   header=[],
                   dtype='float',
                   delimiter=',',
                   skipline=0,
                   comment='#'):
    """
    Method to import data from ascii file (tabular)
    """

    # Open input ascii file
    with open(ascii_file, 'r') as f:

      # Ignore initial line(s) if necessary
      for i in range(0, skipline):
        f.readline()

      # Import header if not passed
      if not header:
        line = f.readline()
        self.header = line.strip().split(delimiter)
      else:
        self.header= header

      # Loop over lines
      for line in f:

         # Skip comments, if any
        if line[0] != comment:
          value = line.strip().split(delimiter)

          # Loop over data values
          data = []
          for i, k in enumerate(self.header):
            data.append(CastValue(value[i],dtype))

          self.AddElement(data)

      f.close()
      return

    # Warn user if model file does not exist
    print 'File not found.'


  def Export(self, ascii_file,
                   write_header='yes',
                   delimiter=','):
    """
    Method to export data object into an ascii file.
    """

    with open(ascii_file, 'w') as f:

      # Write header
      if write_header == 'yes':
        header = delimiter.join(self.header)
        f.write(header + '\n')

      # Write data (loop over rows)
      for i, item in enumerate(self.data):
        data = delimiter.join([str(item[j]) for j in self.header])

        if i < (self.Size()[0]-1):
          f.write(data + '\n')
        else:
          f.write(data)

      f.close()
      return

    # Warn user if model file does not exist
    print 'File not found.'


  def Append(self, new_table):
    """
    Method to merge two data structures consecutively.
    (Header structure must be identical).
    """

    if self.header == new_table.header:
      for i in range(0,new_table.Size()[0]):
        self.data.append(new_table.data[i])

    else:
      print 'Error: headers do not match...'


  def Extract(self, key, dtype='float'):
    """
    Method to extract data values by key.
    Data type can be specified.
    """

    values = []

    for item in self.data:
      value = CastValue(item[key], dtype)
      values.append(value)

    return values


  def Filter(self, key, filter_key):
    """
    Method to filter the data table by key value.
    Value can be a string to amtch (* and ? allowed)
    or a numerical range (as a list of floats).
    In output it is returned a new table.
    """

    NewTab = AsciiTable(self.header)

    # String matching
    if type(filter_key) is str:

      for item in self.data:
        if fnm.fnmatch(item[key],filter_key):
           NewTab.data.append(item)

    # Filter by value
    if type(filter_key) is list:

      for item in self.data:
        if float(item[key]) >= filter_key[0] and \
           float(item[key]) < filter_key[1]:
          NewTab.data.append(item)

    return NewTab


def CastValue(value, dtype='float'):

  if dtype == 'string':
    value = str(value)
  if dtype == 'int':
    value = int(value)
  if dtype == 'float':
    value = float(value)

  return value
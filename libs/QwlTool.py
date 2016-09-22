import numpy as np
import scipy.optimize as spo

def QwlSolver(th, vs, dn, fr):
  '''
  This function solves the quarter-wavelength problem
  (Boore 2003) and return the frequency-dependent
  depth, velocity, density and amplification factor.

  Parameters:

    th = vector of n thickness (m)
    vs = vector of n S-wave velocites (m/s)
    dn = vector of n densities (gr/m3)
    fr = vector of discrete frequencies (Hz)

  Output:

    qwth = vector of quarter-wavelength depths
    qwvs = vector of quarter-wavelength velocities
    qwdn = vector of quarter-wavelength densities
    ampf = vector of quarter-wavelength amp. factors
  '''

  # Initialisation
  fnum = len(fr)
  lnum = len(th)

  th = np.array(th)
  vs = np.array(vs)
  dn = np.array(dn)

  qwth = np.zeros(fnum)
  qwvs = np.zeros(fnum)
  qwdn = np.zeros(fnum)
  ampf = np.zeros(fnum)

  # Rock reference (last layer)
  refv = vs[-1]
  refd = dn[-1]

  for nf in range(fnum):

    # Upper depth bound for the search
    ubnd = np.max(vs)/(4.*fr[nf])

    # Search for quarter-wavelength depth
    qwth[nf] = spo.fminbound(QwlFitFunc,0.,ubnd,args=(lnum,th,1./vs,fr[nf]))

    # Computing average velocity (note: slowness is used)
    qwvs[nf] = 1./QwlForward(qwth[nf],lnum,th,1./vs)

    # Computing average density (for amplification function)
    qwdn[nf] = QwlForward(qwth[nf],lnum,th,dn)

    # Computing amplification function
    ampf[nf] = np.sqrt((refd*refv)/(qwdn[nf]*qwvs[nf]))

  return qwth, qwvs, qwdn, ampf


def QwlFitFunc(z, lnum, th, sl, fr):
  '''
  Misfit function (simple L1 norm)
  '''

  qwsl = QwlForward(z, lnum, th, sl)
  obj = np.abs(z-(1/(4.*fr*qwsl)))
  
  return obj


def QwlForward(z, lnum, th, par):
  '''
  Search function to compute depth-weighted
  average of a generic parameter.
  '''

  ztot = 0.
  sump = 0.
  islast = False

  for nl in range(lnum):

    if ((ztot + th[nl]) < z) and (nl != lnum-1):
      sump += (th[nl]*par[nl])

    else:
      if islast == False:
        sump += (z - ztot)*par[nl]
        islast = True

    ztot += th[nl]

  return sump/z

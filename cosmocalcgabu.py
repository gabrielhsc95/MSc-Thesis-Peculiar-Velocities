#!/usr/bin/env python
"""
Calculate useful values for a given cosmology.  This module uses code adapted
from `CC.py`_ (`James Schombert`_) which is a Python version of the
`Cosmology Calculator`_ (`Ned Wright`_).

The following values are calculated:

    ====  ===================================  ===========
    Name  Value                                Units
    ====  ===================================  ===========
    z     Input redshift
    H0    Hubble constant
    WR    Omega(radiation)
    WK    Omega curvaturve = 1-Omega(total)
    WM    Omega matter
    WV    Omega vacuum
    DTT   Time from z to now                   Gyr
    age   Age of Universe                      Gyr
    zage  Age of Universe at redshift z        Gyr
    DCMR  Comoving radial distance             Gyr Mpc cm
    VCM   Comoving volume within redshift      Gpc3
    DA    Angular size distance                Gyr Mpc cm
    DL    Luminosity distance                  Gyr Mpc cm
    PS    Plate scale - distance per arcsec    kpc cm
    ====  ===================================  ===========

.. _`James Schombert`: http://abyss.uoregon.edu/~js/
.. _`CC.py`: http://www.astro.ucla.edu/~wright/CC.python
.. _`Ned Wright`: http://www.astro.ucla.edu/~wright/intro.html
.. _`Cosmology Calculator`: http://www.astro.ucla.edu/~wright/CosmoCalc.html

:Copyright: Smithsonian Astrophysical Observatory (2009)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

import numpy as np

# Define a few constants
cm_per_pc = 3.0856775813057289536e+18
c = 299792.458                         # velocity of light in km/sec
km_per_ly = 3600*24*365.25*c           # km per light-year
Tyr = 977.8                            # coefficent for converting 1/H into Gyr
arcsec_per_rad = 206264.806
_outvals_str = ('z H0 WM WV WK WR',
                'DA DA_Gyr DA_Mpc DA_cm',
                'DL DL_Gyr DL_Mpc DL_cm',
                'DCMR DCMR_Gyr DCMR_Mpc DCMR_cm',
                'PS_kpc PS_cm',
                'DTT DTT_Gyr',
                'VCM VCM_Gpc3',
                'age age_Gyr',
                'zage zage_Gyr',)
_outvals = (' '.join(_outvals_str)).split()

def cosmocalcgabu(z, H0=71, WM=0.27, WV=None):
    """
    Calculate useful values for the supplied cosmology.

    This routine returns a dictionary of values in the form ``<name>: <value>``,
    where the values are supplied in "natural" units for cosmology, e.g. 1/H0.
    In addition various useful unit conversions are done and stored in the
    dictionary as ``<name>_<unit>: <value>``.  E.g. angular size distance::

      'DA': 0.38250549415474988,
      'DA_Gyr': 5.2678010166833023,
      'DA_Mpc': 1615.1022857909447,
      'DA_cm': 4.9836849147807571e+27

    Example::

     >>> from cosmocalc import cosmocalc
     >>> from pprint import pprint
     >>> pprint(cosmocalc(3, H0=75, WM=.25))
     {'DA': 0.39103776375786625,
      'DA_Gyr': 5.0980896720325548,
      'DA_Mpc': 1563.0689649039205,
      'DA_cm': 4.8231268630387788e+27,
      'DCMR': 1.564151055031465,
      'DCMR_Gyr': 20.392358688130219,
      'DCMR_Mpc': 6252.2758596156818,
      'DCMR_cm': 1.9292507452155115e+28,
      'DL': 6.25660422012586,
      'DL_Gyr': 81.569434752520877,
      'DL_Mpc': 25009.103438462727,
      'DL_cm': 7.717002980862046e+28,
      'DTT': 0.84826379084317027,
      'DTT_Gyr': 11.059097795819358,
      'H0': 75,
      'PS_cm': 2.3383178917293232e+22,
      'PS_kpc': 7.5779721961095019,
      'VCM': 1.2756009121294902,
      'VCM_Gpc3': 1023.7714254161302,
      'WK': 0.0,
      'WM': 0.25,
      'WR': 7.4044444444444448e-05,
      'WV': 0.74992595555555552,
      'age': 1.0133755371756261,
      'age_Gyr': 13.211714670004362,
      'z': 3,
      'zage': 0.16511174633245579,
      'zage_Gyr': 2.1526168741850036}

    :param z: redshift
    :param H0: Hubble constant (default = 71)
    :param WM: Omega matter (default = 0.27)
    :param WV: Omega vacuum (default = 1.0 - WM - 0.4165/(H0*H0))

    :rtype: dictionary of cosmology values (name_unit = value)
    """
    z = np.array(z)

    np.where(z > 100, z, z/299792.458) # Values over 100 are in km/s

    z = z.reshape((z.size, 1))

    if WV is None:
        WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda

    h = H0/100.0
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1.0 - WM - WR - WV
    az = 1.0/(1.0 + 1.0*z)

    n=1000         # number of points in integrals
    points = np.arange(n)
    a = az*(points + 0.5)/n
    adot = np.sqrt(WK + (WM/a) + (WR/np.power(a, 2)) + (WV*np.power(a, 2)))
    age = np.sum(1.0/adot, axis=1)
    zage = (az[:].T*age[:])/n

    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    a = az + (1 - az)*(points + 0.5)/n
    adot = np.sqrt(WK + (WM/a) + (WR/np.power(a, 2)) + (WV*np.power(a, 2)))
    DTT = np.sum(1.0/adot, axis=1)
    DCMR = np.sum(1.0/(a*adot), axis=1)

    DTT = ((1.0 - az)[:].T*DTT[:])/n
    DCMR = ((1.0 - az)[:].T*DCMR[:])/n
    age = DTT + zage

    # tangential comoving distance
    # comoving volume computation
    x = np.sqrt(np.absolute(WK))*DCMR
    x = x.reshape((x.size,))
    yes = np.where(x > 0.1)
    if WK > 0:
        tan_ratio_yes = np.sinh(x[yes])/x[yes]
        vol_ratio_yes = (np.sinh(2.0*x[yes])/4.0 - x[yes]/2.0)/(np.power(x[yes], 3)/3.0)
    else:
        tan_ratio_yes = np.sin(x[yes])/x[yes]
        vol_ratio_yes = (x[yes]/2.0 - np.sin(2.0*x[yes])/4.0)/(np.power(x[yes], 3)/3.0)

    no = np.where(x <= 0.1)
    y = np.power(x[no], 2)
    if WK < 0:
        y = -y
    tan_ratio_no = 1.0 + y/6.0 + np.power(y, 2)/120.0
    vol_ratio_no = 1.0 + y/5.0 + (2.0/105.0)*np.power(y, 2)

    tan_ratio = np.zeros(x.size)
    tan_ratio[yes] = tan_ratio_yes
    tan_ratio[no] = tan_ratio_no

    vol_ratio = np.zeros(x.size)
    vol_ratio[yes] = vol_ratio_yes
    vol_ratio[no] = vol_ratio_no

    DCMT = tan_ratio*DCMR

    VCM = vol_ratio*np.power(DCMR, 3)/3.0 # <<<
    VCM_Gpc3 = 4.0*np.pi*np.power((0.001*c/H0), 3)*VCM

    DA = az[:].T*DCMT[:]
    DL = (1/np.power(az, 2))[:].T*DA[:]

    # Now convert to some more useful units
    Gyr = lambda x: Tyr/H0*x
    Mpc = lambda x: c/H0*x
    cm = lambda x: Mpc(x)*1e6*cm_per_pc

    # Reshaping data for better prints
    z = z.reshape((z.size,))
    DA = DA.reshape((DA.size,))
    DL = DL.reshape((DL.size,))
    DCMR = DCMR.reshape((DCMR.size,))
    DTT = DTT.reshape((DTT.size,))
    age = age.reshape((age.size,))
    zage = zage.reshape((zage.size,))
    VCM = VCM.reshape((VCM.size,))
    VCM_Gpc3 = VCM_Gpc3.reshape((VCM_Gpc3.size,))


    DA_Gyr = Gyr(DA)
    DA_Mpc = Mpc(DA)
    DA_cm = cm(DA)

    DL_Gyr = Gyr(DL)
    DL_Mpc = Mpc(DL)
    DL_cm = cm(DL)

    DCMR_Gyr = Gyr(DCMR)
    DCMR_Mpc = Mpc(DCMR)
    DCMR_cm = cm(DCMR)

    DTT_Gyr = Gyr(DTT)
    age_Gyr = Gyr(age)
    zage_Gyr = Gyr(zage)

    PS_kpc = Mpc(DA)*1000/arcsec_per_rad
    PS_cm = PS_kpc*cm_per_pc*1000



    localvals = locals()
    return dict((x, localvals[x]) for x in _outvals)

def get_options():
    """
    cosmocalc.py [options] redshift [name_unit [name_unit2 ...]]

    Allowed ``name_unit`` values::

      DA DA_Gyr DA_Mpc DA_cm
      DL DL_Gyr DL_Mpc DL_cm
      DCMR DCMR_Gyr DCMR_Mpc DCMR_cm
      PS_kpc PS_cm
      DTT DTT_Gyr
      VCM VCM_Gpc3
      age age_Gyr
      zage zage_Gyr
      H0 WM WV WK WR z

    If no ``name_unit`` values are supplied then all the above will be printed."""
    from optparse import OptionParser

    parser = OptionParser(get_options.__doc__)
    parser.set_defaults()
    parser.add_option("--H0",
                      default=None,
                      type='float',
                      help="Hubble constant")
    parser.add_option("--WM",
                      default=None,
                      type='float',
                      help="")
    parser.add_option("--WV",
                      default=None,
                      type='float',
                      help="")
    opt, args = parser.parse_args()
    return opt, args, parser

def cosmocalc_DL_Mpc(z, H0=71, WM=0.27, WV=None, n=1000):
    z = np.array(z)

    z = z.reshape((z.size, 1))

    if WV is None:
        WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda

    h = H0/100.0
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1.0 - WM - WR - WV
    az = 1.0/(1.0 + 1.0*z)

    points = np.arange(n)

    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    a = az + (1 - az)*(points + 0.5)/n
    adot = np.sqrt(WK + (WM/a) + (WR/np.power(a, 2)) + (WV*np.power(a, 2)))
    DCMR = np.sum(1.0/(a*adot), axis=1)
    DCMR = ((1.0 - az)[:].T*DCMR[:])/n

    x = np.sqrt(np.absolute(WK))*DCMR
    x = x.reshape((x.size,))
    yes = np.where(x > 0.1)
    if WK > 0:
        tan_ratio_yes = np.sinh(x[yes])/x[yes]
    else:
        tan_ratio_yes = np.sin(x[yes])/x[yes]

    no = np.where(x <= 0.1)
    y = np.power(x[no], 2)
    if WK < 0:
        y = -y
    tan_ratio_no = 1.0 + y/6.0 + np.power(y, 2)/120.0

    tan_ratio = np.zeros(x.size)
    tan_ratio[yes] = tan_ratio_yes
    tan_ratio[no] = tan_ratio_no

    DCMT = tan_ratio*DCMR

    DA = az[:].T*DCMT[:]
    DL = (1/np.power(az, 2))[:].T*DA[:]

    DL_Mpc = c/H0*DL

    return DL_Mpc[0]


def main():
    opt, args, parser = get_options()

    if len(args) < 1:
        parser.error('Need a redshift')

    kwargs = dict((key, val) for (key, val) in opt.__dict__.items() if val is not None)
    z = float(args[0])
    cc = cosmocalc(z, **kwargs)
    try:
        outlines = []
        for outkey in (args[1:] or _outvals):
            outlines.append(outkey + ' = ' + str(cc[outkey]))
        print('\n'.join(outlines))
    except KeyError:
        parser.error(outkey + ' is not a valid output name_unit')

if __name__ == '__main__':
    main()

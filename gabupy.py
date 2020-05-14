import numpy as np

class Galaxy:
    """
    A class to store the necessary values of the different catalogs

    Parameters
    ----------
    num           : int
                    Unique integer number assigned in order for each instance
    name          : string
                    Mainly in PGC, but if there isn't in 2MASS
    dist          : float
                    Luminosity distance in Mpc
    RA            : float
                    Right ascension in deg
    dec           : float
                    Declination in deg
    v2MTF_K       : float
                    Value in the 2Mass Tully-Fisher catalog for K Band
    e2MTF_K       : float
                    Error in the measument of the 2MTF for K Band
    v2MTF_H       : float
                    Value in the 2Mass Tully-Fisher catalog for H Band
    e2MTF_H       : float
                    Error in the measument of the 2MTF for H Band
    v2MTF_J       : float
                    Value in the 2Mass Tully-Fisher catalog for J Band
    e2MTF_J       : float
                    Error in the measument of the 2MTF for J Band
    v2MTF_WF15_K  : float
                    Value in the 2Mass Tully-Fisher catalog for K Band using WF15 estimator
    e2MTF_WF15_K  : float
                    Error in the measument of the 2MTF for K Band using WF15 estimator
    v2MTF_WF15_H  : float
                    Value in the 2Mass Tully-Fisher catalog for H Band using WF15 estimator
    e2MTF_WF15_H  : float
                    Error in the measument of the 2MTF for H Band using WF15 estimator
    v2MTF_WF15_J  : float
                    Value in the 2Mass Tully-Fisher catalog for J Band using WF15 estimator
    e2MTF_WF15_J  : float
                    Error in the measument of the 2MTF for J Band using WF15 estimator
    vCF3          : float
                    Value in the Cosmicflows-3 catalog
    eCF3          : float
                    Error estimated for CF3 (100%)
    vCF3_mod      : float
                    Value in the Cosmicflows-3 catalog
    eCF3_mod      : float
                    Error estimated for CF3 (100%)
    vSFIpp        : float
                    Value in the SFI++ catalog
    eSFIpp        : float
                    Error in the measument of the SFI++
    v6dFGS        : float
                    Value in the 6dFGS catalog
    e6dFGS        : float
                    Error estimated for 6dFGS (100%)
    vRFGC_D       : float
                    Value in the Parnovsky RFGC catalog in the Dipole model
    eRFGC_D       : float
                    Error estimated for Parnovsky RFGC (100%)
    vRFGC_Q       : float
                    Value in the Parnovsky RFGC catalog in the Quadrupole model
    eRFGC_Q       : float
                    Error estimated for Parnovsky RFGC (100%)
    vRFGC_O       : float
                    Value in the Parnovsky RFGC catalog in the Octopole model
    eRFGC_O       : float
                    Error estimated for Parnovsky RFGC (100%)
    vBORG         : float
                    Value in the BORG model
    eBORG         : float
                    Error in the BORG model
    vCarrick      : float
                    Value in the Carrick model
    eCarrick      : float
                    Error estimated for Carrick (100%)
    vMLE          : float
                    Combined value of all the previous data
    eMLE          : float
                    Combined error of all the previous data
    """

    v2MTF_H = None
    e2MTF_H = None
    v2MTF_J = None
    e2MTF_J = None
    v2MTF_K = None
    e2MTF_K = None
    v2MTF_WF15_H = None
    e2MTF_WF15_H = None
    v2MTF_WF15_J = None
    e2MTF_WF15_J = None
    v2MTF_WF15_K = None
    e2MTF_WF15_K = None
    vCF3 = None
    eCF3 = None
    vCF3_mod = None
    eCF3_mod = None
    vSFIpp = None
    eSFIpp = None
    v6dFGS = None
    e6dFGS = None
    vRFGC_D = None
    eRFGC_D = None
    vRFGC_Q = None
    eRFGC_Q = None
    vRFGC_O = None
    eRFGC_O = None
    vBORG = None
    eBORG = None
    vCarrick = None
    eCarrick = None
    vMLE = None
    eMLE = None

    def __init__(self, num, name, dist, RA, dec):
        self.num = num
        self.name = name
        self.dist = dist
        self.RA = RA
        self.dec = dec

    def properties(self):
        """
        Print all the value store in galaxy instance

        It prints in multiple lines the value of each variable in the galaxy
        instance, in case there is no value it prints None

        TODO
        ----------
        it is printing an extra None
        """
        text = ("#: {num}\n" +
               "Name: {name}\n" +
               "dist: {dist}\n" +
               "RA = {RA}\n" +
               "dec = {dec}\n" +
               "v2MTF_H = {v2MTF_H}\n" +
               "e2MTF_H = {e2MTF_H}\n" +
               "v2MTF_J = {v2MTF_J}\n" +
               "e2MTF_J = {e2MTF_J}\n" +
               "v2MTF_K = {v2MTF_K}\n" +
               "e2MTF_K = {e2MTF_K}\n" +
               "v2MTF_WF15_H = {v2MTF_WF15_H}\n" +
               "e2MTF_WF15_H = {e2MTF_WF15_H}\n" +
               "v2MTF_WF15_J = {v2MTF_WF15_J}\n" +
               "e2MTF_WF15_J = {e2MTF_WF15_J}\n" +
               "v2MTF_WF15_K = {v2MTF_WF15_K}\n" +
               "e2MTF_WF15_K = {e2MTF_WF15_K}\n" +
               "vCF3 = {vCF3}\n" +
               "eCF3 = {eCF3}\n" +
               "vCF3_mod = {vCF3_mod}\n" +
               "eCF3_mod = {eCF3_mod}\n" +
               "vSFIpp = {vSFIpp}\n" +
               "eSFIpp = {eSFIpp}\n" +
               "v6dFGS = {v6dFGS}\n" +
               "e6dFGS = {e6dFGS}\n" +
               "vRFGC_D = {vRFGC_D}\n" +
               "eRFGC_D = {eRFGC_D}\n" +
               "vRFGC_Q = {vRFGC_Q}\n" +
               "eRFGC_Q = {eRFGC_Q}\n" +
               "vRFGC_O = {vRFGC_O}\n" +
               "eRFGC_O = {eRFGC_O}\n" +
               "vBORG = {vBORG}\n" +
               "eBORG = {eBORG}\n" +
               "vCarrick = {vCarrick}\n" +
               "eCarrick = {eCarrick}\n" +
               "vMLE = {vBORG}\n" +
               "eMLE = {eMLE}")

        print(text.format(num=self.num, name=self.name, dist=self.dist, RA=self.RA, dec=self.dec,
                          v2MTF_H=self.v2MTF_H, e2MTF_H=self.e2MTF_H, v2MTF_J=self.v2MTF_J, e2MTF_J=self.e2MTF_J,
                          v2MTF_K=self.v2MTF_K, e2MTF_K=self.e2MTF_K, v2MTF_WF15_H=self.v2MTF_H,
                          e2MTF_WF15_H=self.e2MTF_H, v2MTF_WF15_J=self.v2MTF_J, e2MTF_WF15_J=self.e2MTF_J,
                          v2MTF_WF15_K=self.v2MTF_K, e2MTF_WF15_K=self.e2MTF_K, vCF3=self.vCF3, eCF3=self.eCF3,
                          vCF3_mod=self.vCF3_mod, eCF3_mod=self.eCF3_mod, vSFIpp=self.vSFIpp, eSFIpp=self.eSFIpp,
                          v6dFGS=self.v6dFGS, e6dFGS=self.e6dFGS, vRFGC_D=self.vRFGC_D, eRFGC_D=self.eRFGC_D,
                          vRFGC_Q=self.vRFGC_Q, eRFGC_Q=self.eRFGC_Q, vRFGC_O=self.vRFGC_O, eRFGC_O=self.eRFGC_O,
                          vBORG=self.vBORG, eBORG=self.eBORG, vCarrick=self.vCarrick, eCarrick=self.eCarrick,
                          vMLE=self.vMLE, eMLE=self.vMLE))


class GLADE2:
    """
    A class to store in a more efficent way the GLADE Catalog

    PGC:        PGC number
    GWGC:       Name in the GWGC catalog
    HyperLEDA:  Name in the HyperLEDA catalog
    TwoMASS:    Name in the 2MASS XSC catalog
    SDSS-DR12:  Name in the SDSS-DR12 QSO catalog
    flag1:      Object type
                Q: the source is from the SDSS-DR12 QSO catalog
                C: the source is a globular cluster
                G: the source is from another catalog and not identified as a globular cluster
    RA:         Right ascention [deg]
    dec:        Declination [deg]
    dist:       Luminosity distance [Mpc]
    dist_err:   Error of distance [Mpc]
    z:          Redshift
    B:          Apparent B magnitude
    B_err:      Error of apparent B magnitude
    B_abs:      Absolute B magnitude
    J:          Apparent J magnitude
    J_err:      Error of apparent J magnitude
    H:          Apparent H magnitude
    H_err:      Error of apparent H magnitude
    K:          Apparent K magnitude
    K_err:      Error of apparent K magnitude
    flag2:      Luminosity distance measurement
                0: the galaxy had neither measured distance nor measured redshift value
                1: the galaxy had measured redshift value, from which we have calculated distance using the following
                   cosmological parameters: H_0=70 km s^-1 Mpc^-1 \Omega_M =0.27 and \Omega_{\Lambda}=0.73
                2: the galaxy had measured distance value from which we have calculated redshift using the following
                   cosmological parameters: H_0=70 km s^-1 Mpc^-1 \Omega_M =0.27 and \Omega_{\Lambda}=0.73
                3: The measured photometric redshift of the galaxy has been changed to spectroscopic redshift,
                   from which we have calculated distance using the following cosmological parameters:
                   H_0=70 km s^-1 Mpc^-1 \Omega_M =0.27 and \Omega_{\Lambda}=0.73
    flag3:      Peculiar velocity correction
                0: velocity field correction has not been applied to the object
                1: we have subtracted the radial velocity of the object
    """
    def __init__(self, PGC, GWGC, HyperLEDA, TwoMASS, SDSS_DR12, flag1, RA, dec, dist, dist_err, z, B, B_err, B_abs,
                 J, J_err, H, H_err, K, K_err, flag2, flag3):
        self.PGC = PGC
        self.GWGC = GWGC
        self.HyperLEDA = HyperLEDA
        self.TwoMASS = TwoMASS
        self.SDSS_DR12 = SDSS_DR12
        self.flag1 = flag1
        self.RA = RA
        self.dec = dec
        self.dist = dist
        self.dist_err = dist_err
        self.z = z
        self.B = B
        self.B_err = B_err
        self.B_abs = B_abs
        self.J = J
        self.J_err = J_err
        self.H = H
        self.H_err = H_err
        self.K = K
        self.K_err = K_err
        self.flag2 = flag2
        self.flag3 = flag3


class GLADE3:
    """
    A class to store in a more efficent way the GLADE Catalog

    PGC:        PGC number
    GWGC:       Name in the GWGC catalog
    HyperLEDA:  Name in the HyperLEDA catalog
    TwoMASS:    Name in the 2MASS XSC catalog
    SDSS-DR12:  Name in the SDSS-DR12 QSO catalog
    flag1:      Object type
                Q: the source is from the SDSS-DR12 QSO catalog
                C: the source is a globular cluster
                G: the source is from another catalog and not identified as a globular cluster
    RA:         Right ascention [deg]
    dec:        Declination [deg]
    dist:       Luminosity distance [Mpc]
    dist_err:   Error of distance [Mpc]
    z:          Redshift
    B:          Apparent B magnitude
    B_err:      Error of apparent B magnitude
    B_abs:      Absolute B magnitude
    J:          Apparent J magnitude
    J_err:      Error of apparent J magnitude
    H:          Apparent H magnitude
    H_err:      Error of apparent H magnitude
    K:          Apparent K magnitude
    K_err:      Error of apparent K magnitude
    flag2:      Luminosity distance measurement
                0: the galaxy had neither measured distance nor measured redshift value
                1: the galaxy had measured redshift value, from which we have calculated distance using the following
                   cosmological parameters: H_0=70 km s^-1 Mpc^-1 \Omega_M =0.27 and \Omega_{\Lambda}=0.73
                2: the galaxy had measured distance value from which we have calculated redshift using the following
                   cosmological parameters: H_0=70 km s^-1 Mpc^-1 \Omega_M =0.27 and \Omega_{\Lambda}=0.73
                3: The measured photometric redshift of the galaxy has been changed to spectroscopic redshift,
                   from which we have calculated distance using the following cosmological parameters:
                   H_0=70 km s^-1 Mpc^-1 \Omega_M =0.27 and \Omega_{\Lambda}=0.73
    flag3:      Peculiar velocity correction
                0: velocity field correction has not been applied to the object
                1: we have subtracted the radial velocity of the object
    wiseX:      IDK
    wiseID:     Name in the WISE catalog
    w1:         IDK
    w2:         IDK
    b_corr:     IDK
    r_corr:     IDK
    """
    def __init__(self, PGC, GWGC, HyperLEDA, TwoMASS, SDSS_DR12, flag1, RA, dec, dist, dist_err, z, B, B_err, B_abs,
                 J, J_err, H, H_err, K, K_err, flag2, flag3, wiseX, wiseID, w1, w2, b_corr, r_corr):
        self.PGC = PGC
        self.GWGC = GWGC
        self.HyperLEDA = HyperLEDA
        self.TwoMASS = TwoMASS
        self.SDSS_DR12 = SDSS_DR12
        self.flag1 = flag1
        self.RA = RA
        self.dec = dec
        self.dist = dist
        self.dist_err = dist_err
        self.z = z
        self.B = B
        self.B_err = B_err
        self.B_abs = B_abs
        self.J = J
        self.J_err = J_err
        self.H = H
        self.H_err = H_err
        self.K = K
        self.K_err = K_err
        self.flag2 = flag2
        self.flag3 = flag3
        self.wiseX = wiseX
        self.wiseID = wiseID
        self.w1 = w1
        self.w2 = w2
        self.b_corr = b_corr
        self.r_corr = r_corr


def interpolate(v000, v001, v010, v011, v100, v101, v110, v111, x, x0, x1, y, y0, y1, z, z0, z1):
    """
    TODO
    """
    avg_v_0 = v000*(x1 - x) + v001*(x - x0)
    avg_v_1 = v010*(x1 - x) + v011*(x - x0)
    avg_v_2 = v100*(x1 - x) + v101*(x - x0)
    avg_v_3 = v110*(x1 - x) + v111*(x - x0)

    avg_v_4 = avg_v_0*(y1 - y) + avg_v_1*(y - y0)
    avg_v_5 = avg_v_2*(y1 - y) + avg_v_3*(y - y0)

    avg_v = avg_v_4*(z1 - z) + avg_v_5*(z - z0)

    return(avg_v)


def mle(v, e):
    """
    Use the Maximum Likelihood Estimation to combine different values.
    """
    e = np.power(e, 2)
    avg_v = 0.0
    for i in range(e.size):
            avg_v += v[i]/e[i]

    avg_e = 0.0
    for i in range(e.size):
            avg_e += 1/e[i]

    avg_v = avg_v/avg_e
    avg_e = 1/np.sqrt(avg_e)

    return np.array([avg_v, avg_e])


def null2none(datum, type_datum):
    """
    To convert the string null in the txt file to None how it used in Python
    """
    if type_datum == 's':
        if datum == 'null':
            datum = None

    if type_datum == 'f':
        if datum == 'null':
            datum = None
        else:
            datum = np.float32(datum)

    if type_datum == 'i':
        if datum == 'null':
            datum = None
        else:
            datum = np.uint8(datum)

    return datum

def cosmocalc_DL_Mpc(z, H0=71, WM=0.27, WV=None, n=1000):
    c = 299792.458

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

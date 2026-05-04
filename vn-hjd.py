#!/usr/bin/python
########################################################################
# program: hjd.py
# author: Vitaly Neustroev
# version: 0.2
# date: May 12, 2021
# description:  
#    
#  Convert geocentric Julian date to Heliocentric Julian date.
#  The file must have two columns: time(JD, MJD, or truncated JD) & something else.
#
#  Program uses the local functions. Uncomment if want to use the PyAstronomy function
########################################################################

import sys
from numpy import *
from sys import stdin
#from PyAstronomy import pyasl  # uncomment if want to use the PyAstronomy function

#from astrolibpy import helio_jd

####################################################################################################

def bprecess(ra0, dec0, mu_radec=None, parallax=None, rad_vel=None, epoch=None):
    """
     NAME:
           BPRECESS
     PURPOSE:
           Precess positions from J2000.0 (FK5) to B1950.0 (FK4)
     EXPLANATION:
           Calculates the mean place of a star at B1950.0 on the FK4 system from
           the mean place at J2000.0 on the FK5 system.

     CALLING SEQUENCE:
           bprecess, ra, dec, ra_1950, dec_1950, [ MU_RADEC = , PARALLAX =
                                           RAD_VEL =, EPOCH =   ]

     INPUTS:
           RA,DEC - Input J2000 right ascension and declination in *degrees*.
                   Scalar or N element vector

     OUTPUTS:
           RA_1950, DEC_1950 - The corresponding B1950 right ascension and
                   declination in *degrees*.    Same number of elements as
                   RA,DEC but always double precision.

     OPTIONAL INPUT-OUTPUT KEYWORDS
           MU_RADEC - 2xN element double precision vector containing the proper
                      motion in seconds of arc per tropical *century* in right
                      ascension and declination.
           PARALLAX - N_element vector giving stellar parallax (seconds of arc)
           RAD_VEL  - N_element vector giving radial velocity in km/s

           The values of MU_RADEC, PARALLAX, and RADVEL will all be modified
           upon output to contain the values of these quantities in the
           B1950 system.  The parallax and radial velocity will have a very
           minor influence on the B1950 position.

           EPOCH - scalar giving epoch of original observations, default 2000.0d
               This keyword value is only used if the MU_RADEC keyword is not set.
     NOTES:
           The algorithm is taken from the Explanatory Supplement to the
           Astronomical Almanac 1992, page 186.
           Also see Aoki et al (1983), A&A, 128,263

           BPRECESS distinguishes between the following two cases:
           (1) The proper motion is known and non-zero
           (2) the proper motion is unknown or known to be exactly zero (i.e.
                   extragalactic radio sources).   In this case, the reverse of
                   the algorithm in Appendix 2 of Aoki et al. (1983) is used to
                   ensure that the output proper motion is  exactly zero. Better
                   precision can be achieved in this case by inputting the EPOCH
                   of the original observations.

           The error in using the IDL procedure PRECESS for converting between
           B1950 and J1950 can be up to 12", mainly in right ascension.   If
           better accuracy than this is needed then BPRECESS should be used.

           An unsystematic comparison of BPRECESS with the IPAC precession
           routine (http://nedwww.ipac.caltech.edu/forms/calculator.html) always
           gives differences less than 0.15".
     EXAMPLE:
           The SAO2000 catalogue gives the J2000 position and proper motion for
           the star HD 119288.   Find the B1950 position.

           RA(2000) = 13h 42m 12.740s      Dec(2000) = 8d 23' 17.69''
           Mu(RA) = -.0257 s/yr      Mu(Dec) = -.090 ''/yr

           IDL> mu_radec = 100D* [ -15D*.0257, -0.090 ]
           IDL> ra = ten(13, 42, 12.740)*15.D
           IDL> dec = ten(8, 23, 17.69)
           IDL> bprecess, ra, dec, ra1950, dec1950, mu_radec = mu_radec
           IDL> print, adstring(ra1950, dec1950,2)
                   ===> 13h 39m 44.526s    +08d 38' 28.63"

     REVISION HISTORY:
           Written,    W. Landsman                October, 1992
           Vectorized, W. Landsman                February, 1994
           Treat case where proper motion not known or exactly zero  November 1994
           Handling of arrays larger than 32767   Lars L. Christensen, march, 1995
           Converted to IDL V5.0   W. Landsman   September 1997
           Fixed bug where A term not initialized for vector input
                W. Landsman        February 2000
          Converted to python 			Sergey Koposov july 2010   
    """

    scal = True
    if isinstance(ra0, ndarray):
        ra = ra0
        dec = dec0
        n = ra.size
        scal = False
    else:
        n = 1
        ra = array([ra0])
        dec = array([dec0])

    if rad_vel is None:   
        rad_vel = zeros(n)
    else:
        if not isinstance(rad_vel, ndarray):
            rad_vel = array([rad_vel],dtype=float)
        if rad_vel.size != n:   
            raise Exception('ERROR - RAD_VEL keyword vector must be of the same length as RA and DEC')

    if (mu_radec is not None):   
        if (array(mu_radec).size != 2 * n):   
            raise Exception('ERROR - MU_RADEC keyword (proper motion) be dimensioned (2,' + strtrim(n, 2) + ')')
        mu_radec = mu_radec * 1.

    if parallax is None:   
        parallax = zeros(n)
    else:   
        if not isinstance(parallax, ndarray):
            parallax = array([parallax],dtype=float)

    if epoch is None:   
        epoch = 2000.0e0

    radeg = 180.e0 / pi
    sec_to_radian = lambda x : deg2rad(x/3600.)

    m = array([array([+0.9999256795e0, -0.0111814828e0, -0.0048590040e0, -0.000551e0, -0.238560e0, +0.435730e0]),
              array([+0.0111814828e0, +0.9999374849e0, -0.0000271557e0, +0.238509e0, -0.002667e0, -0.008541e0]),
   array([+0.0048590039e0, -0.0000271771e0, +0.9999881946e0, -0.435614e0, +0.012254e0, +0.002117e0]),
   array([-0.00000242389840e0, +0.00000002710544e0, +0.00000001177742e0, +0.99990432e0, -0.01118145e0, -0.00485852e0]),
   array([-0.00000002710544e0, -0.00000242392702e0, +0.00000000006585e0, +0.01118145e0, +0.99991613e0, -0.00002716e0]),
   array([-0.00000001177742e0, +0.00000000006585e0, -0.00000242404995e0, +0.00485852e0, -0.00002717e0, +0.99996684e0])])

    a_dot = 1e-3 * array([1.244e0, -1.579e0, -0.660e0])           #in arc seconds per century

    ra_rad = deg2rad(ra)
    dec_rad = deg2rad(dec)
    cosra = cos(ra_rad)
    sinra = sin(ra_rad)
    cosdec = cos(dec_rad)
    sindec = sin(dec_rad)

    dec_1950 = dec * 0.
    ra_1950 = ra * 0.

    for i in range(n):

    # Following statement moved inside loop in Feb 2000.
        a = 1e-6 * array([-1.62557e0, -0.31919e0, -0.13843e0])        #in radians

        r0 = array([cosra[i] * cosdec[i], sinra[i] * cosdec[i], sindec[i]])

        if (mu_radec is not None):   

            mu_a = mu_radec[i,0]
            mu_d = mu_radec[i,1]
            r0_dot = array([-mu_a * sinra[i] * cosdec[i] - mu_d * cosra[i] * sindec[i], mu_a * cosra[i] * cosdec[i] - mu_d * sinra[i] * sindec[i], mu_d * cosdec[i]]) + 21.095e0 * rad_vel[i] * parallax[i] * r0

        else:   
            r0_dot = array([0.0e0, 0.0e0, 0.0e0])

        r_0 = concatenate((r0, r0_dot))
        r_1 = transpose(dot(transpose(m), transpose(r_0)))

        # Include the effects of the E-terms of aberration to form r and r_dot.

        r1 = r_1[0:3]
        r1_dot = r_1[3:6]

        if mu_radec is None:   
            r1 = r1 + sec_to_radian ( r1_dot * (epoch - 1950.0e0) / 100. )
            a = a + sec_to_radian ( a_dot * (epoch - 1950.0e0) / 100. )

        x1 = r_1[0]   ;   y1 = r_1[1]    ;  z1 = r_1[2]
        rmag = sqrt(x1 ** 2 + y1 ** 2 + z1 ** 2)


        s1 = r1 / rmag    ; s1_dot = r1_dot / rmag

        s = s1
        for j in arange(0, 3):
            r = s1 + a - ((s * a).sum()) * s
            s = r / rmag
        x = r[0]          ; y = r[1]     ;  z = r[2]
        r2 = x ** 2 + y ** 2 + z ** 2
        rmag = sqrt(r2)

        if mu_radec is not None:   
            r_dot = s1_dot + a_dot - ((s * a_dot).sum()) * s
            x_dot = r_dot[0]  ; y_dot = r_dot[1]  ;  z_dot = r_dot[2]
            mu_radec[i,0] = (x * y_dot - y * x_dot) / (x ** 2 + y ** 2)
            mu_radec[i,1] = (z_dot * (x ** 2 + y ** 2) - z * (x * x_dot + y * y_dot)) / (r2 * sqrt(x ** 2 + y ** 2))

        dec_1950[i] = arcsin(z / rmag)
        ra_1950[i] = arctan2(y, x)

        if parallax[i] > 0.:   
            rad_vel[i] = (x * x_dot + y * y_dot + z * z_dot) / (21.095 * parallax[i] * rmag)
            parallax[i] = parallax[i] / rmag

    neg = (ra_1950 < 0)
    if neg.any() > 0:   
        ra_1950[neg] = ra_1950[neg] + 2.e0 * pi

    ra_1950 = rad2deg(ra_1950)
    dec_1950 = rad2deg(dec_1950)

    # Make output scalar if input was scalar
    if scal:
        return ra_1950[0],dec_1950[0]
    else:
        return ra_1950, dec_1950


####################################################################################################

def xyz(date, equinox=None):
    """
     NAME:
           XYZ
     PURPOSE:
           Calculate geocentric X,Y, and Z  and velocity coordinates of the Sun
     EXPLANATION:
           Calculates geocentric X,Y, and Z vectors and velocity coordinates
           (dx, dy and dz) of the Sun.   (The positive X axis is directed towards
           the equinox, the y-axis, towards the point on the equator at right
           ascension 6h, and the z axis toward the north pole of the equator).
           Typical position accuracy is <1e-4 AU (15000 km).

     CALLING SEQUENCE:
           XYZ, date, x, y, z, [ xvel, yvel, zvel, EQUINOX = ]

     INPUT:
           date: reduced julian date (=JD - 2400000), scalar or vector

     OUTPUT:
           x,y,z: scalars or vectors giving heliocentric rectangular coordinates
                     (in A.U) for each date supplied.    Note that sqrt(x^2 + y^2
                     + z^2) gives the Earth-Sun distance for the given date.
           xvel, yvel, zvel: velocity vectors corresponding to X, Y and Z.

     OPTIONAL KEYWORD INPUT:
           EQUINOX: equinox of output. Default is 1950.

     EXAMPLE:
           What were the rectangular coordinates and velocities of the Sun on
           Jan 22, 1999 0h UT (= JD 2451200.5) in J2000 coords? NOTE:
           Astronomical Almanac (AA) is in TDT, so add 64 seconds to
           UT to convert.

           IDL> xyz,51200.5+64.d/86400.d,x,y,z,xv,yv,zv,equinox = 2000

           Compare to Astronomical Almanac (1999 page C20)
                       X  (AU)        Y  (AU)     Z (AU)
           XYZ:      0.51456871   -0.76963263  -0.33376880
           AA:       0.51453130   -0.7697110   -0.3337152
           abs(err): 0.00003739    0.00007839   0.00005360
           abs(err)
               (km):   5609          11759         8040

           NOTE: Velocities in AA are for Earth/Moon barycenter
                 (a very minor offset) see AA 1999 page E3
                      X VEL (AU/DAY) YVEL (AU/DAY)   Z VEL (AU/DAY)
           XYZ:      -0.014947268   -0.0083148382    -0.0036068577
           AA:       -0.01494574    -0.00831185      -0.00360365
           abs(err):  0.000001583    0.0000029886     0.0000032077
           abs(err)
            (km/sec): 0.00265        0.00519          0.00557

     PROCEDURE CALLS:
           PRECESS_XYZ
     REVISION HISTORY
           Original algorithm from Almanac for Computers, Doggett et al. USNO 1978
           Adapted from the book Astronomical Photometry by A. Henden
           Written  W. Landsman   STX       June 1989
           Correct error in X coefficient   W. Landsman HSTX  January 1995
           Added velocities, more terms to positions and EQUINOX keyword,
              some minor adjustments to calculations
              P. Plait/ACC March 24, 1999
    """

    picon = pi / 180.0e0
    t = (date - 15020.0e0) / 36525.0e0         #Relative Julian century from 1900

    # NOTE: longitude arguments below are given in *equinox* of date.
    #   Precess these to equinox 1950 to give everything an even footing.
    #   Compute argument of precession from equinox of date back to 1950
    pp = (1.396041e0 + 0.000308e0 * (t + 0.5e0)) * (t - 0.499998e0)

    # Compute mean solar longitude, precessed back to 1950
    el = 279.696678e0 + 36000.76892e0 * t + 0.000303e0 * t * t - pp

    # Compute Mean longitude of the Moon
    c = 270.434164e0 + 480960.e0 * t + 307.883142e0 * t - 0.001133e0 * t * t - pp

    # Compute longitude of Moon's ascending node
    n = 259.183275e0 - 1800.e0 * t - 134.142008e0 * t + 0.002078e0 * t * t - pp

    # Compute mean solar anomaly
    g = 358.475833e0 + 35999.04975e0 * t - 0.00015e0 * t * t

    # Compute the mean jupiter anomaly
    j = 225.444651e0 + 2880.0e0 * t + 154.906654e0 * t * t

    # Compute mean anomaly of Venus
    v = 212.603219e0 + 58320.e0 * t + 197.803875e0 * t + 0.001286e0 * t * t

    # Compute mean anomaly of Mars
    m = 319.529425e0 + 19080.e0 * t + 59.8585e0 * t + 0.000181e0 * t * t

    # Convert degrees to radians for trig functions
    el = el * picon
    g = g * picon
    j = j * picon
    c = c * picon
    v = v * picon
    n = n * picon
    m = m * picon

    # Calculate X,Y,Z using trigonometric series
    x = 0.999860e0 * cos(el) - 0.025127e0 * cos(g - el) + 0.008374e0 * cos(g + el) + 0.000105e0 * cos(g + g + el) + 0.000063e0 * t * cos(g - el) + 0.000035e0 * cos(g + g - el) - 0.000026e0 * sin(g - el - j) - 0.000021e0 * t * cos(g + el) + 0.000018e0 * sin(2.e0 * g + el - 2.e0 * v) + 0.000017e0 * cos(c) - 0.000014e0 * cos(c - 2.e0 * el) + 0.000012e0 * cos(4.e0 * g + el - 8.e0 * m + 3.e0 * j) - 0.000012e0 * cos(4.e0 * g - el - 8.e0 * m + 3.e0 * j) - 0.000012e0 * cos(g + el - v) + 0.000011e0 * cos(2.e0 * g + el - 2.e0 * v) + 0.000011e0 * cos(2.e0 * g - el - 2.e0 * j)


    y = 0.917308e0 * sin(el) + 0.023053e0 * sin(g - el) + 0.007683e0 * sin(g + el) + 0.000097e0 * sin(g + g + el) - 0.000057e0 * t * sin(g - el) - 0.000032e0 * sin(g + g - el) - 0.000024e0 * cos(g - el - j) - 0.000019e0 * t * sin(g + el) - 0.000017e0 * cos(2.e0 * g + el - 2.e0 * v) + 0.000016e0 * sin(c) + 0.000013e0 * sin(c - 2.e0 * el) + 0.000011e0 * sin(4.e0 * g + el - 8.e0 * m + 3.e0 * j) + 0.000011e0 * sin(4.e0 * g - el - 8.e0 * m + 3.e0 * j) - 0.000011e0 * sin(g + el - v) + 0.000010e0 * sin(2.e0 * g + el - 2.e0 * v) - 0.000010e0 * sin(2.e0 * g - el - 2.e0 * j)


    z = 0.397825e0 * sin(el) + 0.009998e0 * sin(g - el) + 0.003332e0 * sin(g + el) + 0.000042e0 * sin(g + g + el) - 0.000025e0 * t * sin(g - el) - 0.000014e0 * sin(g + g - el) - 0.000010e0 * cos(g - el - j)

    #Precess_to new equator?
    if equinox is not None:   
        x, y, z = precess_xyz(x, y, z, 1950, equinox)

    xvel = -0.017200e0 * sin(el) - 0.000288e0 * sin(g + el) - 0.000005e0 * sin(2.e0 * g + el) - 0.000004e0 * sin(c) + 0.000003e0 * sin(c - 2.e0 * el) + 0.000001e0 * t * sin(g + el) - 0.000001e0 * sin(2.e0 * g - el)

    yvel = 0.015780 * cos(el) + 0.000264 * cos(g + el) + 0.000005 * cos(2.e0 * g + el) + 0.000004 * cos(c) + 0.000003 * cos(c - 2.e0 * el) - 0.000001 * t * cos(g + el)

    zvel = 0.006843 * cos(el) + 0.000115 * cos(g + el) + 0.000002 * cos(2.e0 * g + el) + 0.000002 * cos(c) + 0.000001 * cos(c - 2.e0 * el)

    #Precess to new equator?

    if equinox is not None:   
        xvel, yvel, zvel = precess_xyz(xvel, yvel, zvel, 1950, equinox)

    return x, y, z, xvel, yvel, zvel


####################################################################################################

def helio_jd(date, ra, dec, b1950=False, time_diff=False):
    """
     NAME:
          HELIO_JD
     PURPOSE:
          Convert geocentric (reduced) Julian date to heliocentric Julian date
     EXPLANATION:
          This procedure correct for the extra light travel time between the Earth
          and the Sun.

           An online calculator for this quantity is available at
           http://www.physics.sfasu.edu/astro/javascript/hjd.html
     CALLING SEQUENCE:
           jdhelio = HELIO_JD( date, ra, dec, /B1950, /TIME_DIFF)

     INPUTS
           date - reduced Julian date (= JD - 2400000), scalar or vector, MUST
                   be double precision
           ra,dec - scalars giving right ascension and declination in DEGREES
                   Equinox is J2000 unless the /B1950 keyword is set

     OUTPUTS:
           jdhelio - heliocentric reduced Julian date.  If /TIME_DIFF is set, then
                     HELIO_JD() instead returns the time difference in seconds
                     between the geocentric and heliocentric Julian date.

     OPTIONAL INPUT KEYWORDS
           /B1950 - if set, then input coordinates are assumed to be in equinox
                    B1950 coordinates.
           /TIME_DIFF - if set, then HELIO_JD() returns the time difference
                    (heliocentric JD - geocentric JD ) in seconds

     EXAMPLE:
           What is the heliocentric Julian date of an observation of V402 Cygni
           (J2000: RA = 20 9 7.8, Dec = 37 09 07) taken June 15, 1973 at 11:40 UT?

           IDL> juldate, [1973,6,15,11,40], jd      ;Get geocentric Julian date
           IDL> hjd = helio_jd( jd, ten(20,9,7.8)*15., ten(37,9,7) )

           ==> hjd = 41848.9881

     Wayne Warren (Raytheon ITSS) has compared the results of HELIO_JD with the
     FORTRAN subroutines in the STARLINK SLALIB library (see
     http://star-www.rl.ac.uk/).
                                                      Time Diff (sec)
          Date               RA(2000)   Dec(2000)  STARLINK      IDL

     1999-10-29T00:00:00.0  21 08 25.  -67 22 00.  -59.0        -59.0
     1999-10-29T00:00:00.0  02 56 33.4 +00 26 55.  474.1        474.1
     1940-12-11T06:55:00.0  07 34 41.9 -00 30 42.  366.3        370.2
     1992-02-29T03:15:56.2  12 56 27.4 +42 10 17.  350.8        350.9
     2000-03-01T10:26:31.8  14 28 36.7 -20 42 11.  243.7        243.7
     2100-02-26T09:18:24.2  08 26 51.7 +85 47 28.  104.0        108.8
     PROCEDURES CALLED:
           bprecess, xyz

     REVISION HISTORY:
           Algorithm from the book Astronomical Photometry by Henden, p. 114
           Written,   W. Landsman       STX     June, 1989
           Make J2000 default equinox, add B1950, /TIME_DIFF keywords, compute
           variation of the obliquity      W. Landsman   November 1999
           Converted to python 	Sergey Koposov July 2010
    """

    #Because XYZ uses default B1950 coordinates, we'll convert everything to B1950

    if not b1950:
        ra1, dec1 = bprecess(ra, dec)
    else:   
        ra1 = ra
        dec1 = dec


    delta_t = (array(date).astype(float) - 33282.42345905e0) / 36525.0e0
    epsilon_sec = poly1d([44.836e0, -46.8495, -0.00429, 0.00181][::-1])(delta_t)
    epsilon = deg2rad(23.433333e0 + epsilon_sec / 3600.0e0)
    ra1 = deg2rad(ra1)
    dec1 = deg2rad(dec1)

    x, y, z, tmp, tmp, tmp = xyz(date)

    #Find extra distance light must travel in AU, multiply by 1.49598e13 cm/AU,
    #and divide by the speed of light, and multiply by 86400 second/year

    time = -499.00522e0 * (cos(dec1) * cos(ra1) * x + (tan(epsilon) * sin(dec1) + cos(dec1) * sin(ra1)) * y)
    if time_diff:   
        return time
    else:   
        return array(date).astype(float) + time / 86400.0e0



##########################################################################

def WriteHJD(nn,aa,bb,output_file_path):
    """
    Write two columns of data to an external ASCII text file
    """   
    output_file = output_file_path.rstrip('\n') 
    outfile = open(output_file,"w")  
    for i in range (0, nn):
        outfile.write(' %15.7f \t %8.4f \n' %  (aa[i],bb[i]))
    outfile.close()             

#######################################

def enter_float(): 
    """
    Enter a floating point number and check its validity
    """   
    number_float = None        
    while not number_float:
        try:
            s =stdin.readline()
            number_float = float(s)
            if number_float == 0:
                break
        except ValueError:
            print('Invalid Number.  Enter number. ')
    return number_float
    
########################################################################


print("****************************************************")
print()
print("hjd.py inputfile outputfile RA(D.dd) Dec(D.dd) addJD")
print()
print("****************************************************")
print()
if (len(sys.argv) > 1):
  FileName = sys.argv[1]
else:
  print()
  print("Enter the input light curve filename: ")
  output_file_path = stdin.readline()
  FileName = output_file_path.rstrip('\n')

if (len(sys.argv) > 2):
    FileNameOutput = sys.argv[2]
else:
    print()
    print("Enter the output light curve filename: ")
    output_file_path = stdin.readline()
    FileNameOutput = output_file_path.rstrip('\n')

if (len(sys.argv) > 3):
    try:
      RA = float(sys.argv[3])
    except ValueError:
      print("Invalid RA. Enter Right ascension of object for epoch 2000.0 (degrees):")
      RA=enter_float()
else:
    print()
    print("Enter Right ascension of object for epoch 2000.0 (degrees):")
    RA=enter_float()

if (len(sys.argv) > 4):
    try:
      Dec = float(sys.argv[4])
    except ValueError:
      print("Invalid Dec. Enter Declination of object for epoch 2000.0 (degrees):")
      Dec=enter_float()
else:
    print()
    print("Enter Declination of object for epoch 2000.0 (degrees):")
    Dec=enter_float()

data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
Time=data0[0,:]
Flux=data0[1,:]
NumData = len(Time)
JD=Time[0]/2400000
MJD=Time[0]/50000

if (JD > 1) and (JD < 2):
    print("Time is assumed to be in JDs")
elif (MJD > 1) and (MJD < 2):
    print("Time is assumed to be in MJDs. To convert it to JDs, 2400000.5 will be added.")
    Time += 2400000.5
elif (Time[0] >= 0) and (Time[0] < 10000):
    print("Time is assumed to be in truncated JDs. To convert it to JDs, 2450000 will be added.")
    Time += 2450000.0
else:
    print("Unknown time format.")

print("After the conversion, the first time point is ", Time[0])

if (len(sys.argv) > 5):
    try:
      addJD = float(sys.argv[5])
    except ValueError:
      print("Invalid addJD. Enter a constant to be added to Times in order to have full JDs:")
      addJD=enter_float()
else:
    print()
    print("Enter a constant to be added to Times in order to have full JDs:")
    addJD=enter_float()

Time1 = zeros_like(Time)
for j in range(NumData):
    hjd= helio_jd(Time[j]+addJD-2400000.,RA,Dec)
    #hjd= pyasl.helio_jd(Time[j]+addJD-2400000.,RA,Dec)
    Time1[j]=hjd+2400000.

WriteHJD(NumData,Time1,Flux,FileNameOutput)
print("\nThe results are given in JDs.")

########################################################################


#! /usr/bin/env python3
# ver. 20210608
#from __future__ import print_function, division
from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt
import os
import sys
from sys import stdin 
from numpy import * 
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic



################################################################################################
def crosscorrRV(w, f, tw, tf, rvmin, rvmax, drv, mode="doppler", skipedge=0, edgeTapering=None):
    """
      Cross-correlate a spectrum with a template.

      The algorithm implemented here works as follows: For
      each RV shift to be considered, the wavelength axis
      of the template is shifted, either linearly or using
      a proper Doppler shift depending on the `mode`. The
      shifted template is then linearly interpolated at
      the wavelength points of the observation
      (spectrum) to calculate the cross-correlation function.
      VN: additionally, it also calculates chi2.

      Parameters
      ----------
      w : array
          The wavelength axis of the observation.
      f : array
          The flux axis of the observation.
      tw : array
          The wavelength axis of the template.
      tf : array
          The flux axis of the template.
      rvmin : float
          Minimum radial velocity for which to calculate
          the cross-correlation function [km/s].
      rvmax : float
          Maximum radial velocity for which to calculate
          the cross-correlation function [km/s].
      drv : float
          The width of the radial-velocity steps to be applied
          in the calculation of the cross-correlation
          function [km/s].
      mode : string, {lin, doppler}, optional
          The mode determines how the wavelength axis will be
          modified to represent a RV shift. If "lin" is specified,
          a mean wavelength shift will be calculated based on the
          mean wavelength of the observation. The wavelength axis
          will then be shifted by that amount. If "doppler" is
          specified (the default), the wavelength axis will
          properly be Doppler-shifted.
      skipedge : int, optional
          If larger zero, the specified number of bins will be
          skipped from the begin and end of the observation. This
          may be useful if the template does not provide sufficient
          coverage of the observation.
      edgeTapering : float or tuple of two floats
          If not None, the method will "taper off" the edges of the
          observed spectrum by multiplying with a sine function. If a float number
          is specified, this will define the width (in wavelength units)
          to be used for tapering on both sides. If different tapering
          widths shall be used, a tuple with two (positive) numbers
          must be given, specifying the width to be used on the low- and
          high wavelength end. If a nonzero 'skipedge' is given, it
          will be applied first. Edge tapering can help to avoid
          edge effects (see, e.g., Gullberg and Lindegren 2002, A&A 390).

      Returns
      -------
      dRV : array
          The RV axis of the cross-correlation function. The radial
          velocity refer to a shift of the template, i.e., positive
          values indicate that the template has been red-shifted and
          negative numbers indicate a blue-shift of the template.
          The numbers are given in km/s.
      CC : array
          The cross-correlation function.
    """
    if not _ic.check["scipy"]:
        raise(PE.PyARequiredImport("This routine needs scipy (.interpolate.interp1d).", \
                               where="crosscorrRV", \
                               solution="Install scipy"))
    import scipy.interpolate as sci
    # Copy and cut wavelength and flux arrays
    w, f = w.copy(), f.copy()
    if skipedge > 0:
        w, f = w[skipedge:-skipedge], f[skipedge:-skipedge]

    if edgeTapering is not None:
        # Smooth the edges using a sine
        if isinstance(edgeTapering, float):
            edgeTapering = [edgeTapering, edgeTapering]
        if len(edgeTapering) != 2:
            raise(PE.PyAValError("'edgeTapering' must be a float or a list of two floats.", \
                           where="crosscorrRV"))
        if edgeTapering[0] < 0.0 or edgeTapering[1] < 0.0:
            raise(PE.PyAValError("'edgeTapering' must be (a) number(s) >= 0.0.", \
                           where="crosscorrRV"))
        # Carry out edge tapering (left edge)
        indi = np.where(w < w[0]+edgeTapering[0])[0]
        f[indi] *= np.sin((w[indi] - w[0])/edgeTapering[0]*np.pi/2.0)
        # Carry out edge tapering (right edge)
        indi = np.where(w > (w[-1]-edgeTapering[1]))[0]
        f[indi] *= np.sin((w[indi] - w[indi[0]])/edgeTapering[1]*np.pi/2.0 + np.pi/2.0)

    # Speed of light in km/s
    c = 299792.458
    # Check order of rvmin and rvmax
    if rvmax <= rvmin:
        raise(PE.PyAValError("rvmin needs to be smaller than rvmax.",
                         where="crosscorrRV", \
                         solution="Change the order of the parameters."))
    # Check whether template is large enough
    if mode == "lin":
        meanWl = np.mean(w)
        dwlmax = meanWl * (rvmax/c)
        dwlmin = meanWl * (rvmin/c)
        if (tw[0] + dwlmax) > w[0]:
            raise(PE.PyAValError("The minimum wavelength is not covered by the template for all indicated RV shifts.", \
                           where="crosscorrRV", \
                           solution=["Provide a larger template", "Try to use skipedge"]))
        if (tw[-1] + dwlmin) < w[-1]:
            raise(PE.PyAValError("The maximum wavelength is not covered by the template for all indicated RV shifts.", \
                           where="crosscorrRV", \
                           solution=["Provide a larger template", "Try to use skipedge"]))
    elif mode == "doppler":
        # Ensure that the template covers the entire observation for all shifts
        maxwl = tw[-1] * (1.0+rvmin/c)
        minwl = tw[0] * (1.0+rvmax/c)
        if minwl > w[0]:
            raise(PE.PyAValError("The minimum wavelength is not covered by the template for all indicated RV shifts.", \
                           where="crosscorrRV", \
                           solution=["Provide a larger template", "Try to use skipedge"]))
        if maxwl < w[-1]:
            raise(PE.PyAValError("The maximum wavelength is not covered by the template for all indicated RV shifts.", \
                           where="crosscorrRV", \
                           solution=["Provide a larger template", "Try to use skipedge"]))
    else:
        raise(PE.PyAValError("Unknown mode: " + str(mode), \
                         where="crosscorrRV", \
                         solution="See documentation for available modes."))
    # Calculate the cross correlation
    drvs = np.arange(rvmin, rvmax, drv)
    cc = np.zeros(len(drvs))
    chi2 = np.zeros(len(drvs))
    for i, rv in enumerate(drvs):
        if mode == "lin":
            # Shift the template linearly
            fi = sci.interp1d(tw+meanWl*(rv/c), tf)
        elif mode == "doppler":
            # Apply the Doppler shift
            fi = sci.interp1d(tw*(1.0 + rv/c), tf)
        # Shifted template evaluated at location of spectrum
        cc[i] = np.sum((f-mean(f)) * (fi(w)-mean(fi(w)))) / (np.sum((f-mean(f))**2) * np.sum((fi(w)-mean(fi(w)))**2))**0.5 
        chi2[i] = np.sum((f-fi(w))**2)
    return drvs, cc, chi2

##########################################################################
def find_nearest(array, value):
    '''
    Find the nearest value inside an array.
    :param array: array
    :param value: desired value (float)
    :return: nearest value and its index
    '''

    idx = (abs(array - value)).argmin()
    return idx

##################################################################

def WriteData2(nn,aa,bb,output_file_path):
    """
    Write two columns of data to an external ASCII text file
    """
    output_file = output_file_path.rstrip('\n')
    outfile = open(output_file,"w")
    for i in range (0, nn):
        outfile.write(' %12.6f \t %12.6e \n' %  (aa[i],bb[i]))
#        outfile.write(' %12.6f \t %12.6f \n' %  (aa[i],bb[i]))
    outfile.close()

########################################################################

def WriteData3(nn,aa,bb,cc,output_file_path):
    """
    Write three columns of data to an external ASCII text file
    """
    outfile = open(output_file_path,"w")
    for i in range (0, nn):
        outfile.write(' %s \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        #outfile.write(' %12.6f \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close() 


#######################################################################################
def printHeader():
    """ Print the Header """
    print(" ")
    print("*******************************************************************************")
    print("                cross-corr: ver. 20210608  (c) Vitaly Neustroev")

def printUsage():
    """ Print the usage info """
    print(" ")
    print("Usage:")
    print("cross-cor.py Template ObsSpectrum RV- RV+ RVstep [HellCorr]")
    print(" ")
    print("Parameters:")
    print("  RV- and RV+: the RV range of the cross-correlation (in km/s).")
    print("  HelCorr: if given then corrections for heliocentric motion will be applied.")
    print(" ")
    print("*******************************************************************************")
    print(" ")

##########################################################################


printHeader()
cl=2.997e5
pi=3.1415926
isHelCorr = False


if (len(sys.argv) > 6):
    try:
        HelCorr=float(sys.argv[6])
        isHelCorr = True
    except:
        printUsage()
        exit()
else:
    print("No heliocentric correction will be applied.")
    #HelCorr = input("Enter the RV step of the cross-correlation: ")
    #HelCorr = float(HelCorr)

if (len(sys.argv) < 2):    
    printUsage()

if (len(sys.argv) > 1):
    TempFile = sys.argv[1]
else:
    print(" ")
    sys.stdout.write("Enter the filename of the template: ")
    sys.stdout.flush()
    output_file_path = stdin.readline()
    TempFile = output_file_path.rstrip('\n')    
if not os.path.exists(TempFile): 
    print("The file",TempFile, "doesn't exist. Stop executing!")
    exit()
else:
    data1 = loadtxt(TempFile, unpack=True, skiprows=0)
    WaveTemp = data1[0,:]
    FluxTemp = data1[1,:] #*1e17

if (len(sys.argv) > 2):
    SpecName = sys.argv[2]
else:
    print(" ")
    sys.stdout.write("Enter the filename of the spectrum: ")
    sys.stdout.flush()   
    output_file_path = stdin.readline()
    SpecName = output_file_path.rstrip('\n')
if not os.path.exists(SpecName): 
    print("The file",SpecName, "doesn't exist. Stop executing!")
    exit()
else:
    data2 = loadtxt(SpecName, unpack=True, skiprows=0)
    WaveSpec = data2[0,:]
    FluxSpec = data2[1,:] #*1e17
    try:
       SpecErr = data2[2,:]
       isErr = True
    except IndexError:
       isErr = False 

if isHelCorr:
    WaveSpec = cl*WaveSpec/(HelCorr+cl)
    WriteData2(len(WaveSpec),WaveSpec,FluxSpec,SpecName+'.helcorr')

TempNew = interp(WaveSpec,WaveTemp,FluxTemp)
TempAve1=average(TempNew)

idx_beg=find_nearest(WaveTemp,WaveSpec[0])
idx_end=find_nearest(WaveTemp,WaveSpec[-1])
TempAve=average(FluxTemp[idx_beg:idx_end])
SpecAve=average(FluxSpec)
FluxTemp=FluxTemp/TempAve1*SpecAve

if (len(sys.argv) > 3):
    try:
        RVneg=float(sys.argv[3])
    except:
        printUsage()
        exit()
else:
    RVneg = input("Enter the NEGATIVE RV range of the cross-correlation: ")
    RVneg = float(RVneg)

if (len(sys.argv) > 4):
    try:
        RVplus=float(sys.argv[4])
    except:
        printUsage()
        exit()
else:
    RVplus = input("Enter the POSITIVE RV range of the cross-correlation: ")
    RVplus = float(RVplus)

if (len(sys.argv) > 5):
    try:
        RVstep=float(sys.argv[5])
    except:
        printUsage()
        exit()
else:
    RVstep = input("Enter the RV step of the cross-correlation: ")
    RVstep = float(RVstep)


# Carry out the cross-correlation.
rv, cc, chi2 = crosscorrRV(WaveSpec, FluxSpec, WaveTemp, FluxTemp, RVneg, RVplus, RVstep, skipedge=0)

# Find the index of maximum cross-correlation function and minimum of chi2
maxind = np.argmax(cc)
minind = np.argmin(chi2)

print("Chi2 function is minimized at dRV = ", rv[minind], " km/s")
print("CCF  function is maximized at dRV = ", rv[maxind], " km/s")
print("CCF  function is = ", cc[maxind])
print()
if abs(rv[maxind]-rv[minind]) > 10.:
    print("[  Attention!  ]")
    print("dRV from CCF and Chi2 are different by ",rv[maxind]-rv[minind], "km/s")
    print()

if rv[maxind] > 0.0:
  print("A red-shift with respect to the template")
else:
  print("A blue-shift with respect to the template")

FluxSpecShifted, WaveSpecShifted = pyasl.dopplerShift(WaveSpec, FluxSpec, -rv[maxind], edgeHandling="fillValue",fillValue=1.0)

if isErr:
    WriteData3(len(WaveSpecShifted),WaveSpecShifted,FluxSpec,SpecErr,SpecName+'.shifted')   
else:
    WriteData2(len(WaveSpecShifted),WaveSpecShifted,FluxSpec,SpecName+'.shifted')

plt.figure(figsize=(12, 8))
plt.title("Cross-correlation is maximized at dRV = "+"%4.1f" % rv[maxind] +" km/s" )
plt.plot(rv, cc, 'bp-')
plt.plot(rv[maxind], cc[maxind], 'ro')
#plt.ylabel("Flux", size=14)
plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=14)          
plt.show()

# Plot template and data
plt.figure(figsize=(12, 8))
plt.title("Template (blue), original (blue) and shifted data (red)")
plt.plot(WaveSpecShifted, FluxSpec, 'r-')
plt.plot(WaveTemp, FluxTemp, 'k--')
plt.plot(WaveSpec, FluxSpec, 'b-')
plt.xlim(0.99*WaveSpecShifted[0],WaveSpecShifted[-1]*1.01)
plt.ylim(0.9*min(FluxSpec),max(FluxSpec)*1.1)
plt.ylabel("Flux", size=14)
plt.xlabel("Wavelength ($\mathrm{\AA}$)", size=14)

plt.show()

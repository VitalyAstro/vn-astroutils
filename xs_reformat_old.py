#!/usr/bin/python
#import numpy
from numpy import *
import sys
from sys import stdin
#import numpy as np
from pylab import *
#import pyfits
from astropy.io import fits as pyfits


##########################################################################
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
        outfile.write(' %12.6f \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close()
#########################################################################
def WriteData4(nn,aa,bb,cc,dd,output_file_path):
    """
    Write three columns of data to an external ASCII text file
    """
    outfile = open(output_file_path,"w")
    for i in range (0, nn):
        outfile.write(' %12.6f \t %12.6e \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i],dd[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close()
#########################################################################
def MedStd (x, k):
    """
    Median average
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = zeros ((len (x), k), dtype=x.dtype)
    y[:,k2] = x
    for i in range (k2):
        j = k2 - i
        y[j:,i] = x[:-j]
        y[:j,i] = x[0]
        y[:-j,-(i+1)] = x[j:]
        y[-j:,-(i+1)] = x[-1]
#    sort(a,axis=1)
    return std (y, axis=1, ddof=1)
#########################################################################
def MedFilt (x, k):
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = zeros ((len (x), k), dtype=x.dtype)
    y[:,k2] = x
    for i in range (k2):
        j = k2 - i
        y[j:,i] = x[:-j]
        y[:j,i] = x[0]
        y[:-j,-(i+1)] = x[j:]
        y[-j:,-(i+1)] = x[-1]
    # sort(a,axis=1)
    return median (y, axis=1)
#########################################################################
def Rebin2Col(X,Y,output_file_path):
    """
    Write three columns of data to an external ASCII text file
    """
    NumData = len(X)
    NumBins = NumData/Bins
    Xn=[]
    Yn=[]
    for i in range(NumBins):
        Xn.append(sum(X[i*Bins:(i+1)*Bins])/1.0/Bins)
        Yn.append(sum(Y[i*Bins:(i+1)*Bins])/1.0/Bins)
    WriteData2(len(Xn),Xn,Yn,output_file_path)
    return Xn,Yn
#########################################################################

def CleanSpec(FileName):
    """
    Write three columns of data to an external ASCII text file
    """
#    data0 = loadtxt(FileName, usecols=[0,1], unpack=True,skiprows=0)

    hdulist = pyfits.open(FileName)

    # get to the data part (in extension 1)
    scidata = hdulist[1].data

    Wavelength = scidata[0][0]
    Flux = scidata[0][1]
    arr2 = scidata[0][2]
    arr3 = scidata[0][3]
    arr4 = scidata[0][4]

    #data0 = loadtxt(FileName, usecols=[0,1,2], unpack=True,skiprows=535)
    data0 = loadtxt(FileName, usecols=[0,1,2], unpack=True,skiprows=0)
    Length=len(data0[0,:])
    EndSkip=2350
    EndSkip=0
    Wavelength=data0[0,:Length-EndSkip]
    Flux=data0[1,:Length-EndSkip]
    FluxNew=data0[1,:Length-EndSkip]
    FluxErr=data0[2,:Length-EndSkip]
    Length=len(Wavelength)

    FluxMed=MedFilt(Flux,Median1)
    F4std=sort(abs(Flux-FluxMed))
    AveDev=std(F4std[0:Length*0.9])
    Max=max(FluxMed)
    # FluxStd=MedStd(Flux,Median2)
    for i in range(Length):
        if abs(FluxNew[i]-FluxMed[i]) > Sigmas*AveDev:
            FluxNew[i] = FluxMed[i]
        if FluxNew[i] < 0.0:
            FluxNew[i] = FluxMed[i]
    if Clip0:
        clip(FluxNew,0.0, Max,out=FluxNew)

    if Plot:
        plot(Wavelength, Flux, "b-",Wavelength, FluxNew, "b-", Wavelength, FluxMed, "r-")
        xlabel('Wavelength')
        ylabel('Flux')
        show()

    WriteData4(Length,Wavelength,FluxNew,FluxErr,FluxMed,FileName+'.new')

    if Rebin:
        Xn,Yn = Rebin2Col(Wavelength,FluxNew,FileName+'.bin')
    else:
        Xn=Wavelength
        Yn=FluxNew

    if Interp:
        WaveInterp = linspace(InterW1, InterW2, InterPoints)
        FluxInterp = interp(WaveInterp, Xn,Yn)
#        FluxInterp = interp(WaveInterp, Wavelength,FluxNew)
        WriteData2(len(WaveInterp),WaveInterp,FluxInterp,FileName+'.interp')

#########################################################################


def ReadSpec(FileName):

    def coupling_flux_and_fluxerr_arrays():
    
        other_flux = []
        other_fluxerr = []
        matched_fluxerr = []
    
        # First, find all arrays that have to do with flux (using TUTYP)
        for i in other_arrays: # i is the FITS index (starts from 1, not from 0)
           if 'data.fluxaxis.value' in utype[i-1]:
               other_flux.append(i) 
           if 'data.fluxaxis.accuracy.staterror' in utype[i-1]:
               other_fluxerr.append(i) 
    
        print("Coupling flux and fluxerror arrays:")
        print("    - Coupled:           %s - %s" % (name[iflux-1], name[ierr-1]))
        
        for i in other_flux:
           i_namespace_length = utype[i-1].find(':')
           i_namespace = utype[i-1][:i_namespace_length]
           matched_couple = 0
           # for this flux, seek the corresponding fluxerr
           for j in other_fluxerr:
               j_namespace_length = utype[j-1].find(':')
               j_namespace = utype[j-1][:j_namespace_length]
               if i_namespace == j_namespace:
                   print("    - Coupled:           %s - %s - namespace: %s " % (name[i-1], name[j-1], i_namespace))
                   matched_couple=1
                   matched_fluxerr.append(j)
           if matched_couple == 0:
               # no matching fluxerr could be found, reporting:
               print("    - Uncoupled flux:    %s - utype: %s" % (name[i-1], Utype[i-1]))
    
        # Above the unmatched flux arrays were reported;
        # here we report the unmatched fluxerr arrays (if any).
    
        unmatched_fluxerr = list(set(other_fluxerr) - set(matched_fluxerr))
        for i in unmatched_fluxerr:
            print("    - Uncoupled fluxerr: %s - utype: %s" % (name[i-1], Utype[i-1]))
    
        print("")
    
    ######### End of coupling_flux_and_fluxerr_arrays()  #############################################


# Open file, get pointers to the primary header unit, the header of the first extension, the data part:
        
    hdulist = pyfits.open(FileName)   # Open FITS file
    
    phu    = hdulist[0].header        # Primary Header Unit: metadata
    scihead = hdulist[1].header       # Header of the first FITS extension: metadata
    scidata = hdulist[1].data         # Data in the first FITS extension: the spectrum
    
    # Checking some compliance
    # 1. keyword PRODCATG must be present in primary header unit
    try:
        prodcatg = phu['PRODCATG']
    except:
        errorMessage = 'Keyword PRODCATG not found in primary header.\nFile not compliant with the ESO Science Data Product standard.'
        print('FileName = %s   NOT COMPLIANT' % filename)
        print(errorMessage)
#        exit(1)
    
    # 2. value of keyword PRODCATG must match SCIENCE.SPECTRUM*
    if not prodcatg.startswith('SCIENCE.SPECTRUM'):
        errorMessage = "Expected header keyword: PRODCATG = 'SCIENCE.SPECTRUM'\nFound: PRODCATG = '%s'\nFile not compliant with the 1d spectrum specifications\nof the ESO Science Data Product standard." % prodcatg
        print('FileName = %s   NOT COMPLIANT' % filename)
        print(errorMessage)
        exit(1)
    
    # 3. Various keywords must be defined, among them the ones here below:
    try:
        #origfile=phu['ORIGFILE']    # Original filename as assigned by the data provider
        origfile = phu['OBJECT']    # Object name (the python code corrected by V. Neustroev)
        MJD =  scihead['TMID']          # MJD time of observation (the python code corrected by V. Neustroev)
        DISPELEM = phu['DISPELEM']  # ARM (UVB/Vis/NIR) (the python code corrected by V. Neustroev)
        instrume = phu['INSTRUME']  # Name of the instrument
        wavelmin = phu['WAVELMIN']  # Minimum wavelength in nm
        wavelmax = phu['WAVELMAX']  # Maximum wavelength in nm
        respower = phu['SPEC_RES']  # Spectral resolving power (lambda / delta_lambda)
        snr      = phu['SNR']       # Signal to Noise Ratio
        specaxisucd  = scihead['TUCD1'] # Gives the type of spectral axis (see SPECTRAL AXIS below)
    except:
       errorMessage='File not compliant with the 1D spectrum specifications of the ESO Science Data Product standard; some of the mandatory keywords were not found in primary header unit'
       print('FileName = %s   NOT COMPLIANT' % filename)
       print('ERROR = %s' % (errorMessage))
       exit(1)
    
    # SPECTRAL AXIS: could be either wavelength, frequency, or energy;
    # if wavelength, the distinction between wavelength in air or in vacuum is provided by the presence of the obs.atmos token in the TUCD1.
    # the variable spectype will carry to whole info.
    spectype = None
    if specaxisucd.startswith('em.wl'):
        if specaxisucd == 'em.wl':
            spectype = 'wavelength in vacuum (TUCD1=%s)' % specaxisucd
        elif specaxisucd == 'em.wl;obs.atmos':
            spectype = 'wavelength in air (TUCD1=%s)' % specaxisucd
        else:
            spectype = 'wavelength (TUCD1=%s)' % specaxisucd
    elif specaxisucd.startswith('em.freq'):
        spectype = 'frequency (TUCD1=%s)' % specaxisucd
    elif specaxisucd.startswith('em.ener'):
        spectype = 'energy (TUCD1=%s)' % specaxisucd
    
    # Report main characteristics of the spectrum:
    #print('************************************************************************************************************************')
    print('filename=%s   ORIGFILE=%s'  % (FileName,origfile))
    print('Instrume=%s   Wmin=%snm   Wmax=%snm   R=%s   SNR=%s'  % (instrume,wavelmin,wavelmax,respower,snr))
    print('Spectral axis: %s' % (spectype))
    print('------------------------------------------------------------------------------------------------------------------------')
    
    # Check VO compliance (ESO SDP is based on the VO standard):
    try:
        voclass=scihead['VOCLASS']
    except:
        print('File %s is not a valid VO 1d spectrum (VOCLASS keyword missing)' % (FileName))
        exit(1)
    
    # TFIELDS is a required FITS binary table keyword
    try:
        tfields=int(scihead['TFIELDS'])
    except:
        print('File %s is not a valid ESO SDP 1d spectrum (TFIELDS keyword missing)' % (filename))
        exit(1)
    
    #################################
    # METADATA PART
    #################################
    
    # Reading name, unit, utype for each column (array) in the FITS binary table (extension 1).
    
    name = []
    unit = []
    utype= [] # lowercase utype string: for case-insensitive matches
    Utype= [] # original utype, with case preserved: for display
    
    print("AVAILABLE ARRAYS:")
    print ("name            index  UNIT                               UTYPE")
    for i in range(1, tfields+1):
        thisname = scihead['TTYPE'+str(i)]
        try:
           thisunit = scihead['TUNIT'+str(i)]
        except:
           thisunit=""
        try:
           thisutype=scihead['TUTYP'+str(i)]
        except:
           thisutype='no_utype_assigned:field_not_part_of_the_standard'
        print ("%-15s %2d     %-34s [%-s]" % (thisname, i, thisunit, thisutype))
        name.append(thisname)
        unit.append(thisunit)
        utype.append(thisutype.lower())
        Utype.append(thisutype)
    
    print('------------------------------------------------------------------------------------------------------------------------')
    
    # Recognising the main scientific arrays (spectral, flux and flux error) and the "other" ones.
    # A 1D spectrum can contain several flux (and fluxerror) arrays, but one is defined to be the best.
    # The best one can be recognised by the (lowercased) utype which is either "spectrum.data.fluxaxis.value" or "spec:data.fluxaxis.value".
    
    other_arrays = []  # It will contain the indeces of the fields not considered main arrays. FITS indeces starts from 1!
    
    # Getting the indexes of the FITS columns
    # for the main spectral array (ispec), flux array (iflux), and flux_error (ierr) array:
    
    for i in range(1, tfields+1):
    
         # Remember that the index of Python arrays starts from 0, while the FITS index from 1.
         tutyp=utype[i-1]
    
         # The ESO Science Data Product standard format
         # prescribes the spectral axis to be stored in column 1;
         # there would be no need to look for it, but we need the other_arrays anyway.
    
         # The TUTYPn keywords follow either the Spectrum Data Model standard v1.1 for spectra with a single flux array,
         # or the Spectral Data Model standard v2.0 for spectra with any number of flux arrays
         # These data model standards are available from the International Virtual Observatory Alliance
         # web site at: http://ivoa.net/documents/
    
         if tutyp == 'spectrum.data.spectralaxis.value':
             ispec = i
         elif tutyp == 'spec:data.spectralaxis.value':
             ispec = i
         elif tutyp == 'spectrum.data.fluxaxis.value':
             iflux = i
         elif tutyp == 'spec:data.fluxaxis.value':
             iflux = i
         elif tutyp == 'spectrum.data.fluxaxis.accuracy.staterror':
             ierr  = i
         elif tutyp == 'spec:data.fluxaxis.accuracy.staterror':
             ierr  = i
         #elif tutyp == 'spec:Data.FluxAxis.Accuracy.QualityStatus':
             #iqual  = i
         #elif tutyp == 'spec:Data.FluxAxis.Accuracy.QualityStatus':
             #iqual  = i
         else:
             # Storing the indeces of other, not considered main, arrays:
             other_arrays.append( i )
    
    # --------------------------------------------------------------------------------------------------------
    # Checking if other flux and fluxerr arrays exist, and coupling them by the namespace in the utype
    # E.g.: eso:Data.FluxAxis.Value and eso:Data.FluxAxis.Accuracy.StatError form a couple, in the sense
    # that the second array is the error related to the flux stored in the first array.
    
    coupling_flux_and_fluxerr_arrays()
    
    # --------------------------------------------------------------------------------------------------------
    
    
    # Number of points in the spectrum: NELEM
    NELEM = scihead['NELEM']
     
    
    #################################
    # DATA PART and plots
    #################################
    
    print("\nThe spectrum has %d points\n" % (NELEM))
    
    # Main arrays:
    spec = np.array(scidata[0][ispec - 1])
    flux = np.array(scidata[0][iflux - 1])
    err  = np.array(scidata[0][ierr - 1])
    user_index = 4
    qual = np.array(scidata[0][user_index - 1])
    NewFileName = 'MJD'+ str(round(MJD,7)) + '_' + DISPELEM + '.txt'
    WriteData4(len(spec),spec*10,flux,err,qual,NewFileName)


########################################################################################
#                               MAIN
########################################################################################
        

print "*****************************************"
print " "
print "xs_reformat.py FileList"
print " "
print "*****************************************"
print " "

if (len(sys.argv) > 1):
  FileList = sys.argv[1]
else:
  print " "
  print "Enter the filename of the list of spectra: "
  output_file_path = stdin.readline()
  FileList = output_file_path.rstrip('\n')


Median1 = 21
Median2 = 21
Sigmas = 5.0
Clip0 = 1
#Clip0 = False
Interp = 1
Interp = False
Rebin = 1
#Rebin = False
Plot = 1
Plot = False

if Rebin:
    Bins = int(raw_input("Enter the binning factor (integer number of points): "))
if Interp:
    InterPoints = int(raw_input("Enter the integer number of points in the interpolated spectrum: "))
    InterW1 = int(raw_input("Enter the first wavelength in the interpolated spectrum: "))
    InterW2 = int(raw_input("Enter the last wavelength in the interpolated spectrum: "))

infile = open(FileList, "r")
lines = infile.readlines()
infile.close()

#Test
# FileName=FileList
# data0 = loadtxt(FileName, usecols=[0,1,2], unpack=True,skiprows=0)
# Wavelength=data0[0,:]
# Length=len(Wavelength)
# Flux=data0[1,:]
# FluxNew=data0[1,:]
# FluxErr=data0[2,:]
#
# FluxMed=MedFilt(Flux,Median1)
# F4std=sort(abs(Flux-FluxMed))
# AveDev=std(F4std[0:Length*0.9])
# print "Std=",AveDev
# Max=max(FluxMed)
# print "MaxMed=", Max
# # FluxStd=MedStd(Flux,Median2)
# for i in range(Length):
#     if abs(FluxNew[i]-FluxMed[i]) > Sigmas*AveDev:
#         FluxNew[i] = FluxMed[i]
#     if FluxNew[i] < 0.0:
#         FluxNew[i] = FluxMed[i]
# if Clip0:
#     clip(FluxNew,0.0, Max,out=FluxNew)
#
# plot(Wavelength, Flux, "b-",Wavelength, FluxNew, "b-", Wavelength, FluxMed, "r-")
# xlabel('Wavelength')
# ylabel('Flux')
# show()
#
# WriteData4(Length,Wavelength,FluxNew,FluxErr,FluxMed,FileName+'.new')
#Test


for line in lines:
    FileName = line.rstrip('\n')
    ReadSpec(FileName)
    #CleanSpec(FileName)
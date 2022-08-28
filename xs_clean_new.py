#!/usr/bin/python
#import numpy

from numpy import *
import sys
from sys import stdin
#import numpy as np
from pylab import *
from scipy import signal
from scipy import interpolate
import os
import os.path



# ==============================================================================
def find_nearest(array, value):
    '''
    Find the nearest value inside an array.
    :param array: array
    :param value: desired value (float)
    :return: nearest value and its index
    '''

    idx = (np.abs(array - value)).argmin()

    return idx


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

def Rebin(X,Y,Yerr,isErr,Bins):
    """
    Rebinning the XYYerr dataset merging N_bins together.
    """
    NumData = len(X)
    NumBins = int(NumData/Bins)
    Xn=[]
    Yn=[]
    Err=[]    
    if isErr:
        Weights = 1./Yerr**2
    for i in range(NumBins):
        Xn.append(np.mean(X[i*Bins:(i+1)*Bins]))
        if isErr:
            Ytemp,Etemp = weighted_avg_and_std(Y[i*Bins:(i+1)*Bins],Weights[i*Bins:(i+1)*Bins])
            Yn.append(Ytemp)
            Err.append(Etemp)
        else:
            Yn.append(np.mean(Y[i*Bins:(i+1)*Bins]))
            Err.append(np.std(Y[i*Bins:(i+1)*Bins]))
    return Xn,Yn,Err

#########################################################################

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = 1. / np.sum(weights)
    return (average, math.sqrt(variance))

#########################################################################


def SpecCut(aa,bb,cc,W1,W2):
    """
    A spectrum cut
    """
    idx_beg=find_nearest(aa,W1)
    idx_end=find_nearest(aa,W2)   
    wave = []
    flux = []    
    fluxerr = []
    wave.extend(aa[idx_beg:idx_end+1])
    flux.extend(bb[idx_beg:idx_end+1])
    fluxerr.extend(cc[idx_beg:idx_end+1])
    
    return wave,flux,fluxerr

#########################################################################

def FirstSpec(FileName,FileNameSky):
    """
    Reading of the wavelengths from the first spectrum and Testing for the presense of columns.
    """
    global isClipMask, isClean, isWeighted, isSky, Sky, dLam, SkyStrengthLimit
    try:
        data0 = loadtxt(FileName, unpack=True,skiprows=0)
        Wavelength=data0[0,:]
        Length=len(Wavelength)
        Flux=data0[1,:]
    except:
        print("Something wrong with the file ",FileName)
        print("Please check whether it is a spectrum or not. Exiting...")
        exit(-1)
    try:
       FluxErr=data0[2,:]
       isErr = True
       print ("The first spectrum has the 3rd column. We assume that ALL the files include error data.")
    except IndexError:
       isErr = False
       isWeighted = False
       isClean = False
       #SpecErr = data1[1,:]
    try:
       Mask=data0[3,:]
       isMask = 1
       print ("The first spectrum has the 4th column. We assume that ALL the files include mask data.")   
    except IndexError:
       isMask = 0
       isClipMask = False
       print(isClipMask)

    if isSky:
        try:
            data1 = loadtxt(FileNameSky, unpack=True,skiprows=0)
            WaveSky=data1[0,:]
            StrengthSky=data1[1,:]
            Sky = np.zeros_like(Wavelength)
            for i in range(len(WaveSky)):
                # idx_beg=find_nearest(Wavelength,WaveSky[i]-dLam)
                # idx_end=find_nearest(Wavelength,WaveSky[i]+dLam)
                # Sky[idx_beg:idx_end+1] = 2048       
                if (StrengthSky[i] <= SkyStrengthLimit):
                    #print(WaveSky[i],Wavelength[0])
                    continue
                else:
                    idx_beg=find_nearest(Wavelength,WaveSky[i]-dLam)
                    idx_end=find_nearest(Wavelength,WaveSky[i]+dLam)
                    Sky[idx_beg:idx_end+1] = 2048            
            Sky[0]  = 0
            Sky[-1] = 0                
        except FileNotFoundError:
            print ("The file with sky lines is not found. No sky lines will be cut off.")
            isSky = False
    return Wavelength,isErr,isMask


#########################################################################

def CleanSpectra(FileNames,FirstWave,isErr,isMask):
    """
    Clean spectra and write them to text files
    """
    
    for line in FileNames:
        FileName = line.rstrip('\n')
        try:
            data0 = loadtxt(FileName, unpack=True,skiprows=0)
            Wavelength=data0[0,:]
            Length=len(Wavelength)
            Flux=data0[1,:]
            print("\n",FileName, end = '')
        except:
            print("Something wrong with the file ",FileName)
            print("Please check whether it is a spectrum or not. Skipping...")
            continue        
        WaveMask=[]
        FluxMask=[]
        FluxErrMask=[]        
        WaveClip=[]
        FluxClip=[]
        FluxErrClip=[]        
        if isErr: 
            FluxErr=data0[2,:]
            if isMask and isClipMask:
                Mask=data0[3,:]
                for i in range(Length):
                    if Mask[i] == 0:
                        WaveMask.append(Wavelength[i])
                        FluxMask.append(Flux[i])
                        FluxErrMask.append(FluxErr[i])
                FluxMaskNew = interp(FirstWave,WaveMask,FluxMask)
                FluxErrMaskNew = interp(FirstWave,WaveMask,FluxErrMask)                                
                if not isInterp:
                    WriteData3(len(FirstWave),FirstWave,FluxMaskNew,FluxErrMaskNew,FileName+'.mask')
                    print(' --> .mask', end = '')
                else:
                    WaveCut,FluxCut,ErrCut = SpecCut(FirstWave,FluxMaskNew,FluxErrMaskNew,InterW1,InterW2)
                    WriteData3(len(WaveCut),WaveCut,FluxCut,ErrCut,FileName+'.crop-mask')                
                    print(' --> .crop-mask', end = '')
            else:
                FluxMaskNew = interp(FirstWave,Wavelength,Flux)
                FluxErrMaskNew = interp(FirstWave,Wavelength,FluxErr)
            if isClean:
                FluxMed=signal.medfilt(FluxMaskNew,Median1)
                for i in range(len(FirstWave)):
                    if abs(FluxMaskNew[i]-FluxMed[i]) <= Sigmas*FluxErrMaskNew[i]:
                        WaveClip.append(FirstWave[i])
                        FluxClip.append(FluxMaskNew[i])
                        FluxErrClip.append(FluxErrMaskNew[i])
                FluxClipNew = interp(FirstWave,WaveClip,FluxClip)
                FluxErrClipNew = interp(FirstWave,WaveClip,FluxErrClip)
                if not isInterp:               
                    WriteData3(len(FirstWave),FirstWave,FluxClipNew,FluxErrClipNew,FileName+'.clean')
                    print(' --> .clean', end = '')
                else:
                    WaveCut,FluxCut,ErrCut = SpecCut(FirstWave,FluxClipNew,FluxErrClipNew,InterW1,InterW2)
                    WriteData3(len(WaveCut),WaveCut,FluxCut,ErrCut,FileName+'.crop-clean')                
                    print(' --> .crop-clean', end = '')
            else:
                FluxClipNew = FluxMaskNew
                FluxErrClipNew = FluxErrMaskNew
        else:
            FluxClipNew = Flux
            FluxErrClipNew = np.ones_like(Flux)
            
        if isSky:
            WaveClip=[]
            FluxClip=[]
            FluxErrClip=[]
            for i in range(Length):
                if Sky[i] == 0:
                    WaveClip.append(Wavelength[i])
                    FluxClip.append(FluxClipNew[i])
                    FluxErrClip.append(FluxErrClipNew[i])
            f1 = interpolate.interp1d(WaveClip,FluxClip)
            f2 = interpolate.interp1d(WaveClip,FluxErrClip)
            FluxClipNew = f1(FirstWave)
            FluxErrClipNew = f2(FirstWave)
            if isInterp:
                WaveCut,FluxCut,ErrCut = SpecCut(FirstWave,FluxClipNew,FluxErrClipNew,InterW1,InterW2)
                if isErr:
                    WriteData3(len(WaveCut),WaveCut,FluxCut,ErrCut,FileName+'.crop-sky')
                else:
                    WriteData2(len(WaveCut),WaveCut,FluxCut,FileName+'.crop-sky')
                print(' --> .crop-sky', end = '')
            else:
                if isErr:
                    WriteData3(len(FirstWave),FirstWave,FluxClipNew,FluxErrClipNew,FileName+'.sky')
                else:
                    WriteData2(len(FirstWave),FirstWave,FluxClipNew,FileName+'.sky')
                print(' --> .sky', end = '')

        if isBin:
            if isInterp:
                WaveBin,FluxBin,ErrBin = Rebin(array(WaveCut),array(FluxCut),array(ErrCut),isErr,Bins)
            else:
                WaveBin,FluxBin,ErrBin = Rebin(FirstWave,FluxClipNew,FluxErrClipNew,isErr,Bins)
            WriteData3(len(WaveBin),WaveBin,FluxBin,ErrBin,FileName+'.bin')
            print(' --> .bin', end = '')       
    return 

#########################################################################



def print_header():
    print ("")
    print ("********************************************************************************")
    print ("**                              xs_clean.py                                   **")
    print ("** A little utility to clean TXT-files of spectra obtained with ESO/X-shooter **")
    print ("**                    from bad data points and outliers                       **")
    print ("**                             2021-June-07                                    **")
    print ("**                           Vitaly Neustroev                                 **")
    print ("********************************************************************************")
    print ("")

def usage():
    print_header()
    "Usage function"
    print ("Usage: %s [options] FileName/@FileList" % sys.argv[0])
    print (" ")
    print ("FileName is a spectrum filename; @FileList is a list of spectra.")
    print ("Each spectrum can consist of 4 columns: Wavelength, Flux, FluxErr QualityMask")
    print ("                           (for example, spectra obtained with ESO/X-shooter)")
    print ("Options: -cqslmpb")
    print ("     -h: Help")
    print ("     -c: The spectra will NOT be cleaned    [default: will be cleaned]")
    print ("     -q: The spectra will NOT be clipped    [default: will be clipped using the QualityMask]")
    print ("     -s: The spectra will NOT be clipped off the sky lines [default: will be clipped using") 
    print ("                                 the wavelengths of skylines from the file 'skylines.txt']")
    print ("     +s: Will be asked to enter the Strength limit of the skylines and the range of wavelengths") 
    print ("                                          around skylines to be cut off [default: 0.1 and 1.25]")
    print ("     -l: Will be asked to enter the number of standard deviations at the given wavelength") 
    print ("                                        to use for the clipping limit [default value is 5]")
    print ("     -m: Will be asked to enter the size of the median filter window [default value is 11]")
    print ("     -p: Will be asked to enter the first and the last wavelengths in cropped spectra")
    print ("            (the cropping will be applied after the cleaning and clipping procedures)")
    print ("     -b: Will be asked to enter the the binning factor (integer number of points)")
    print ("                             (the rebinning will be applied after all procedures)")
    print ("")
    sys.exit(-1)


##########################################################################

#os.chdir("/home/benj/OwnCloud/My_Papers/BW_Scl/Data/Xsh/NIR/NewReduction/Telluric_Cleaned/Phased/Test/")
#os.chdir("/home/benj/OwnCloud/Scripts/T/")

global isClipMask, isClean, isSky, isInterp, isBin, isErr, Sky, dLam, SkyStrengthLimit
FileList=False
isClipMask = True
isClean = True
isInterp = False
isSky = True
isBin = False
lines = []
Sky = []
HomeDir = "/scisoft/Other_Soft/Files4scripts/"
FileNameSky = HomeDir+"skylines.txt"

dLam = 1.25
SkyStrengthLimit = 0.1
Median1 = 11
Sigmas = 5.0

s1 = 0
spec_txt = "spectra"

if len(sys.argv) == 1:
    usage()

for i in range(len(sys.argv)-1):
    CmdLinePar = sys.argv[i+1]
    if (CmdLinePar[0] != '-') and (CmdLinePar[0] != '+'):
        if not FileList:
            FileName = CmdLinePar
            if FileName[0] == '@':
                s1 = 1
                FileList = True
            if not os.path.isfile(FileName[s1:]):
                print("The File ",FileName," doesn't exist. Exiting...")
                exit(-1)               
            try:
                if FileList:
                    infile = open(FileName[s1:], "r")
                    lines = infile.readlines()
                    infile.close()
                    print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
                else:
                    lines.append(FileName[s1:])
                    FileList = True
                    spec_txt = "spectrum"
            except:
                print("Something wrong with the FileName/FileList ",FileName[s1:])               

    elif CmdLinePar[0] == '-':
        if ('h' or 'H') in CmdLinePar[1:]:
            usage()
            exit()
        if ('c' or 'C') in CmdLinePar[1:]:
            isClean = False
        if ('q' or 'Q') in CmdLinePar[1:]:
            isClipMask = False
        if ('s' or 'S') in CmdLinePar[1:]:
            isSky = False
        if ('l' or 'L') in CmdLinePar[1:]:
            Sigmas = float(input("Enter the number of standard deviations for the clipping limit [default value is 5]: "))
        if ('m' or 'M') in CmdLinePar[1:]:
            Median1 = int(input("Enter the size of the median filter window (the default value is 11): "))
        if ('p' or 'P') in CmdLinePar[1:]:
            InterW1 = float(input("Enter the first wavelength in the cropped spectra: "))
            InterW2 = float(input("Enter the last wavelength in the cropped spectra: "))
            isInterp = True
        if ('b' or 'B') in CmdLinePar[1:]:
            Bins = int(input("Enter the binning factor (integer number of points): "))
            isBin = True

    elif CmdLinePar[0] == '+':
        if ('s' or 'S') in CmdLinePar[1:]:
            isSky = True
            SkyStrengthLimit = float(input("Enter the Strength limit of the skylines to be cut off [default value is 0.1]: "))
            SkyStrengthLimit = float(input("Enter the range of wavelengths around skylines to be cut off [default value is 1.25 A]: "))

while (not FileList):
    #s1 = -1
    FileName = input("Enter the spectrum or @list filename: ")
    if FileName[0] == '@':
        s1 = 1
        FileList = True
    if not os.path.isfile(FileName[s1:]):
        print("The File ",FileName," doesn't exist. Exiting...")
        exit(-1)
    try:
        if FileList:
            infile = open(FileName[s1:], "r")
            lines = infile.readlines()
            infile.close()
            print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
        else:
            #lines = FileName[s1:]
            lines.append(FileName[s1:])
            FileList = True
            spec_txt = "spectrum"
    except:
        print("Something wrong with the FileName/FileList ",FileName[s1:])   
#NumOfSpec = len(lines)

WaveFirst, isErr, isMask = FirstSpec(lines[0].rstrip('\n'),FileNameSky)

if isClipMask:
    print('\nThe '+spec_txt+' will be masked.')
else:
    print('\nThe '+spec_txt+' will NOT be masked.')
if isClean:
    print('The '+spec_txt+' will be cleaned with Median',Median1,' and sigma',Sigmas)
else:
    print('The '+spec_txt+' will NOT be cleaned.')
if isSky:
    print('The sky lines will be cut off.')
else:
    print('The sky lines will NOT be cut off.')  
if isInterp:
    print('The '+spec_txt+' will be cropped between',InterW1,'and',InterW2)
else:
    print('The '+spec_txt+' will NOT be cropped.')
if isBin:
    print('The '+spec_txt+' will be rebinned merging ', Bins,' bins together.')
else:
    print('The '+spec_txt+' will NOT be rebinned.')

CleanSpectra(lines,WaveFirst,isErr,isMask)
print('\nDone.')

#!/usr/bin/python
#import numpy

from numpy import *
import sys
from sys import stdin
import numpy as np
from pylab import *
from scipy import signal
from scipy import interpolate
from astropy.stats import sigma_clip
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
def WriteResults(nn,aa,bb,cc,dd,output_file_path):
    """
    Write three columns of data to an external ASCII text file
    """
    outfile = open(output_file_path,"w")
    outfile.write('#    FileName  \t   Phase \t+/- PhaseStd \t Number_of_spectra\n')
    for i in range (0, nn):
        outfile.write(' %s \t %8.4f \t %8.4f \t %d\n' %  (aa[i],bb[i],cc[i],dd[i]))
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

def ReadSpectra(FileNames,FirstWave,isErr,isMask):
    """
    Write three columns of data to an external ASCII text file
    """
    
    AllFlux=[]
    AllFluxErr=[]
    for line in FileNames:
        FileName = line.rstrip('\n')
        data0 = loadtxt(FileName, unpack=True,skiprows=0)
        Wavelength=data0[0,:]
        Length=len(Wavelength)
        Flux=data0[1,:]
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
                    #if FluxNew[i] < 0.0:
                        #FluxNew[i] = FluxMed[i]
                FluxClipNew = interp(FirstWave,WaveClip,FluxClip)
                FluxErrClipNew = interp(FirstWave,WaveClip,FluxErrClip)
            else:
                FluxClipNew = FluxMaskNew
                FluxErrClipNew = FluxErrMaskNew
        else:
            FluxClipNew = Flux
            FluxErrClipNew = np.ones_like(Flux)
            
        #if isSky:
            #WaveClip=[]
            #FluxClip=[]
            #FluxErrClip=[]
            #for i in range(len(FirstWave)):
                #if Sky[i] == 0:
                    #WaveClip.append(Wavelength[i])
                    #FluxClip.append(FluxClipNew[i])
                    #FluxErrClip.append(FluxErrClipNew[i])
                ##else:
                    ##print(Wavelength[i])
            #f1 = interpolate.interp1d(WaveClip,FluxClip)
            #f2 = interpolate.interp1d(WaveClip,FluxErrClip)
            #FluxClipNew = f1(FirstWave)
            #FluxErrClipNew = f2(FirstWave)
            ## FluxMaskNew = interp(FirstWave,WaveClip,FluxClip)
            ## FluxErrMaskNew = interp(FirstWave,WaveClip,FluxErrClip)

        AllFlux.append(FluxClipNew)
        AllFluxErr.append(FluxErrClipNew)
        
    #print(array(AllFlux))

    #if Interp:
        #wa,fl,flerr = SpecCut(Wavelength,FluxClipNew,FluxErrClipNew,InterW1,InterW2)
        #WriteData3(len(wa),wa,fl,flerr,FileName+'.cut')

    return array(AllFlux),array(AllFluxErr)

#########################################################################


def FirstSpec(FileName):
    """
    Reading of the wavelengths from the first spectrum and Testing for the presense of columns.
    """
    global isClipMask, isClean, isWeighted, isInterp, InterW
    data0 = loadtxt(FileName, unpack=True,skiprows=0)
    Wavelength=data0[0,:]
    Length=len(Wavelength)
    Flux=data0[1,:]
    try:
       FluxErr=data0[2,:]
       isErr = 1
       print ("The first spectrum has the 3rd column. It will be assumed that it is the error data.")
    except IndexError:
       isErr = 0
       isWeighted = False
       isClean = False
       #SpecErr = data1[1,:]
    try:
       Mask=data0[3,:]
       isMask = 1
       print ("The first spectrum has the 4th column. It will be assumed that the files include mask data.")   
    except IndexError:
       isMask = 0
       isClipMask = False
       
    #print("\nFirstSpec--->",len(Wavelength),len(Sky))    
    if isInterp:               
        idx_beg=find_nearest(Wavelength,InterW[0])
        idx_end=find_nearest(Wavelength,InterW[1])+1   
        Wavelength = Wavelength[idx_beg:idx_end]
        #Sky = Sky[idx_beg:idx_end]
    #else:
        #idx_beg=0
        #idx_end=-1
        ##NewWave = Wavelength
    #print("FirstSpec--->",len(Wavelength),len(Sky))
    
    return Wavelength,isErr,isMask


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

def renorm(num,values, errors):
    """
    Renormalize the spectra to have the same average values.
    num -- the number of spectra.
    values, errors -- Numpy ndarrays with the same shape.
    """

    normval = np.empty_like(values)
    normerr = np.empty_like(errors)
    
    TotalAver = mean(values)
    for i in range(num):
        normval[i,:] = values[i,:] / mean(values[i,:]) * TotalAver
        normerr[i,:] = errors[i,:] / mean(values[i,:]) * TotalAver
    
    return (normval, normerr)


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

def weighted_avg_std_sigmaclip(values, weights, isSigmaClip, sigma=5):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    Perform sigma-clipping on the provided data.
    """
    if isSigmaClip:
        #print("sigma=",sigma)
        maskedarray = sigma_clip(values, sigma=sigma, maxiters=None, cenfunc=mean, masked=True, copy=False)
        #print(len(maskedarray),maskedarray.size,"\n",maskedarray)        
        #new_vals = []
        new_weights = []
        new_vals = maskedarray.compressed()        
        if len(new_vals)< maskedarray.size:
            #print(values)
            #print(len(maskedarray),maskedarray.size,"\n",maskedarray)        
            #print(len(new_vals),new_vals.size,"\n",new_vals)
            for i in range(maskedarray.size):
                if not maskedarray[:][i]:
                    new_weights.append(weights[i])
            print("Before sigma-clipping: ",len(values)," After: ",len(new_vals))
            exit()
        else:
            new_weights = weights
        average = np.average(new_vals, weights=new_weights)
        variance = 1. / np.sum(new_weights)
    else:
        average = np.average(values, weights=weights)
        variance = 1. / np.sum(weights)
    return (average, math.sqrt(variance))

#########################################################################


def print_header():
    print ("")
    print ("*****************************************************************************************************")
    print ("**                                       spec_phasefolding.py                                      **")
    print ("**                       A utility to perform phase-folding of spectral TXT-files.                 **")
    print ("**              The spectra can be first cleaned from bad data points and outliers and then        **")
    print ("**                        averaged using a weighted or ordinary arithmetic mean.                   **")
    print ("**                                           2021-May-21                                           **")
    print ("**                                        Vitaly Neustroev                                         **")
    print ("*****************************************************************************************************")
    print ("")

def print_help():
    print ("Help: %s -h\n" % sys.argv[0])

def usage():
    print_header()
    "Usage function"
    print ("Usage: %s [options] [@]FileList [ResultFile]" % sys.argv[0])
    print (" ")
    print ("FileList is a list of spectra with filenames in the 1st column and phases in the 2nd column.")
    print ("Each spectrum can consist of 4 columns: Wavelength, Flux, FluxErr QualityMask")
    print ("                          (for example, spectra obtained with ESO/X-shooter).")
    print ("The folding of spectra can be based on a phase interval (default) or by adding lines with")
    print ("                             the ['*' or '#'] symbol in the first column to the FileList.")
    print ("\nOptions: =[Fi|Fs] -cnqw +lmpbdD")
    print ("     =Fs: The folding of spectra based on marks in the FileList.")
    print ("     =Fi: The folding of spectra based on a phase interval [default method]. if initiated,")
    print ("                the user will be asked to enter the phase interval [default value is 0.1].")
    print("")
    print ("     -h: Help (this message)")
    print ("     -c: The spectra will NOT be cleaned    [default: will be cleaned]")
    print ("     -n: The spectra will NOT be normalized [default: will be normalized]")
    print ("     -q: The spectra will NOT be clipped    [default: will be clipped using the QualityMask]")
    #print ("     -s: The spectra will NOT be clipped off the sky lines [default: will be clipped using") 
    #print ("                                 the wavelengths of skylines from the file 'skylines.txt']")
    print ("     -w: An ordinary arithmetic average will be used instead of a weighted average")
    print ("     +l: Will be asked to enter the number of standard deviations at the given wavelength") 
    print ("                                        to use for the clipping limit [default value is 5]")
    print ("     +m: Will be asked to enter the size of the median filter window [default value is 11]")
    print ("     +p: Will be asked to enter the first and the last wavelengths in cropped spectra")
    print ("            (the cropping will be applied after the cleaning and clipping procedures)")
    print ("     +b: Will be asked to enter the binning factor (integer number of points)")
    print ("                         (the rebinning will be applied after all procedures)")
    print ("     +d: When folding the spectra, sigma-clipping will be performed for each wavelength") 
    print ("                     [default number of standard deviations to use for the sigma-clipping limit is 5]")
    print ("     +D: Will be asked to enter the number of standard deviations to use for the sigma-clipping limit") 
    print ("                                                                            when folding the spectra.")
    print ("")
    sys.exit(-1)


##########################################################################


#os.chdir("/home/benj/OwnCloud/My_Papers/BW_Scl/Data/Xsh/NIR/NewReduction/Telluric_Cleaned/Phased/Test/")
#os.chdir("/home/benj/OwnCloud/Scripts/T/")
#os.chdir("/home/benj/Work/Analysis/BW_Scl/PythonTests")

#print_header()

global isClipMask, isClean, isWeighted, isInterp, InterW
FileList=False
ResultFileName=False
isClipMask = True
isClean = True
isRenorm = True
isWeighted = True
isInterp = False
isBin = False
isSigmaClip = False
InterW = [None,None]
HomeDir = "/scisoft/Other_Soft/Files4scripts/"

s1 = 0
Median1 = 11
Sigmas = 5.0
SigmaClip = 5.0
FoldingType = "PhaseIntervals"
Ph_step = 0.10


if len(sys.argv) == 1:
    usage()
else:
    print_header()
    print_help()

for i in range(len(sys.argv)-1):
    CmdLinePar = sys.argv[i+1]
    # print('CmdParameter=',CmdLinePar,'First symbol: ',CmdLinePar[0])
    # print(CmdLinePar[0])
    if ((CmdLinePar[0] != '-') and (CmdLinePar[0] != '+')) and (CmdLinePar[0] != '='):
        if not FileList:
            FileName = CmdLinePar
            if FileName[0] == '@':
                s1 = 1
            try:
                infile = open(FileName[s1:], "r")
                lines = infile.readlines()
                infile.close()
                #print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
                FileList=True                
            except:
                print("Something wrong with the FileList ",FileName[s1:])
        else:
            FoldedFileNames = CmdLinePar.rstrip('\n')
            if not os.path.isfile(FoldedFileNames):                
                ResultFileName = False
            else:
                while not os.path.isfile(FoldedFileNames):
                    print("The File ",FoldedFileNames," already exists.")
                    FoldedFileNames = input("Enter the file name for the result: ")
                ResultFileName = True            

    elif CmdLinePar[0] == '-':
        if ('h' or 'H') in CmdLinePar[1:]:
            usage()
            exit()
        if ('c' or 'C') in CmdLinePar[1:]:
            isClean = False
        if ('q' or 'Q') in CmdLinePar[1:]:
            isClipMask = False
        #if ('s' or 'S') in CmdLinePar[1:]:
            #isSky = False
        if ('n' or 'N') in CmdLinePar[1:]:
            isRenorm = False
        if ('w' or 'W') in CmdLinePar[1:]:
            isWeighted = False
        if ('l' or 'L') in CmdLinePar[1:]:
            Sigmas = float(input("Enter the number of standard deviations for the clipping limit [default value is 5]: "))
        if ('m' or 'M') in CmdLinePar[1:]:
            Median1 = int(input("Enter the size of the median filter window (the default value is 11): "))
        if ('p' or 'P') in CmdLinePar[1:]:
            InterW[0] = float(input("Enter the first wavelength in the cropped spectra: "))
            InterW[1] = float(input("Enter the last wavelength in the cropped spectra: "))
            isInterp = True
        if ('b' or 'B') in CmdLinePar[1:]:
            Bins = int(input("Enter the binning factor (integer number of points): "))
            isBin = True
    elif CmdLinePar[0] == '+':
        if ('h' or 'H') in CmdLinePar[1:]:
            usage()
            exit()
        if ('c' or 'C') in CmdLinePar[1:]:
            isClean = False
        if ('q' or 'Q') in CmdLinePar[1:]:
            isClipMask = False
        #if ('s' or 'S') in CmdLinePar[1:]:
            #isSky = False
        if ('n' or 'N') in CmdLinePar[1:]:
            isRenorm = False
        if ('w' or 'W') in CmdLinePar[1:]:
            isWeighted = False
        if ('l' or 'L') in CmdLinePar[1:]:
            Sigmas = float(input("Enter the number of standard deviations for the clipping limit [default value is 5]: "))
        if ('m' or 'M') in CmdLinePar[1:]:
            Median1 = int(input("Enter the size of the median filter window (the default value is 11): "))
        if ('p' or 'P') in CmdLinePar[1:]:
            InterW[0] = float(input("Enter the first wavelength in the cropped spectra: "))
            InterW[1] = float(input("Enter the last wavelength in the cropped spectra: "))
            isInterp = True
        if ('b' or 'B') in CmdLinePar[1:]:
            Bins = int(input("Enter the binning factor (integer number of points): "))
            isBin = True
        if ('d' or 'D') in CmdLinePar[1:]:
            Bins = float(input("Enter the number of standard deviations to use for the sigma-clipping limit: "))
            isSigmaClip = True
    elif CmdLinePar[0] == '=':
        if CmdLinePar[1:3].lower() == 'fi':
            Ph_step = float(input("Enter the phase interval for phase-folding: "))
            FoldingType = "PhaseIntervals"
        elif CmdLinePar[1:3].lower() == 'fs':
            FoldingType = "SymbolBreaks"
        else:
            print(CmdLinePar[0:3]," --> unknown option. Phase intervals of",Ph_step,"will be used for phase-folding.")
            
print("\n")
while (not FileList):
    FileName = input("Enter the file name of a list of spectra: ")
    while not os.path.isfile(FileName):
        print("The File ",FileName," doesn't exist.")
        FileName = input("Enter the file name of a list of spectra: ")
    try:
        infile = open(FileName[s1:], "r")
        lines = infile.readlines()
        infile.close()
        #print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
        FileList=True                
    except:
        print("Something wrong with the FileList ",FileName[s1:])
#NumOfSpec = len(lines)

if not ResultFileName:
    i = 1
    FoldedFileNames = FileName[s1:] + ".folded"
    while os.path.isfile(FoldedFileNames):
        FoldedFileNames = FileName[s1:] + '.' + str(i) + ".folded"
        i+=1
ResultFileName = True            



specnum=0
while (lines[specnum][0] == '*') or (lines[specnum][0] == '#'):
    specnum += 1
#print(lines[specnum].rstrip('\n'))
temptxt=lines[specnum].split()
#print(temptxt[0])
#print(float(temptxt[1])+1.0)
WaveFirst, isErr, isMask = FirstSpec(temptxt[0])
#WaveFirst, isErr, isMask = FirstSpec(lines[specnum].rstrip('\n'),FileNameSky)

#exit()

if FoldingType == "SymbolBreaks":
    print("\nThe folding of spectra will be based on marks in the FileList.")
elif FoldingType == "PhaseIntervals":
    print("\nPhase intervals of",Ph_step,"will be used for phase-folding.")

if isClipMask:
    print('\nThe spectra will be masked.')
else:
    print('\nThe spectra will NOT be masked.')

#if isSky:
    #print('The sky lines will be cut off.')
#else:
    #print('The sky lines will NOT be cut off.')
    
if isClean:
    print('The spectra will be cleaned with Median',Median1,' and sigma',Sigmas)
else:
    print('The spectra will NOT be cleaned.')

if isRenorm:
    print('The spectra will be normalized before averagin.')
else:
    print('The spectra will NOT be normalized before averagin.')

if isWeighted:
    print('The resulting spectra will be the WEIGHTED average.')
else:
    print('The resulting spectra will be the ordinary arithmetic average.')

if isSigmaClip:
    print('When folding the spectra, sigma-clipping will be performed for each wavelength with sigma=',SigmaClip)
else:
    print('When folding the spectra, sigma-clipping will NOT be performed.')
    
if isInterp:
    print('The resulting spectra will be cropped between',InterW[0],'and',InterW[1])
else:
    print('The resulting spectra will NOT be cropped.')

if isBin:
    print('The resulting spectra will be rebinned merging ', Bins,' bins together.')
else:
    print('The resulting spectra will NOT be rebinned.')

print('\nThe list of resulting spectra with the corresponded phases will be written in the file ',FoldedFileNames)


ResNames  = []
ResPhases = []
ResStd    = []
ResNum    = []

specnum = 0
filenames = []
phases = []
FluxAve = np.empty_like(WaveFirst)
FluxErr = np.empty_like(WaveFirst)
print("\nPhase folding:")
if FoldingType == "SymbolBreaks":
    while specnum <= len(lines):
        if (specnum < len(lines)) and ((lines[specnum][0] != '*') and (lines[specnum][0] != '#')):
            temptxt=lines[specnum].split()
            #print(temptxt[0])
            #print(float(temptxt[1]))
            filenames.append(temptxt[0])
            phases.append(float(temptxt[1]))
        elif len(filenames)>0:
            AllFlux,AllFluxErr = ReadSpectra(filenames,WaveFirst,isErr,isMask)
            ph_mean=mean(phases)
            ph_std = std(phases)
            phstr="{:.4f}".format(ph_mean)
            ResultFile = "ph_"+phstr+".dat"
            for i in range(len(WaveFirst)):
                if min((AllFluxErr[:,i]))>0.0:
                    FluxAve[i],FluxErr[i] = weighted_avg_std_sigmaclip(AllFlux[:,i], 1./(AllFluxErr[:,i])**2, isSigmaClip, SigmaClip)
                    #FluxAve[i],FluxErr[i] = weighted_avg_and_std(AllFlux[:,i], 1./(AllFluxErr[:,i])**2)
            if isBin:
                WaveBin,FluxBin,ErrBin = Rebin(WaveFirst,FluxAve,FluxErr,isErr,Bins)            
                #WaveBin,FluxBin,ErrBin = Rebin(WaveFirst,FluxAve,FluxErr,isErr,Bins)   
            else:
                WaveBin,FluxBin,ErrBin = WaveFirst,FluxAve,FluxErr
                
            if isErr:
                WriteData3(len(WaveBin),WaveBin,FluxBin,ErrBin,ResultFile)
            else:
                WriteData2(len(WaveBin),WaveBin,FluxBin,ResultFile)       
            #if isErr:
                #WriteData3(len(WaveFirst),WaveFirst,FluxAve,FluxErr,ResultFile)
            #else:
                #WriteData2(len(WaveFirst),WaveFirst,FluxAve,ResultFile)       
            print(ResultFile," --> ",len(filenames),"spectra averaged",phases)
            print("phase=","{:.4f}".format(ph_mean),"+/-","{:.4f}".format(ph_std))
            ResNames.append(ResultFile)
            ResPhases.append(ph_mean)
            ResStd.append(ph_std)
            ResNum.append(len(filenames))
            filenames = []
            phases = []
        specnum += 1
elif FoldingType == "PhaseIntervals":
    for line in lines:
        if (line[0] != '*') and (line[0] != '#'):
            temptxt=line.split()
            filenames.append(temptxt[0])
            phases.append(float(temptxt[1]))
    i=0
    while i<len(phases)-1:
        T1=phases[i]
        j=find_nearest(array(phases), T1+Ph_step)
        if phases[j]>T1+Ph_step:
            j=j-1
        AllFlux,AllFluxErr = ReadSpectra(filenames[i:j+1],WaveFirst,isErr,isMask)
        ph_mean= mean(phases[i:j+1])
        ph_std = std(phases[i:j+1])
        phstr="{:.4f}".format(ph_mean)
        ResultFile = "ph_"+phstr+".dat"
        for k in range(len(WaveFirst)):
            if min((AllFluxErr[:,k]))>0.0:
                FluxAve[k],FluxErr[k] = weighted_avg_std_sigmaclip(AllFlux[:,k], 1./(AllFluxErr[:,k])**2, isSigmaClip, SigmaClip)
        if isBin:
            WaveBin,FluxBin,ErrBin = Rebin(WaveFirst,FluxAve,FluxErr,isErr,Bins)            
        else:
            WaveBin,FluxBin,ErrBin = WaveFirst,FluxAve,FluxErr           
        if isErr:
            WriteData3(len(WaveBin),WaveBin,FluxBin,ErrBin,ResultFile)
        else:
            WriteData2(len(WaveBin),WaveBin,FluxBin,ResultFile)       
        #print(ResultFile," --> ",len(phases[i:j+1]),"spectra averaged:",phases[i:j+1])
        print("ph=","{:.4f}".format(ph_mean),"±","{:.4f}".format(ph_std),"-->",len(phases[i:j+1]),"spectra averaged:",phases[i:j+1])
        ResNames.append(ResultFile)
        ResPhases.append(ph_mean)
        ResStd.append(ph_std)        
        ResNum.append(len(phases[i:j+1]))
        i=j+1
    
WriteResults(len(ResNames),ResNames,ResPhases,ResStd,ResNum,FoldedFileNames)
print("\nDone.")        
exit()


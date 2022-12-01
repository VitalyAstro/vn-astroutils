#!/usr/bin/env python3
"""spec-aver: ver. 20221130  (c) Vitaly Neustroev"""
DateVer='1.07:  2022-Dec-01'

from numpy import *
import sys
from sys import stdin
import numpy as np
from pylab import *
from scipy import signal
from scipy import interpolate
import os
import os.path
import glob

# ==============================================================================
def print_history():
    print('2022-Dec-01: Patterns can be used instead of list of files.')
    print('2022-Nov-30: Different standard deviations of the average spectrum can be used.')
    print('2022-Nov-25: Weighted standard deviations of the average spectrum are now calculated correctly.')

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
def WriteData5(nn,aa,bb,cc,dd,ee,output_file_path):
    """
    Write five columns of data to an external ASCII text file
    """
    outfile = open(output_file_path,"w")
    for i in range (0, nn):
        outfile.write(' %12.6f \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i],dd[i],ee[i]))
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

def SpecSkyCut(aa,bb,cc,W1,W2):
    """
    A spectrum cut off the sky lines.
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
    Read the spectra from ASCII text files
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

        if isSky:
            WaveClip=[]
            FluxClip=[]
            FluxErrClip=[]
            for i in range(Length):
                if Sky[i] == 0:
                    WaveClip.append(Wavelength[i])
                    FluxClip.append(FluxClipNew[i])
                    FluxErrClip.append(FluxErrClipNew[i])
                #else:
                    #print(Sky[i],Wavelength[i])
            #print(WaveClip)
            f1 = interpolate.interp1d(WaveClip,FluxClip)
            f2 = interpolate.interp1d(WaveClip,FluxErrClip)
            FluxClipNew = f1(FirstWave)
            FluxErrClipNew = f2(FirstWave)
            # FluxMaskNew = interp(FirstWave,WaveClip,FluxClip)
            # FluxErrMaskNew = interp(FirstWave,WaveClip,FluxErrClip)

        AllFlux.append(FluxClipNew)
        AllFluxErr.append(FluxErrClipNew)

    return array(AllFlux),array(AllFluxErr)

#########################################################################


def FirstSpec(FileName,FileNameSky):
    """
    Reading of the wavelengths from the first spectrum and Testing for the presense of columns.
    """
    global isClipMask, isClean, isWeighted, isSky, Sky, dLam, SkyStrengthLimit
    data0 = loadtxt(FileName, unpack=True,skiprows=0)
    Wavelength=data0[0,:]
    Length=len(Wavelength)
    Flux=data0[1,:]
    try:
       FluxErr=data0[2,:]
       isErr = 1
       print ("The first spectrum has the 3rd column. We assume that ALL the files include error data.")
    except IndexError:
       isErr = 0
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
                #if (WaveSky[i] < Wavelength[0]) or (WaveSky[i] > Wavelength[-1]):
                if (StrengthSky[i] <= SkyStrengthLimit):
                    #print(WaveSky[i],Wavelength[0])
                    continue
                else:
                    idx_beg=find_nearest(Wavelength,WaveSky[i]-dLam)
                    idx_end=find_nearest(Wavelength,WaveSky[i]+dLam)
                    Sky[idx_beg:idx_end+1] = 2048
                    #print(idx_beg,idx_end+1,Wavelength[idx_beg:idx_end+1],Sky[idx_beg:idx_end+1])
            Sky[0]  = 0
            Sky[-1] = 0
            print('\nThe sky lines will be cut off.')
        except:
            print ("\nThe file",FileNameSky,"with sky lines is NOT found. NO sky lines will be cut off.")
            isSky = False
    else:
        print('\nThe sky lines will NOT be cut off.')
    return Wavelength,isErr,isMask


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
        normval[i,:] = values[i,:] / median(values[i,:]) * TotalAver
        normerr[i,:] = errors[i,:] / median(values[i,:]) * TotalAver

    return (normval, normerr)


#########################################################################

#def weighted_avg_and_std_old(values, weights):
    #"""
    #Return the weighted average and standard deviation.
    #values, weights -- Numpy ndarrays with the same shape.
    #"""
    #average = np.average(values, weights=weights)
    #variance = 1. / np.sum(weights)
    #return (average, math.sqrt(variance))

##########################################################################

#def weighted_avg_and_std(values, weights, HowErr):
    #"""
    #Return the weighted average and standard deviation.
    #values, weights -- Numpy ndarrays with the same shape.
    #"""
    #average = np.average(values, weights=weights)
    #NZ=np.count_nonzero(weights)
    #if HowErr == 2:
        #variance = 1. / np.sum(weights)
    #else:
        #variance = np.sum(weights*(values-average)**2) / (NZ-1) / np.sum(weights)
    #return (average, math.sqrt(variance))

##########################################################################

#def weighted_std(values, weights, HowErr):
    #"""
    #Return the weighted standard deviation.
    #values, weights -- Numpy ndarrays with the same shape.
    #"""
    #average = np.average(values, weights=weights)
    #NZ=np.count_nonzero(weights)
    #if HowErr == 2:
        #variance = 1. / np.sum(weights)
    #else:
        #variance = np.sum(weights*(values-average)**2) / (NZ-1) / np.sum(weights)
    #return math.sqrt(variance)




def print_header():
    print ("")
    print ("***********************************************************************************************")
    print ("**                                     spec_aver.py                                          **")
    print ("**                               ver. %s                                     **" % DateVer)
    print ("** A utility to combine TXT-files of spectra using a weighted or ordinary arithmetic average **")
    print ("**            The spectra can be first cleaned from bad data points and outliers             **")
    print ("**                                   Vitaly Neustroev                                        **")
    print ("***********************************************************************************************")
    print ("")

def usage():
    print_header()
    "Usage function"
    print ("Usage: %s [options] @FileList/pattern [ResultFile]" % sys.argv[0])
    print (" ")
    print ("@FileList is a list of spectra.")
    print ("If no FileList is given then a list will be created from the files matching the pattern.")
    print ("   Sometimes a pattern is not parsed correctly. Put it with single quotes, e.g. '*.txt'.")
    print ("\nEach spectrum can consist of 4 columns: Wavelength, Flux, FluxErr QualityMask")
    print ("                           (for example, spectra obtained with ESO/X-shooter)")
    print ("The resulting spectrum will have the 3rd column with FluxErr which by default is calulated using")
    print ("                            the formula SQRT(sum(weights*(values-average)^2)/sum(weights)/(N-1))")
    print ("The option -e will switch it to SQRT(1/sum(weights)) while with -E the errors will be calculated")
    print ("                      from the resulting spectrum (default approach for the data with no errors)")
    print ("Options: -cnqwlmseE +s")
    print ("     -h: Help")
    print ("     -H: Summary of changes")
    print ("     -c: The spectra will NOT be cleaned    [default: will be cleaned]")
    print ("     -n: The spectra will NOT be normalized [default: will be normalized]")
    print ("     -q: The spectra will NOT be clipped    [default: will be clipped using the QualityMask]")
    print ("     -s: The spectra will NOT be clipped off the sky lines [default: will be clipped using")
    print ("                                 the wavelengths of skylines from the file 'skylines.txt']")
    print ("     +s: Will be asked to enter the Strength limit of the skylines and the range of wavelengths")
    print ("                                          around skylines to be cut off [default: 0.1 and 1.25]")
    print ("     -w: An ordinary arithmetic average will be used instead of a weighted average")
    print ("     -l: Will be asked to enter the number of standard deviations at the given wavelength")
    print ("                                        to use for the clipping limit [default value is 5]")
    print ("     -m: Will be asked to enter the size of the median filter window [default value is 21]")
    print ("     -e or -E: the method of error calculation (see above)")
    print ("")
    sys.exit(-1)


##########################################################################


if __name__ == '__main__':

    #try:
       #os.chdir("d:\\VBN\\OwnCloud\\My_Papers\\SSS122222\\VLT\\2019\\Tests")
       #retval = os.getcwd()
       #print ("Windows directory changed successfully: %s" % retval)
    #except:
       #try:
           #os.chdir("/home/benj/OwnCloud/My_Papers/SSS122222/VLT/2019/Tests")
           ### Check current working directory.
           #retval = os.getcwd()
           #print ("Linux directory changed successfully: %s" % retval)
       #except:
           #print ("The directory could not change.")
    
    global isClipMask, isClean, isWeighted, isSky, Sky, dLam, SkyStrengthLimit
    FileList=False
    isTemplate = False
    ResultFileName=False
    isClipMask = True
    isClean = True
    isRenorm = True
    isWeighted = True
    isInterp = False
    isSky = True
    Sky = []
    ScriptDir = "/scisoft/Other_Soft/Files4scripts/"
    if os.path.isfile('./skylines.txt'):
        FileNameSky = "./skylines.txt"
    else:
        FileNameSky = ScriptDir+"skylines.txt"
    
    dLam = 1.25
    SkyStrengthLimit = 0.1
    s1 = 0
    Median1 = 21
    Sigmas = 5.0
    HowError = 1
    
    if len(sys.argv) == 1:
        usage()   
    
    for i in range(len(sys.argv)-1):
        CmdLinePar = sys.argv[i+1]
        if (CmdLinePar[0] != '-') and (CmdLinePar[0] != '+'):
            if not FileList:
                FileName = CmdLinePar
                if FileName[0] == '@':
                    s1 = 1
                    try:
                        infile = open(FileName[s1:], "r")
                        lines = infile.readlines()
                        infile.close()
                        print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
                        FileList=True
                    except:
                        print("Something wrong with the FileList ",FileName[s1:])
                else:
                    isTemplate = True
                    lines=glob.glob(FileName)
                    if len(lines)>0:
                        print("\n",len(lines),"file(s) are found.")
                        FileList=True
                    else:
                        print("No file(s) matching the template are found.")
            else:
                ResultFile = CmdLinePar.rstrip('\n')
                if not os.path.isfile(ResultFile):
                    ResultFileName = True
                else:
                    while not os.path.isfile(ResultFile):
                        print("The File ",ResultFile," already exists.")
                        ResultFile = input("Enter the file name for the result: ")
                    ResultFileName = True
        elif CmdLinePar[0] == '-':
            if ('h') in CmdLinePar[1:]:
                usage()
                exit()
            if ('H') in CmdLinePar[1:]:
                print_history()
                exit()
            if ('c' or 'C') in CmdLinePar[1:]:
                isClean = False
            if ('q' or 'Q') in CmdLinePar[1:]:
                isClipMask = False
            if ('s' or 'S') in CmdLinePar[1:]:
                isSky = False
            if ('n' or 'N') in CmdLinePar[1:]:
                isRenorm = False
            if ('w' or 'W') in CmdLinePar[1:]:
                isWeighted = False
            if ('l' or 'L') in CmdLinePar[1:]:
                Sigmas = float(input("Enter the number of standard deviations for the clipping limit [default value is 5]: "))
            if ('m' or 'M') in CmdLinePar[1:]:
                Median1 = int(input("Enter the size of the median filter window (the default value is 21): "))
            if ('e') in CmdLinePar[1:]:
                HowError = 2
            if ('E') in CmdLinePar[1:]:
                HowError = 3
    
        elif CmdLinePar[0] == '+':
            if ('s' or 'S') in CmdLinePar[1:]:
                isSky = True
                SkyStrengthLimit = float(input("Enter the Strength limit of the skylines to be cut off [default value is 0.1]: "))
                SkyStrengthLimit = float(input("Enter the range of wavelengths around skylines to be cut off [default value is 1.25 A]: "))
    
    while (not FileList):
        FileName = input("Enter the file name of a list of spectra: ")
        while not os.path.isfile(FileName):
            print("The File ",FileName," doesn't exist.")
            FileName = input("Enter the file name of a list of spectra: ")
        try:
            infile = open(FileName[s1:], "r")
            lines = infile.readlines()
            infile.close()
            print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
            FileList=True
        except:
            print("Something wrong with the FileList ",FileName[s1:])
    NumOfSpec = len(lines)
    
    if not ResultFileName:
        i = 1
        if isTemplate:
            ResultFile = "AveragedSpectrum"
        else:
            ResultFile = FileName[s1:]
        ResultFileTemp = ResultFile
        while os.path.isfile(ResultFileTemp + ".dat"):
            ResultFileTemp = ResultFile + '.' + str(i)
    ##        ResultFile = FileName[s1:] + '.' + str(i) + ".dat"
            i+=1
        ResultFile = ResultFileTemp + ".dat"
    ResultFileName = True
    
    
    
    WaveFirst, isErr, isMask = FirstSpec(lines[0].rstrip('\n'),FileNameSky)
    
    if not isErr:
        HowError = 3
    
    if isClipMask:
        print('The spectra will be masked.')
    else:
        print('The spectra will NOT be masked.')
    
    #isSkyFileExists = os.path.isfile(FileNameSky)
    
    #if isSky and isSkyFileExists:
        #print('The sky lines will be cut off.')
    #elif not isSky:
        #print('The sky lines will NOT be cut off.')
    #else:
        #print ("The file with sky lines is not found. No sky lines will be cut off.")
        #isSky = False
    
    if isClean:
        print('The spectra will be cleaned with Median',Median1,' and sigma',Sigmas)
    else:
        print('The spectra will NOT be cleaned.')
    
    if isRenorm:
        print('The spectra will be normalized.')
    else:
        print('The spectra will NOT be normalized.')
    
    if isWeighted:
        print('The resulting spectrum will be the WEIGHTED average.')
    else:
        print('The resulting spectrum will be the ordinary arithmetic average.')
    
    print('\nThe resulting spectrum will be written in the file ',ResultFile)
    
    print('\nThe flux errors in the resulting spectrum will be calulated using the')
    if HowError == 3:
        print('spectrum itself.')
    elif HowError == 2:
        print ("formula SQRT(1/sum(weights)).")
    else:
        print ("formula SQRT(sum(weights*(values-average)^2)/sum(weights)/(N-1)).")
    
    if isInterp:
        InterW1 = int(raw_input("Enter the first wavelength in the cut spectrum: "))
        InterW2 = int(raw_input("Enter the last wavelength in the cut spectrum: "))
    
    
    AllFlux,AllFluxErr = ReadSpectra(lines,WaveFirst,isErr,isMask)
    FluxMed=median(AllFlux)
    if isRenorm:
        AllFlux,AllFluxErr = renorm(NumOfSpec,AllFlux,AllFluxErr)
    FluxAve = np.empty_like(WaveFirst)
    FluxErr = np.empty_like(WaveFirst)
    
    for i in range(len(WaveFirst)):
        if min((AllFluxErr[:,i])) == 0.0:
            AllFluxErr[:,i] = AllFluxErr[:,i] + 0.001*FluxMed
        Weights=1./(AllFluxErr[:,i])**2
        FluxAve[i] = np.average(AllFlux[:,i], weights=Weights)
        NZ=np.count_nonzero(Weights)
        if HowError == 2:
            FluxErr[i] = math.sqrt(1. / np.sum(Weights))
        elif HowError == 1:
            #FluxErr[i] = weighted_std(AllFlux[:,i], 1./(AllFluxErr[:,i])**2, HowError)
            FluxErr[i] = math.sqrt(np.sum(Weights*(AllFlux[:,i]-FluxAve[i])**2) / (NZ-1) / np.sum(Weights))
    if HowError == 3:
        FluxErr = MedStd(FluxAve - MedFilt(FluxAve,Median1), Median1)
    
    WriteData3(len(WaveFirst),WaveFirst,FluxAve,FluxErr,ResultFile)


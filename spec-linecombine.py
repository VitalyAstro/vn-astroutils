#!/usr/bin/python
#import numpy
from numpy import *
import sys
from sys import stdin
#import numpy as np
import math
import statistics
from pylab import *
import os


##########################################################################
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
def WriteData1(nn,aa,output_file_path):
    """
    Write one column of data to an external ASCII text file
    """
    output_file = output_file_path.rstrip('\n')
    outfile = open(output_file,"w")
    for i in range (0, nn):
        outfile.write(' %12.6e \n' %  (aa[i]))
#        outfile.write(' %12.6f \t %12.6f \n' %  (aa[i],bb[i]))
    outfile.close()
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

def Tmp(nn,aa,output_file_path):
    """
    Write one column of data to an external ASCII text file
    """   
    try:
        data0 = loadtxt(FileName, unpack=True,skiprows=0)
        Wavelength=data0[0,:]
        Length=len(Wavelength)
        Flux=data0[1,:]
        print("\n",FileName, end = '')
    except:
        print("Something wrong with the file ",FileName)
        print("Please check whether it is a spectrum or not. Skipping...")
        #continue        
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
        if isErr:
            normerr[i,:] = errors[i,:] / mean(values[i,:]) * TotalAver
    
    return (normval, normerr)

#########################################################################


def print_header():
    print ("")
    print ("****************************************************************************************************")
    print ("**                                   spec_linecombine.py                                          **")
    print ("** A utility to combine spectral lines in spectra using a weighted or ordinary arithmetic average **")
    print ("**                                       2021-June-16                                             **")
    print ("**                                     Vitaly Neustroev                                           **")
    print ("****************************************************************************************************")
    print ("")

########################################################################

def usage():
    print_header()
    "Usage function"
    print ("Usage: %s [options] [@]FileList" % sys.argv[0])
    print (" ")
    print ("FileList is a list of spectra.")
    print ("Each spectrum can consist of 3 columns: Wavelength, Flux, and FluxErr.")
    print ("Options: -hnw")
    print ("     -h: Help")
    print ("     -n: The spectra will NOT be normalized [default: will be normalized]")
    print ("     -w: An ordinary arithmetic average will be used instead of a weighted average")
    print ("")
    sys.exit(-1)


##########################################################################

print_header()
#print " "
#print "**********************************************************************"
#print " "
#print "vel-files-rebin.py FileList"
#print " "
#print "The script reads the files 'spectra.dat' & 'wavelength.dat', prepared"
#print "for the 'Velocity.exe' program, and averages the specified spectra."
#print "The result will be written in the 's.dat' file."
#print " "
#print "**********************************************************************"
#print " "

#os.chdir("/home/benj/OwnCloud/Scripts/")
#os.chdir("/home/benj/OwnCloud/Work/IDL/DopMap/BW_Scl/VIS/PhFold20/NewNames/")
FileList=False
isRenorm=True
isWeighted=True
s1 = 0
lines = []

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
        if ('n' or 'N') in CmdLinePar[1:]:
            isRenorm = False
        if ('w' or 'W') in CmdLinePar[1:]:
            isWeighted = False
        #if ('l' or 'L') in CmdLinePar[1:]:
            #Sigmas = float(input("Enter the number of standard deviations for the clipping limit [default value is 5]: "))


c = 299792.458
SpecLines = [5269.541, 5328.042, 5371.493]   # FeI
#SpecLines = [8498.02, 8542.09]     # CaII
#SpecLines = [8498.02, 8542.09, 8662.14]     # CaII
SpecLines.sort()
MedW = statistics.median_high(SpecLines)
print("The lines to be combined: ",SpecLines)
delV = float(input("Enter the velocity range around the center of each line to be combined (in km/s): "))
#w0 = float(input("Enter the central wavelength of the combined line: "))
print("The central wavelength of the combined line will be ",MedW)
if isRenorm:
    print('\nThe lines will be normalized.')
else:
    print('\nThe lines will NOT be normalized.')

if isWeighted:
    print('The combined line will be the WEIGHTED average if the error data is present.')
else:
    print('The combined line will be the ordinary arithmetic average.')

print('The resulting spectra will be written in the files [Filename].combined\n')

for line in lines:
    Fluxes  = []
    FluxErrors = []
    CurFile = line.rstrip('\n')
    ResultFile = CurFile+".combined"
    try:
        data0 = loadtxt(CurFile, unpack=True,skiprows=0)
        Wavelength=data0[0,:]
        Length=len(Wavelength)
        Flux=data0[1,:]
    except:
        print("Something wrong with the file ",line)
        print("Please check whether it is a spectrum or not. Skipping...")
        #exit()
        continue        
    
    idx= find_nearest(Wavelength,SpecLines[-1])
    dV = (Wavelength[idx]-Wavelength[idx-1])/SpecLines[-1]*c
    idx_low  = idx-math.ceil(delV/dV)
    idx_high = idx+math.ceil(delV/dV)
    Velocities=(Wavelength[idx_low:idx_high]-SpecLines[-1])/SpecLines[-1]*c
    Fluxes.append(Flux[idx_low:idx_high])
    FluxAve = np.zeros_like(Velocities)
    ErrAve  = np.zeros_like(Velocities)
    
    try:
        FluxErr = data0[2,:]
        isErr = 1
        print ("The spectrum",CurFile,"has the 3rd column. It will be assumed that it is the error data.")
        FluxErrors.append(FluxErr[idx_low:idx_high])
    except IndexError:
        isErr = 0
        FluxErr = np.ones_like(Velocities)
    
    for SpecLine in SpecLines[0:-1]:
        FluxSpec = interp(Velocities,(Wavelength-SpecLine)/SpecLine*c,Flux)
        Fluxes.append(FluxSpec)
        if isErr:
            ErrSpec = interp(Velocities,(Wavelength-SpecLine)/SpecLine*c,FluxErr)
            FluxErrors.append(ErrSpec)

    Fluxes = np.asarray(Fluxes)
    FluxErrors = np.asarray(FluxErrors)
    
    if isRenorm:
        Fluxes,FluxErrors = renorm(len(SpecLines),Fluxes,FluxErrors)
    
    for i in range(len(Velocities)):
        if (isErr and min((FluxErrors[:,i]))>0.0) and isWeighted:
            FluxAve[i],ErrAve[i] = weighted_avg_and_std(Fluxes[:,i], 1./(FluxErrors[:,i])**2)
        else:
            FluxAve[i] = mean(Fluxes[:,i])
    
    Velocities = Velocities*MedW/c+MedW   
    if isErr:
        WriteData3(len(Velocities),Velocities,FluxAve,ErrAve,ResultFile)
    else:
        WriteData2(len(Velocities),Velocities,FluxAve,ResultFile)

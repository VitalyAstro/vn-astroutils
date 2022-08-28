#! /usr/bin/env python3
# ver. 20210720
#from __future__ import print_function, division
from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt
import os
import sys
from sys import stdin 
from numpy import * 
#from PyAstronomy.pyaC import pyaErrors as PE
#from PyAstronomy.pyasl import _ic



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
        #outfile.write(' %s \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        outfile.write(' %12.6f \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close() 


#######################################################################################
def printHeader():
    """ Print the Header """
    print(" ")
    print("********************************************************************************")
    print("              spec-barycorr: ver. 20210720  (c) Vitaly Neustroev ")

def printUsage():
    """ Print the usage info """
    print(" ")
    print("Usage:")
    print("Usage: %s FileName/@FileList BaryCorr/W1 [W2 Step BaryCorr]" % sys.argv[0])
    print(" ")
    print("FileName is a spectrum filename; @FileList is a list of spectra.")
    print("Each spectrum can consist of 3 columns: Wavelength, Flux, and FluxErr")    
    print(" ")
    print("Parameters:")
    print("  BaryCorr/W1: if this parameter < 50, it will be assumed to be the barycentric")
    print("               RV correction (in km/s).  Corrected wavelengths will be rebinned")
    print("                                                   back to their original format")
    print("               If >50, then it will be the 1st wavelength in the cropped spectra")
    print("  W2: the last wavelength in the cropped spectra")
    print("  dW: the wavelength step in the cropped spectra")
    print("  BaryCorr: if given then corrections for barycentric motion will be applied.")
    print(" ")
    print("********************************************************************************")
    print(" ")

##########################################################################


printHeader()
cl=2.997e5
pi=3.1415926
isHelCorr  = False
isInterp   = False
isFileList = False

s1 = 0
spec_txt = "spectra"


if (len(sys.argv) < 3):    
    printUsage()
    exit(0)
#else:
    #printHeader()


FileName = sys.argv[1]
if FileName[0] == '@':
    s1 = 1
    isFileList = True
if not os.path.isfile(FileName[s1:]):
    print("The File ",FileName," doesn't exist. Exiting...")
    exit(-1)               
try:
    if isFileList:
        infile = open(FileName[s1:], "r")
        lines = infile.readlines()
        infile.close()
        print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
    else:
        lines = []
        lines.append(FileName[s1:])
        FileList = True
        spec_txt = "spectrum"
except:
    print("Something wrong with the FileName/FileList ",FileName[s1:])               


try:
    Temp=float(sys.argv[2])
    if Temp < 50:
        HelCorr = Temp
        isHelCorr = True
        print("Adopted parameters:")
        print("          BaryCorr:",HelCorr)
        
    else:
        InterW1 = Temp
        isInterp = True
except:
    InterW1 = float(input("Enter the first wavelength in the cropped spectra: "))
    isInterp = True

if isInterp:
    try:
        InterW2=float(sys.argv[3])
    except:
        InterW2 = float(input("Enter the last wavelength in the cropped spectra: "))
    try:
        dW=float(sys.argv[4])
    except:
        dW = float(input("Enter the wavelength step in the cropped spectra: "))

    #InterPoints = int(ceil((InterW2 - InterW1) / dW))+1
    #Temp = InterW1+dW*(InterPoints-1)
    #if Temp>InterW2:
        #InterPoints -=1
        #InterW2 = InterW1+dW*InterPoints
    #WaveInterp = linspace(InterW1, InterW2, InterPoints)
    WaveInterp = arange(InterW1, InterW2+dW/2. , dW)
    InterPoints = len(WaveInterp)

    try:
        HelCorr=float(sys.argv[5])
        isHelCorr = True
    except:
        print("No heliocentric correction will be applied.")
        HelCorr = 0

    print("Adopted parameters:")
    print("          BaryCorr:",HelCorr)
    print("                W1:",InterW1)
    if Temp != InterW2:
        print("                W2:",InterW2," (corrected)")
    else:
        print("                W2:",InterW2)
    print("                dW:",dW)



if (HelCorr == 0) and (not isInterp):
    print("Wavelengths don't change! Alter some parameters. Stop executing!")
    exit(-1)

#if isInterp:
    #WaveInterp = linspace(InterW1, InterW2, InterPoints)

for line in lines:
    FileName = line.rstrip('\n')
    isErr = False
    try:
        data0 = loadtxt(FileName, unpack=True,skiprows=0)
        Wavelength=data0[0,:]
        WaveHelCorr = cl*Wavelength/(-HelCorr+cl)
        Length=len(Wavelength)
        Flux=data0[1,:]
        print("\n",FileName, end = '')
    except:
        print("\nSomething wrong with the file ",FileName)
        print("Please check whether it is a spectrum or not. Skipping...")
        continue   
    try:
       FluxErr=data0[2,:]
       isErr = True
    except IndexError:
       isErr = False

    if isInterp:
        if (WaveHelCorr[0]>WaveInterp[0]) or (WaveHelCorr[-1]<WaveInterp[-1]):
            print(" ---> Warning! The provided wavelength range is larger than original.", end = '')
            #print(" ---> Warning! The provided wavelength range is larger than original. The spectrum will be extrapolated.")        
        #print(" ---> Warning! The provided wavelength range is larger than original.")
        #print(" ---> Warning! The provided wavelength range is larger than original. The spectrum will be extrapolated.")
        FluxInterp = interp(WaveInterp, WaveHelCorr,Flux)
        if isErr:
            FluxErrInterp = interp(WaveInterp, WaveHelCorr,FluxErr)
            WriteData3(InterPoints,WaveInterp,FluxInterp,FluxErrInterp,FileName+'.interp')
        else:
            WriteData2(InterPoints,WaveInterp,FluxInterp,FileName+'.interp')
        print(" --> .interp", end = '')
    else:
        FluxInterp = interp(Wavelength, WaveHelCorr,Flux)        
        if isErr:
            FluxErrInterp = interp(Wavelength, WaveHelCorr,FluxErr)
            WriteData3(Length,Wavelength,FluxInterp,FluxErrInterp,FileName+'.barycorr')
        else:
            WriteData2(Length,Wavelength,FluxInterp,FileName+'.barycorr')
        print(" --> .barycorr", end = '')
print("\nDone.")
    

        

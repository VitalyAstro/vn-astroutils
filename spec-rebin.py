#! /usr/bin/env python3
# ver. 20241115
'''
History:
20241115: Added the binning by groups.
'''
#from __future__ import print_function, division
#from PyAstronomy import pyasl
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

########################################################################

def Grouping(xx,yy,yerr,bins,x1st,xLast):
    """
    Reduce XYYerr data using subgroup stats on X
    """
    if x1st>0:
        x1=find_nearest(xx, x1st-(xx[bins]-xx[0])/2)+1
    else:
        x1=0
    if xLast>0:
        x2=find_nearest(xx, xLast+(xx[bins]-xx[0])/2)+1
    else:
        x2=len(xx)

    i=x1
    xnew=[]
    ynew=[]
    yerrnew=[]
    while i+bins <= x2:
        xnew.append(np.mean(xx[i:i+bins]))
        ynew.append(np.mean(yy[i:i+bins]))
        yerrnew.append(np.std(yy[i:i+bins])/bins**0.5)
        i+=bins
#        j+=1
    if i < x2:
        xnew.append(np.mean(xx[i:x2]))
        ynew.append(np.mean(yy[i:x2]))
        yerrnew.append(yerrnew[-1])
    return xnew,ynew,yerrnew



#######################################################################################
def printHeader():
    """ Print the Header """
    print(" ")
    print("**************************************************************************************")
    print("                spec-rebin: ver. 20241115 (c) Vitaly Neustroev   \n")

def printUsage():
    """ Print the usage info """
    print(" ")
    print("Usage:")
    print("Usage: %s [-bNN]/File-of-wavelengths] FileName/@FileList [W1 W2 dW]" % sys.argv[0])
    print(" ")
    print("FileName is a spectrum filename; @FileList is a list of spectra.")
    print("Each spectrum can consist of 3 columns: Wavelength, Flux, and FluxErr")
    print(" ")
    print("Parameters:")
    print("  File-of-wavelengths: if given then other files will be rebinned to its wavelengths")
    print("                                                    (only first column will be read)")
    print("  -bNN: Binning in groups of NN lines. NN must be integer, for example -b5 or -b12.")
    print("  W1: the 1st wavelength  in the cropped spectra")
    print("  W2: the last wavelength in the cropped spectra")
    print("  dW: the wavelength step in the cropped spectra (not used when -b option is applied.)")
    print(" ")
    print("**************************************************************************************")
    print(" ")

##########################################################################


printHeader()
cl=2.997e5
pi=3.1415926
isBin      = False
isInterp   = False
isFileList = False
isFileWave = True

s1 = 0
spec_txt = "spectra"


if (len(sys.argv) < 3):
    printUsage()
    exit(0)
#else:
    #printHeader()


TestFileName = sys.argv[2]
#print('cmdline 2:',TestFileName)
if (TestFileName[0] == '@') and os.path.isfile(TestFileName[1:]):
    FileName = TestFileName[1:]
    isFileList = True
    FileWave = sys.argv[1]
elif os.path.isfile(TestFileName):
    FileName = TestFileName
    FileWave = sys.argv[1]
else:
    #print(TestFileName,'is not a file')
    isFileWave = False
    TestFileName = sys.argv[1]

if isFileWave:
    FileWave = sys.argv[1]
    if (FileWave[0:2] == '-b'):
        try:
            BinsNumber = int(FileWave[2:])
#            print('\nThe recognised Number of Wavelength Bins is',BinsNumber)
            isBin = True
            isFileWave = False
        except:
            print('\nSomething wrong with the -b parameter. No binning will be used.')
            isFileWave = False
    elif not os.path.isfile(sys.argv[1]):
        print("The File-of-Wavelengths ",sys.argv[1]," doesn't exist. Exiting...")
        exit(-1)
    else:
        try:
            data0 = loadtxt(FileWave, unpack=True,skiprows=0)
            WaveInterp=data0[0,:]
            WavePoints=len(WaveInterp)
            #print("\n",FileName, end = '')
        except:
            print("\nSomething wrong with the File-of-Wavelengths ",FileName)
            exit(-1)
            #print("Please check whether it is a spectrum or not. Skipping...")
            #continue
elif (TestFileName[0] == '@') and os.path.isfile(TestFileName[1:]):
    FileName = TestFileName[1:]
    isFileList = True
elif os.path.isfile(TestFileName):
    FileName = TestFileName
else:
    print("Something wrong with the FileName/FileList ")
    print("Exiting...")
    exit(-1)

try:
    if isFileList:
        infile = open(FileName, "r")
        lines = infile.readlines()
        infile.close()
        print("The FileList ",FileName," includes ",len(lines)," FileNames")
    else:
        lines = []
        lines.append(FileName)
        #print("The FileList ",FileName," includes ",len(lines)," FileNames")
        #FileList = True
        #spec_txt = "spectrum"
except:
    print("Something wrong with the FileName/FileList ",FileName)
    exit(-1)



if isFileWave:
    print("\nFile-of-Wavelengths:",FileWave)
    print(WaveInterp)
elif isBin:
    print('\nThe recognised Number of Wavelength Bins is',BinsNumber)
    try:
        InterW1 = float(sys.argv[3])
    except:
        InterW1 = 0
        InterW2 = 0
    if InterW1 > 0:
        try:
            InterW2=float(sys.argv[4])
        except:
            InterW2 = 0
else:
    print("\nThere is no File-of-Wavelengths.")
    try:
        InterW1 = float(sys.argv[2])
    except:
        InterW1 = float(input("Enter the first Wavelength in the cropped spectra: "))
    try:
        InterW2=float(sys.argv[3])
    except:
        InterW2 = float(input("Enter the last Wavelength in the cropped spectra: "))
    try:
        dW=float(sys.argv[4])
    except:
        dW = float(input("Enter the Wavelength step in the cropped spectra: "))


    print("Adopted parameters:")
    print("                W1:",InterW1)
    print("                W2:",InterW2)
    print("                dW:",dW)

    WaveInterp = arange(InterW1, InterW2+dW/2. , dW)
    WavePoints = len(WaveInterp)
    #print(WaveInterp)


#if isFileList:
    #print("\nFileList:",FileName)
#else:
    #print("\nThere is no FileList.")




for line in lines:
    FileName = line.rstrip('\n')
    isErr = False
    try:
        data0 = loadtxt(FileName, unpack=True,skiprows=0)
        Wavelength=data0[0,:]
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
       FluxErr=data0[1,:]

    if isBin:
        WaveInterp,FluxInterp,FluxErrInterp = Grouping(Wavelength,Flux,FluxErr,BinsNumber,InterW1,InterW2)
        WriteData3(len(WaveInterp),WaveInterp,FluxInterp,FluxErrInterp,FileName+'.bin')
        print(" --> .bin", end = '')
    else:
        if (Wavelength[0]>WaveInterp[0]) or (Wavelength[-1]<WaveInterp[-1]):
            print(" ---> Warning! The provided wavelength range is larger than original.", end = '')
        FluxInterp = interp(WaveInterp, Wavelength,Flux)
        if isErr:
            FluxErrInterp = interp(WaveInterp, Wavelength,FluxErr)
            WriteData3(WavePoints,WaveInterp,FluxInterp,FluxErrInterp,FileName+'.interp')
        else:
            WriteData2(WavePoints,WaveInterp,FluxInterp,FileName+'.interp')
        print(" --> .interp", end = '')
print("\nDone.")



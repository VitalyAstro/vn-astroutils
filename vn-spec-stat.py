#!/usr/bin/python

DateVer='2.0:  2025-March-06'
ScriptDir2 = '/scisoft/Other_Soft/Files4scripts/'

#import numpy
from numpy import *
import sys
from sys import stdin
#import numpy as np
from pylab import *
from statistics import mean
import os
import glob
from pathlib import Path
import math
from math import log10


# ==============================================================================
def print_history():
    print('\nHistory:')
    print('2025-Mar-06: User interface significantly improved.')

##########################################################################


def SpecStat(FileName,W1,W2):
    """
    Write three columns of data to an external ASCII text file
    """
    data0 = loadtxt(FileName, usecols=[0,1], unpack=True,skiprows=0)
    #data0 = loadtxt(FileName, usecols=[0,1,2], unpack=True,skiprows=535)
    Length=len(data0[0,:])
    EndSkip=0
    Wavelength=data0[0,:Length-EndSkip]
    Flux=data0[1,:Length-EndSkip]
    FluxNew=data0[1,:Length-EndSkip]
    #FluxErr=data0[2,:Length-EndSkip]
    Length=len(Wavelength)
    idxW1=find_nearest(Wavelength,W1)
    idxW2=find_nearest(Wavelength,W2)
    return mean(Flux[idxW1:idxW2])

#########################################################################


def TestFileList(FileName):
    """
    Testing if a file has more than 1 column
    """
    c=0
    try:
        data0 = loadtxt(FileName)
        Wavelength=data0[0,:]
        c+=1
        try:
            Flux=data0[1,:]
            c+=1
        except Exception:
            pass
    except Exception:
        pass
    return c

#########################################################################


def SpecPhotometry(SpecName,filter_name):
    """
    Calculation of fluxes in different photometric filters from a spectrum
    """
    data0 = loadtxt(SpecName, usecols=[0,1], unpack=True,skiprows=0)
    Wavelength=data0[0,:]
    Length=len(Wavelength)
    Flux=data0[1,:]

    FilePath = CheckFile(filter_name,ScriptDir2)
    WaveFilter,FilterCurve,NumberOfFilters = ReadFilterCurves(FilePath)

    if Wavelength[0] > WaveFilter[0]:
        #print("Attention! The object spectrum is only partially covered by some of filter curves.")
        #print("Some (bluish) fluxes can be underestimated.")
        Wbegin = Wavelength[0]
    else:
        Wbegin = WaveFilter[0]

    if Wavelength[-1] < WaveFilter[-1]:
        #print("Attention! The object spectrum is only partially covered by some of filter curves.")
        #print("Some (reddish) fluxes can be underestimated.")
        Wend = Wavelength[-1]
    else:
        Wend = WaveFilter[-1]

    Wfirst=math.ceil(Wbegin)
    Wlast =math.floor(Wend)

    WaveNew = linspace(Wfirst, Wlast, Wlast-Wfirst+1)
    SpecNew = interp(WaveNew, Wavelength,Flux)

    SpecFlux=[]
    for i in range(NumberOfFilters):
        #print(WaveFilter[0],FilterCurve[i][0],WaveFilter[-1],FilterCurve[i][-1])
        #print(WaveNew[0],WaveNew[-1])
        #FilterTemp = FilterCurve[i]
        #WaveTemp = WaveFilter
        #print(len(WaveTemp),len(FilterCurve[i][:]))
        FilterNew = interp(WaveNew, WaveFilter[:], FilterCurve[i][:])
        #FilterNew = interp(WaveNew, WaveTemp, FilterTemp)
        SpecFilter = SpecNew * FilterNew
        FilterArea = trapz(FilterNew, dx=1.0)
        SpecFilterArea = trapz(SpecFilter, dx=1.0)
        SpecFlux.append(SpecFilterArea / FilterArea)

    return NumberOfFilters,SpecFlux

#########################################################################

def find_nearest(array, value):
    '''
    Find the nearest value inside an array.
    :param array: array
    :param value: desired value (float)
    :return: nearest value and its index
    '''

    idx = (abs(array - value)).argmin()
    return idx

##########################################################################


def print_menu():       ## Your menu design here
#    print("\n")
#    print(43 * "-" , "Filter sets" , 43 * "-")
    print("1. UVOT filters")
    print("2. Standard Johnson-Cousins UBVRI filters (Bessell)")
    print("3. PANSTARRS-PS1 filters")
    print("4. Sloan/SDSS filters")
    print("5. NIR J-H-Ks-K filters")
    print("6. Johnson-Cousins UBVRI filters (Astrodon)")
    print("7. GAIA filters")
    print(99 * "-")

def print_filter_sets():       ## Your menu design here
    print("\n")
    print(15 * "-" , "Filter sets and the corresponding files of the transmission curves:" , 15 * "-")
    print("1. UVOT filters - 'uvot.dat'")
    print("2. Standard Johnson-Cousins UBVRI filters (Bessell) - 'bessell.dat'")
    print("3. PANSTARRS-PS1 filters - 'panstarrs1.dat'")
    print("4. Sloan/SDSS filters - 'sdss.dat'")
    print("5. NIR J-H-Ks-K filters - 'jhksk.dat'")
    print("6. Johnson-Cousins UBVRI filters (Astrodon) - 'ubvri.dat'")
    print("7. GAIA filters - 'gaia.dat'")
    print("\nThe transmission curves must be located in the current or one of the following folders:")
    print(str(Path.home())+'/.vn-astroutils/')
    print(ScriptDir2)
    print(99 * "-")

##########################################################################

def filter_selection():
    loop=True
    while loop:          ## While loop which will keep going until loop = False
        print()
        print(30 * "-" , "MENU: select a filter set" , 30 * "-")
        print("0. Exit")
        print_menu()    ## Displays menu
        choice = input("Enter your choice [1-7]: ")
        try:
            choice = int(choice)
        except:
            continue

        print()
        if choice==0:
            print("Exit")
            return None,None,None
            ## You can add your code or functions here
        elif choice==1:
            print("UVOT filters have been selected. Vega system.")
            FileName = 'uvot.dat'
            FilterNames = ['w2','m2','w1','uu','bb','vv']
            FilterZPs = [5.36274e-9,4.67904e-9,4.14563e-9,3.63751e-9,6.47559e-9,3.72393e-9]
            #VegaSystem = False
            loop=False
            ## You can add your code or functions here
        elif choice==2:
            print("Johnson-Cousins UBVRI filters (Bessell) have been selected. Vega system.")
            FileName = 'bessell.dat'
            FilterNames = ['U','B','V','R','I']
            FilterZPs = [417.5e-11,632.0e-11,363.1e-11,217.7e-11,112.6e-11]
            loop=False
            ## You can add your code or functions here
        elif choice==3:
            print("PANSTARRS-PS1 filters have been selected. AB system.")
            FileName = 'panstarrs1.dat'
            FilterNames = ['g','r','i','z','y']
            FilterZPs = [4810.,6170.,7520.,8660.,9620.]
            VegaSystem = False
            loop=False
            ## You can add your code or functions here
        elif choice==4:
            print("Sloan/SDSS filters have been selected. AB system.")
            FileName = 'sdss.dat'
            FilterNames = ['u','g','r','i','z']
            FilterZPs = [3551,4686,6165,7481,8931]
            VegaSystem = False
            loop=False
            ## You can add your code or functions here
        elif choice==5:
            print("NIR J-H-Ks-K filters have been selected. Vega system.")
            FileName = 'jhksk.dat'
            FilterNames = ['J','H','Ks','K']
            FilterZPs = [31.47e-11,11.38e-11,3.961e-11,3.961e-11]
            loop=False
            ## You can add your code or functions here
            loop=False # This will make the while loop to end as not value of loop is set to False
        elif choice==6:
            print("Johnson-Cousins UBVRI filters (Astrodon) have been selected. Vega system.")
            FileName = 'ubvri.dat'
            FilterNames = ['U','B','V','R','I']
            FilterZPs = [417.5e-11,632.0e-11,363.1e-11,217.7e-11,112.6e-11]
            loop=False
        elif choice==7:
            print("GAIA filters have been selected. Vega system.")
            FileName = 'gaia.dat'
            FilterNames = ['G   ','Gbp3','Grp3','Grvs']
            FilterZPs = [2.50386e-9,4.07852e-9,1.26902e-9,9.03937e-10]
            loop=False
            ## You can add your code or functions here
        else:
            # Any integer inputs other than values 1-7 we print(an error message
            filter_set = -1
            input("Wrong Filter set selection. Enter 'Enter' to try again..")
    return FileName,FilterNames,FilterZPs


##########################################################################

def ReadFilterCurves(FileName):
    """
    Reading Fiter Curves from an external ASCII text file
    """
    FilterCurve = []
#    FilePath = CheckFile(FileName)
    data0 = loadtxt(FileName, unpack=True,skiprows=0)

    WaveFilter = data0[0,:]
    Length=len(data0[0,:])
    i=0
    while True:
        try:
            FilterCurve.append(data0[i+1,:])
            i+=1
        except IndexError:
            NumberOfFilters=i
            break
    return WaveFilter,FilterCurve,NumberOfFilters

#########################################################################


def CheckFile(FileName,ScriptDir2):
    """
    Testing the presence of a file in the current or 2 other predefined folders
    """

    ScriptDir = str(Path.home())+'/.vn-astroutils/'
    Script2ndDir = ScriptDir2
    if os.path.isfile('./'+FileName):
        FullFileName = './'+FileName
    elif os.path.isfile(ScriptDir+FileName):
        FullFileName = ScriptDir+FileName
    elif os.path.isfile(Script2ndDir+FileName):
        FullFileName = Script2ndDir+FileName
    else:
        print('\nThe file',FileName,' is not found in any of the folders \n./')
        print(ScriptDir)
        print(Script2ndDir)
        print('Exit...')
        exit(-1)
    return FullFileName


#########################################################################


def print_header():
    print("")
    print("***************************************************************************************************")
    print("**                                      vn-spec-stat.py                                          **")
    print("**                                 ver. %s                                      **" % DateVer)
    print("**             A utility to calculate integrated fluxes in a wavelength range and                **")
    print("**                    different photometric filters from a list of spectra                       **")
    print("**                                     Vitaly  Neustroev                                         **")
    print("***************************************************************************************************")
    print()
    print("Usage: %s [options] @FileList/pattern [W1 [W2]]" % sys.argv[0])

def print_usage():
    """ Print the usage info """
    print("")
    print("Parameters:")
    print ("  @FileList is a list of spectra.")
    print ("    If no FileList is given then a list will be created from the files matching the pattern.")
    print ("    A filename given instead of FileList will be considered as a FileList consisting of 1 file.")
    print ("    Sometimes a pattern is not parsed correctly. Put it with single quotes, e.g. '*.txt'.")
    print("     The spectrum must have at least 2 columns: Wavelength and Flux")
    print("  W1: the start wavelength of the range for statistics.")
    print("  W2: the start wavelength of the range for statistics.")
    print ("Options:")
    print("  -f: the list of supported Filter sets")
    print("  -h: the program's history")
    print(" ")
    print("***************************************************************************************************")
    print(" ")



if __name__ == "__main__":
    print_header()
    FileList=False
    isTemplate = False
    s1 = 0
#    FiltName = 'V'

    if len(sys.argv) == 1:
        print_usage()

    for i in range(len(sys.argv)-1):
        CmdLinePar = sys.argv[i+1]
        if (CmdLinePar[0] != '-') and (CmdLinePar[0] != '+'):
            if not FileList:
                FileName = CmdLinePar
                if FileName[0] == '@':
                    s1 = 1
                    try:
                        infile = open(FileName[s1:], "r")
                        FileNames = infile.readlines()
                        infile.close()
                        print("The FileList ",FileName[s1:]," includes ",len(FileNames)," FileNames")
                        FileList=True
                        outfilename = FileName[s1:]
                    except:
                        print("Something wrong with the FileList ",FileName[s1:])
                elif TestFileList(FileName)==0:
                    try:
                        infile = open(FileName[s1:], "r")
                        FileNames = infile.readlines()
                        infile.close()
                        print("The FileList ",FileName[s1:]," includes ",len(FileNames)," FileNames")
                        FileList=True
                        outfilename = FileName[s1:]
                    except:
                        print("Something wrong with the FileList ",FileName[s1:])
                else:
                    isTemplate = True
                    FileNames=glob.glob(FileName)
                    if len(FileNames)>0:
                        print("\n",len(FileNames),"file(s) are found.")
                        FileList=True
                        outfilename = 'SpecStat'
                    else:
                        print("No file(s) matching the template are found.")
        elif CmdLinePar[0] == '-':
            if ('f' or 'F') in CmdLinePar[1:]:
                print_filter_sets()
                exit()
            if ('h') in CmdLinePar[1:]:
                print_history()
                exit()

    if (not FileList):
        FileName = input("Enter the file name of a list of spectra: ")
        FileName = FileName.rstrip('\n')
        if TestFileList(FileName)==0:
            try:
                infile = open(FileName[s1:], "r")
                FileNames = infile.readlines()
                infile.close()
                print("The FileList ",FileName[s1:]," includes ",len(FileNames)," FileNames")
                FileList=True
                outfilename = FileName[s1:]
            except:
                print("Something wrong with the FileList ",FileName[s1:])
                exit()
        else:
            isTemplate = True
            FileNames=glob.glob(FileName)
            if len(FileNames)>0:
                print("\n",len(FileNames),"file(s) are found.")
                FileList=True
                outfilename = 'SpecStat'
            else:
                print("No file(s) matching the template are found.")
                exit()
    NumOfSpec = len(FileNames)


if (len(sys.argv) > 2):
    W1 = float(sys.argv[2])
    print("The start wavelength of the range for statistics: ",W1)
else:
    print()
    W1 = input("Enter the start wavelength of the range for statistics: ")
    W1 = float(W1)
if (len(sys.argv) > 3):
    W2 = float(sys.argv[3])
    print("The end wavelength of the range for statistics:   ",W2)
else:
    W2 = input("Enter the end wavelength of the range for statistics:   ")
    W2 = float(W2)

if (len(sys.argv) > 4):
  FiltName = sys.argv[4]


NumOfFiles=len(FileNames)

outfile = open(outfilename+".Stat","w")
for i in range(NumOfFiles):
    FileName = FileNames[i].rstrip('\n')
    outfile.write(' %s \t %12.6e  \n' %  (FileName,SpecStat(FileName,W1,W2)))
outfile.close()
print('\n*** The integrated flux in the wavelength range is written in the file ',outfilename+".Stat")

FiltName,FilterNames,FilterZPs=filter_selection()
if FiltName:
    outfile2 = open(outfilename+".Phot","w")
    outfile2.write(' %s \t' %  ('#FileName'))
    for j in range(len(FilterNames)):
        outfile2.write(' %12s  \t' %  (FilterNames[j]))
    outfile2.write('\n')
    for i in range(NumOfFiles):
        FileName = FileNames[i].rstrip('\n')
        NumberOfFilters,SpecFlux = SpecPhotometry(FileName,FiltName)
        outfile2.write(' %s \t' %  (FileName))
        for j in range(NumberOfFilters):
            outfile2.write(' %12.6e  \t' %  (SpecFlux[j]))
        outfile2.write('\n')
    print("Number of Filters =",NumberOfFilters)
    outfile2.close()
    print('*** The flux in photometric bands is written in the file ',outfilename+".Phot")

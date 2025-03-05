#!/usr/bin/env python3
"""
A utility to calculate magnitudes and fluxes in different photometric filters from a spectrum
vn-spec2phot.py: ver. 20250305  (c) Vitaly Neustroev
"""
DateVer='2.0:  2025-March-05'
ScriptDir2 = '/scisoft/Other_Soft/Files4scripts/'

import sys
import os
import os.path
#import numpy as np
import math
from numpy import *
from math import log10
from sys import stdin
from numpy import trapz
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pathlib import Path


# ==============================================================================
def print_history():
    print('\nHistory:')
    print('2022-Dec-04: Location of files of transmission curves has been formalized.')
    print('2025-Mar-04: GAIA filters added.')
    print('2025-Mar-05: User interface improved.')

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

def Magnitude(VegaSystem,ZP,Flux):
    """
    Converting Fluxes to Magnitudes
    """
    if VegaSystem:
        Magnitude = -2.5*math.log10(Flux / ZP)
#        print("Vega System")
    else:
#        Magnitude = -2.5*math.log10(Flux) - 21.1
        Magnitude = -5.*math.log10(ZP*1.e-8) -2.5*math.log10(Flux) - 42.408
#        print("AB System")
    return Magnitude


#########################################################################

## Text menu in Python

def print_menu():       ## Your menu design here
    print("\n")
    print(43 * "-" , "Filter sets" , 43 * "-")
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


#########################################################################

def print_header():
    print("")
    print("***************************************************************************************************")
    print("**                                      vn-spec2phot.py                                          **")
    print("**                                 ver. %s                                      **" % DateVer)
    print("** A utility to calculate magnitudes and fluxes in different photometric filters from a spectrum **")
    print("**                                     Vitaly  Neustroev                                         **")
    print("***************************************************************************************************")
    print()
    print("Usage: %s [ -h | -f | -bNN ] FileName" % sys.argv[0])

def print_usage():
    """ Print the usage info """
    print("")
    print("FileName is a spectrum filename.")
    print("The file must have at least 2 columns: Wavelength and Flux")
    print(" ")
    print("Parameters:")
    print("  -bNN: a Filter set. NN must be integer from 1 to 7, for example -b5.")
    print("  -f: the list of supported Filter sets")
    print("  -h: the program's history")
    print(" ")
    print("***************************************************************************************************")
    print(" ")

##########################################################################

#########################################################################

print_header()
is_file = False
VegaSystem = True
loop=True

if len(sys.argv) == 1:
    print_usage()
    print("\nEnter the filename of a spectrum: ")
    output_file_path = stdin.readline()
    SpecFile = output_file_path.rstrip('\n')
    is_file = True
else:
    for i in range(len(sys.argv)-1):
        CmdLinePar = sys.argv[i+1]
        if CmdLinePar[0] == '-' or CmdLinePar[0] == '+':
            if ('h' or 'H') in CmdLinePar[1:]:
                print_history()
                exit()
            if ('f' or 'F') in CmdLinePar[1:]:
                print_filter_sets()
                exit()
            if CmdLinePar[0:2] == '-b':
                try:
                    filter_set = int(CmdLinePar[2:])
                except:
                    print('\nSomething wrong with the -b parameter. No filter set is recognised.')
                    filter_set = -1
        else:
            SpecFile = CmdLinePar
            is_file = True
if not is_file:
    print("\nEnter the filename of a spectrum: ")
    output_file_path = stdin.readline()
    SpecFile = output_file_path.rstrip('\n')
    is_file = True

if not os.path.exists(SpecFile):
    if not SpecFile=='':
        print("The file",SpecFile, "doesn't exist. Stop executing!")
    exit()
else:
    data1 = loadtxt(SpecFile, usecols=[0,1], unpack=True,skiprows=0)
    WaveSpec = data1[0,:]
    Spectrum = data1[1,:]
    try:
        data1 = loadtxt(SpecFile, usecols=[0,1], unpack=True,skiprows=0)
        WaveSpec = data1[0,:]
        Spectrum = data1[1,:]
    except:
        print("Something wrong with the file ",SpecFile)
        print("Exiting...")
        exit(-1)



while loop:          ## While loop which will keep going until loop = False
    if filter_set<0:
        print(30 * "-" , "MENU: select a filter set" , 30 * "-")
        print("0. Exit")
        print_menu()    ## Displays menu
        choice = input("Enter your choice [1-7]: ")
        try:
            choice = int(choice)
        except:
            continue
    else:
        choice = filter_set

    if choice==0:
        print("Exit")
        exit()
        ## You can add your code or functions here
    elif choice==1:
        print("UVOT filters have been selected")
        FileName = 'uvot.dat'
        FilterNames = ['w1','m1','w2']
        FilterZPs = [0,0,0]
        VegaSystem = False
        loop=False
        ## You can add your code or functions here
    elif choice==2:
        print("Johnson-Cousins UBVRI filters (Bessell) have been selected")
        FileName = 'bessell.dat'
        FilterNames = ['U','B','V','R','I']
        FilterZPs = [417.5e-11,632.0e-11,363.1e-11,217.7e-11,112.6e-11]
        loop=False
        ## You can add your code or functions here
    elif choice==3:
        print("PANSTARRS-PS1 filters have been selected")
        FileName = 'panstarrs1.dat'
        FilterNames = ['g','r','i','z','y']
        FilterZPs = [4810.,6170.,7520.,8660.,9620.]
        VegaSystem = False
        loop=False
        ## You can add your code or functions here
    elif choice==4:
        print("Sloan/SDSS filters have been selected")
        FileName = 'sdss.dat'
        FilterNames = ['u','g','r','i','z']
        FilterZPs = [3551,4686,6165,7481,8931]
        VegaSystem = False
        loop=False
        ## You can add your code or functions here
    elif choice==5:
        print("NIR J-H-Ks-K filters have been selected")
        FileName = 'jhksk.dat'
        FilterNames = ['J','H','Ks','K']
        FilterZPs = [31.47e-11,11.38e-11,3.961e-11,3.961e-11]
        loop=False
        ## You can add your code or functions here
        loop=False # This will make the while loop to end as not value of loop is set to False
    elif choice==6:
        print("Johnson-Cousins UBVRI filters (Astrodon) have been selected")
        FileName = 'ubvri.dat'
        FilterNames = ['U','B','V','R','I']
        FilterZPs = [417.5e-11,632.0e-11,363.1e-11,217.7e-11,112.6e-11]
        loop=False
    elif choice==7:
        print("GAIA filters have been selected")
        FileName = 'gaia.dat'
        FilterNames = ['G   ','Gbp3','Grp3','Grvs']
        FilterZPs = [2.50386e-9,4.07852e-9,1.26902e-9,9.03937e-10]
        loop=False
        ## You can add your code or functions here
    else:
        # Any integer inputs other than values 1-7 we print(an error message
        filter_set = -1
        input("Wrong Filter set selection. Enter 'Enter' to try again..")

FilePath = CheckFile(FileName,ScriptDir2)
WaveFilter,FilterCurve,NumberOfFilters = ReadFilterCurves(FilePath)
print("Number of Filters =",NumberOfFilters)

##if not os.path.exists(FileName):
##    print("The file", FileName, "doesn't exist")
##    exit()
##else:
##    WaveFilter,FilterCurve,NumberOfFilters = ReadFilterCurves(FileName)
##    print("Number of Filters =",NumberOfFilters)
###    print(WaveFilter[0],FilterCurve[0][0],FilterCurve[1][0],FilterCurve[2][0],FilterCurve[3][0],FilterCurve[4][0]


#print(WaveSpec[0], WaveSpec[-1], (WaveSpec[-1]-WaveSpec[0])/len(WaveSpec)

Wstep = math.ceil((WaveSpec[-1]-WaveSpec[0])/len(WaveSpec))
if Wstep<1.0:
    Wstep = 1.0
#print(Wstep)
Wstep = 1.0

if WaveSpec[0] > WaveFilter[0]:
    print("\nAttention! The object spectrum is only partially covered by some of filter curves.")
    print("Some (bluish) fluxes can be underestimated.")
    Wbegin = WaveSpec[0]
else:
    Wbegin = WaveFilter[0]

if WaveSpec[-1] < WaveFilter[-1]:
    print("\nAttention! The object spectrum is only partially covered by some of filter curves.")
    print("Some (reddish) fluxes can be underestimated.")
    Wend = WaveSpec[-1]
else:
    Wend = WaveFilter[-1]

Wnumber = int(math.floor((Wend-Wbegin)/Wstep))
#print("Wavelengths",Wbegin, Wend,"Number",Wnumber, (Wend-Wbegin)/Wstep
Wbegin+=((Wend-Wbegin)/Wstep-Wnumber)*Wstep/2
Wend=Wbegin+Wnumber*Wstep

#print("Wavelengths",Wbegin, Wend,"Number",Wnumber, Wstep

WaveNew = linspace(Wbegin, Wend, num=Wnumber+1, endpoint=True)
SpecNew = interp(WaveNew, WaveSpec,Spectrum)
##        FluxInterp = interp(WaveInterp, Wavelength,FluxNew)
        #WriteData2(len(WaveInterp),WaveInterp,FluxInterp,FileName+'.interp')

#print(SpecNew


# Plotting code
f, (ax1, ax2) = plt.subplots(2, figsize=(15,8))
gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
SpecMax=max(SpecNew)
SpecMin=min(SpecNew)
Wmax=max(WaveNew)
Wmin=min(WaveNew)
FiltMax=0


ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

#ax1.plot(WaveSpec, Spectrum, color="blue", lw=1.5, label="1 $\mathrm{\AA}\ $ Sampling")
ax1.plot(WaveNew, SpecNew, color="red", lw=1.5, label="The spectrum")

LamEff=[]
LamMean=[]
SpecFlux=[]
for i in range(NumberOfFilters):
    FilterNew = interp(WaveNew, WaveFilter[:], FilterCurve[i][:])
    SpecFilter = SpecNew * FilterNew
    SpecFilterLambda = SpecFilter * WaveNew
    FilterLambda = FilterNew * WaveNew
    FilterArea = trapz(FilterNew, dx=Wstep)
    FilterLambdaArea = trapz(FilterLambda, dx=Wstep)
    SpecFilterArea = trapz(SpecFilter, dx=Wstep)
    SpecFilterLambdaArea = trapz(SpecFilterLambda, dx=Wstep)
    LamEff.append(SpecFilterLambdaArea / SpecFilterArea)
    LamMean.append(FilterLambdaArea / FilterArea)
    SpecFlux.append(SpecFilterArea / FilterArea)
    print("Filter",i+1,":",FilterNames[i],"= %.3f"%(Magnitude(VegaSystem,FilterZPs[i],SpecFlux[i])),"mag   Leff=%.1f"%(LamEff[i]),"  Lmean=%.1f"%(LamMean[i]),"  Flux=%.6e"%(SpecFlux[i]))

    ax2.plot(WaveNew, FilterNew, color="red", lw=1.5, label="The spectrum")
    FiltMax=max(FiltMax,max(FilterNew))


ax1.plot(LamEff, SpecFlux, 'bo', markersize=9, label="Effective wavelength")
ax1.plot(LamMean, SpecFlux, 'ko', markersize=9, label="Mean wavelength")
#ax1.plot(LamEff, SpecFlux, color="black", lw=1.5, label="5 $\mathrm{\AA}\ $ Sampling")
#ax1.scatter(LamEff, SpecFlux, s, c="g", alpha=0.5, marker=r'$\clubsuit$', label="5 $\mathrm{\AA}\ $ Sampling")


#ax2.plot(WaveFilter[:], FilterCurve[NumberOfFilters-1][:], color="blue", lw=1.5, label="1 $\mathrm{\AA}\ $ Sampling")

#ax1.set_ylabel("Flux ($10^{-19}\ \mathrm{W/m^2/\\AA)}$", size=18)
ax1.set_ylabel("Flux ($\\mathrm{erg/cm^2/\\AA)}$", size=18)
ax2.set_ylabel("Filter curves", size=14)
ax2.set_xlabel("Wavelength ($\\mathrm{\\AA}$)", size=18)

ax1.set_xlim(Wmin, Wmax)
ax1.set_ylim(SpecMin, SpecMax)
ax2.set_xlim(Wmin, Wmax)
ax2.set_ylim(0.0, FiltMax)
ax1.legend()

plt.show()


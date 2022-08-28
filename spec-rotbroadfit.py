#!/usr/bin/python
#import numpy
from numpy import *
import numpy as np
import sys
from sys import stdin
#import numpy as np
from pylab import *
# Rotational Broadening

import pymc
# ... and now the funcFit package

from PyAstronomy import funcFit as fuf
from PyAstronomy import modelSuite as ms
from PyAstronomy import pyasl

import scipy.integrate as sci

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
    data0 = loadtxt(FileName, usecols=[0,1,2], unpack=True,skiprows=0)
    Wavelength=data0[0,:]
    Length=len(Wavelength)
    Flux=data0[1,:]
    FluxNew=data0[1,:]
    FluxErr=data0[2,:]

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



print("***********************************************************")
print(" ")
print("spec-rotbroadfit.py FileName FWHM/SpecPower Vbroad LimbDark")
print(" ")
print("FWHM/SpecPower: if this parameter < 100 then FWHM ")
print("                if > 100 then SpecPower")
print(" ")
print("***********************************************************")
print(" ")

#os.chdir("/home/benj/OwnCloud/My_Papers/BW_Scl/MgII/Final/AbsoluteFinal")


c = 299792.458
SpecLines     = [4481.126, 4481.150, 4481.325]     # MgII
#SpecLines     = [4481.13 , 4481.15 , 4481.33 ]     # MgII
LineTransProb = [0.7367  , -0.5643 , 0.5818  ]     # MgII
LineTransProb = 10**array(LineTransProb)
#print(LineTransProb)
#exit()

BroadEps = 0.35
if (len(sys.argv) > 1):
  FileName = sys.argv[1]
else:
  print(" ")
  print("Enter the spectrum filename: ")
  output_file_path = stdin.readline()
  FileName = output_file_path.rstrip('\n')

#if FileName[0]=='@':
    #FileList=True
    #infile = open(FileName[1:], "r")
    #lines = infile.readlines()
    #infile.close()
    #print 'The FileList ',FileName[1:],' includes ',len(lines),' FileNames'
#else:
    #FileList=False
    #lines = []
    #lines.append(FileName)
    #print 'The FileList includes ',len(lines),' FileName(s)'

if (len(sys.argv) > 2):
    FWHM = float(sys.argv[2])
else:
    FWHM = float(raw_input("Enter FWHM: "))

if FWHM>100.:
    SpecPower = True
    print("SpecPower = ", FWHM)
else:
    SpecPower = False
    print("FWHM = ", FWHM)

if (len(sys.argv) > 3):
    Vbroad = float(sys.argv[3])
else:
    Vbroad = float(raw_input("Enter vsini (km/s): "))
    BroadEps = float(raw_input("Enter limb darkening coefficient: "))

#infile = open(FileList, "r")
#lines = infile.readlines()
#infile.close()


#FileName=FileList
data0 = loadtxt(FileName, usecols=[0,1,2], unpack=True,skiprows=0)
Wavelength=data0[0,:]
Length=len(Wavelength)
Flux=data0[1,:]
FluxErr=data0[2,:]

MidW = (Wavelength[-1] - Wavelength[0]) / 2 + Wavelength[0]
dW = 3000 * MidW / c

InterW1=ceil((Wavelength[0] - dW)*100)/100
InterW2=floor((Wavelength[-1] +dW)*100)/100
InterPoints  = int((InterW2 - InterW1)*100+1)
WaveTemplate = linspace(InterW1, InterW2, InterPoints)
FluxTemplate = np.zeros_like(WaveTemplate)

for i in range(len(SpecLines)):
    index = find_nearest(WaveTemplate,SpecLines[i])
    FluxTemplate[index] = LineTransProb[i]

WriteData2(len(WaveTemplate),WaveTemplate,FluxTemplate,FileName+'.test')

#FluxTemplate = pyasl.instrBroadGaussFast(WaveTemplate,FluxTemplate,FWHM, edgeHandling="firstlast", maxsig=5.0)
#WriteData2(len(WaveTemplate),WaveTemplate,FluxTemplate,FileName+'.conv1')


#RotFlux = pyasl.rotBroad(WaveTemplate, FluxTemplate, BroadEps, Vbroad)
#WriteData2(len(WaveTemplate),WaveTemplate,RotFlux,FileName+'.rot')
#plot(WaveTemplate,FluxTemplate, "b-", WaveTemplate, RotFlux, "r-")
#show()

#if SpecPower:
    #SpecTemplate, avefwhm = pyasl.instrBroadGaussFast(WaveTemplate,RotFlux,FWHM,
                                        #edgeHandling="firstlast", fullout=True, maxsig=5.0)
    #print("FWHM used for the Gaussian kernel: ", avefwhm, " A")
#else:
    #SpecTemplate = pyasl.broadGaussFast(WaveTemplate,RotFlux,FWHM/2.35482, edgeHandling="firstlast", maxsig=5.0)

#WriteData2(len(WaveTemplate),WaveTemplate,SpecTemplate,FileName+'.conv')
#plot(WaveTemplate,SpecTemplate, "b-", WaveTemplate, RotFlux, "r-")
#show()



#Fit parameters:
#lineScale - A common scaling of the area of all lines.
#scale - A scaling of the entire spectrum.
#eps - Linear limb-darkening coefficient.
#vrad - Radial velocity [km/s].
#vsini - Projected stellar rotational velocity [km/s].
#A{n} - The amplitudes (area) parameters of the individual Gaussians.
#sig{n} - The standard deviations of the individual Gaussian.
#mu{n} - The position of the individual Gaussians.

# Create our line list with 4 line
#lineList = np.zeros((4, 3))
lineList = np.zeros((3, 2))
# Assign wavelengths (in A)
lineList[0, 0] = 4481.13
lineList[1, 0] = 4481.15
lineList[2, 0] = 4481.33
#lineList[3, 0] = 5007.64
# Assign EWs (in A)
lineList[0, 1] = 5.45380995
lineList[1, 1] = 0.27270933
lineList[2, 1] = 3.81768419
#lineList[3, 1] = 0.12
# Assign depths (0-1)
#lineList[0, 2] = 0.97
#lineList[1, 2] = 0.9
#lineList[2, 2] = 0.99
#lineList[3, 2] = 0.35

#wvl = np.arange(5000., 5010., 0.01)
wvl = Wavelength
#Flux
# Get an instance of the LLGauss class
avefwhm=0.821933
avefwhm=0.98
avefwhm=0.1
llg = ms.LLGauss(lineList,uniSig=avefwhm)
# ... and define some starting value
#llg["xmax"] = 7.0
llg["eps"] = BroadEps
llg["vsini"] = 150
llg["vrad"] = 47.0
llg["lineScale"] = 6.0
#llg["mu"] = 4482.28
#llg["gsig"] = 7200.0

# Have a look at the model parameters
#llg.parameterSummary()

# Use all line strengths for fitting
#llg.thawLineStrengths()
#llg.thaw(["lineScale", "vsini"])
#llg.thaw(["lineScale", "scale", "vsini"])
llg.thaw(["lineScale", "scale", "vrad","vsini"])
#llg.thaw(["lineScale", "scale", "vrad","vsini", "eps"])


llg.parameterSummary()
print(llg.freeParameters())
#lims = {"lineScale": [0.40, 0.46], "scale": [0.95, 1.05], "vrad": [60., 80.], "vsini": [100., 500.], "eps": [0.1, 0.9]}
#llg.setRestriction({"scale": [0.95, 1.05], "vrad": [60., 80.], "vsini": [100., 500.], "eps": [0.1, 0.9]})
#print(llg.getRestrictions())
#exit()
# and fit
llg.fit(Wavelength, Flux)

# ... and show the resulting parameter values.
print("... and show the resulting parameter values.")

llg.parameterSummary()

# Plot the result
#plt.subplot(2, 1, 2)
#plt.errorbar(Wavelength, Flux, yerr=np.ones(len(wvl))*0.01, fmt='bp')
#plt.plot(Wavelength, llg.evaluate(Wavelength), 'r--')
#plt.show()

# Plot the data and the model
WriteData2(len(Wavelength),Wavelength,llg.model,FileName+".fit")
plt.plot(Wavelength, Flux, 'b--')
plt.plot(Wavelength, llg.model, 'r--')
plt.show()


X0 = llg.freeParameters()
print("Starting point for sampling: ", X0)

llg.freeze("scale")
llg.freeze("lineScale")
llg.thaw(["vrad","vsini","eps"])

llg.setRestriction({"vrad": [0., 120.], "vsini": [50., 500.], "eps": [0.1, 0.9]})
llg.fit(Wavelength, Flux)
print("... and show the resulting parameter values.")
llg.parameterSummary()
plt.plot(Wavelength, Flux, 'b--')
plt.plot(Wavelength, llg.model, 'r--')
plt.show()


print("\nHere was the previous exit\n")
#exit()

# Now we specify the limits within which the individual parameters
# can be varied. Actually, you specify the limits of uniform priors
# here.
lims = {"lineScale": [0.040, 0.046], "scale": [0.95, 1.05], "vrad": [40., 120.], "vsini": [100., 500.]}
steps = {"lineScale": 0.005, "scale": 0.005, "vrad": 0.5, "vsini": 10}

# Start the sampling. The resulting Marchov-Chain will be written
# to the file 'mcmcTA.tmp'. In default configuration, pickle
# is used to write that file.
# To save the chain to a compressed 'hdf5'
# file, you have to specify the dbArgs keyword; e.g., use:
#   dbArgs = {"db":"hdf5", "dbname":"mcmcExample.hdf5"}

llg["eps"] = BroadEps
llg.freeze("eps")
llg.thaw(["lineScale", "scale", "vrad","vsini"])
X0 = llg.freeParameters()
print("Starting point for sampling: ", X0)
llg.parameterSummary()

lims = {"lineScale": [0.010, 0.050], "scale": [0.95, 1.05], "vrad": [20., 70.], "vsini": [10., 500.]}
steps = {"lineScale": 0.005, "scale": 0.005, "vrad": 0.5, "vsini": 10}

llg.fitMCMC(Wavelength, Flux, X0, lims, steps, yerr=FluxErr,
           iter=3000, burn=500, thin=1,
           dbfile="mcmcTA.emcee")

exit()            ########################################################################

priors = {"vrad":fuf.FuFPrior("limuniform", lower=0.0, upper=100.), \
"vsini":fuf.FuFPrior("limuniform", lower=0.0, upper=100.), \
"lineScale":fuf.FuFPrior("limuniform", lower=0.0, upper=100.), \
"scale":fuf.FuFPrior("limuniform", lower=0.0, upper=100.)}

# Note that the filename should end in .emcee. Substitute this filename
# in the following examples.
llg.fitEMCEE(Wavelength, Flux, yerr=FluxErr, sampleArgs={"iters":2500}, \
dbfile="mcmcTA.emcee", priors=priors)

# Create an instance of TraceAnalysis
# telling it which file to use
#ta = fuf.TraceAnalysis("mcmcTA.tmp")
ta = fuf.TraceAnalysis("mcmcTA.emcee")

# Have a look at the deviance to check if and when
# the chains reached equilibrium.
ta.plotTrace("deviance")
ta.show()

# Say, we are sure that after 500 iterations, the chain
# reached equilibrium. We use this as the burn-in phase
ta.setBurn(500)

# Have a second look at the deviance, this time considering
# the burn-in. Note that the first 500 iterations are not
# removed from the chain. They are just not considered any
# more.
ta.plotTrace("deviance")
ta.show()


# Create an instance of TraceAnalysis
# telling it which file to use
ta = fuf.TraceAnalysis("mcmcTA.emcee")

# Use the burn-in from the previous example
ta.setBurn(500)

# See which model parameters have been sampled
print("Available parameters: ", ta.availableParameters())

# Access the traces of these parameters
print("Trace for vrad: ", ta["vrad"])

# Access the traces of these parameters
print("Trace for vsini: ", ta["vsini"])

# Calculate mean, median, standard deviation, and
# credibility interval for the available parameters
for p in ta.availableParameters():
  hpd = ta.hpd(p, cred=0.95)
  print("Parameter %5s, mean = % g, median = % g, std = % g, 95%% HPD = % g - % g" \
        % (p, ta.mean(p), ta.median(p), ta.std(p), hpd[0], hpd[1]))




exit()
############################################

#xmax - Maximal extent of the profile
#eps  - Linear limb-darkening coefficient
#A    - Area under the profile (negative for absorption)
#off  - An offset applied to the profile
#lin  - Gradient of a linear term to adjust the 'continuum'
#mu   - Center of the profile (same units as xmax)
#gsig - The standard deviation of a Gaussian with which the
        #rotational profile is convoluted, e.g., to model instrumental resolution.




# Get an instance of the model ...
x = ms.RotBroadProfile()
# ... and define some starting value
x["xmax"] = 7.0
x["eps"] = BroadEps
x["A"] = -0.42
x["off"] = 0.0
x["lin"] = 0.0
x["mu"] = 4482.28
x["gsig"] = 7200.0

# Define a radial velocity axis
#vv = np.linspace(-90., 90., 200)
#Wavelength

# Construct some "data" and ...
#data = x.evaluate(Wavelength)
# ... add noise
#data += np.random.normal(0.0, 1e-3, data.size)

# Fit the model using A, xmax, and eps as free
# parameters ...
x.thaw(["A", "xmax", "off","lin"])
x.fit(Wavelength, Flux)
# ... and show the resulting parameter values.
x.parameterSummary()

# Plot the data and the model
plt.plot(Wavelength, Flux, 'b--')
plt.plot(Wavelength, x.model, 'r--')
plt.show()


exit()
############################################

WaveTemplate = linspace(InterW1, InterW2, InterPoints)


if Interp:
    WaveInterp = linspace(InterW1, InterW2, InterPoints)
    FluxInterp = interp(WaveInterp, Wavelength,Flux)
#        FluxInterp = interp(WaveInterp, Wavelength,FluxNew)
    WriteData2(len(WaveInterp),WaveInterp,FluxInterp,FileName+'.interp')
else:
    WaveInterp = Wavelength
    FluxInterp = Flux

if Rebin:
    Xn,Yn = Rebin2Col(WaveInterp,FluxInterp,FileName+'.bin')
else:
    Xn = WaveInterp
    Yn = FluxInterp
    Xnm = WaveInterp
    Ynm = FluxInterp

FluxMed=MedFilt(FluxInterp,Median3)
WriteData2(len(WaveInterp),WaveInterp,FluxMed,FileName+'.median')

if Rebin:
    Xnm,Ynm = Rebin2Col(WaveInterp,FluxMed,FileName+'.medbin')

if RotBroad:
    RotFlux = pyasl.rotBroad(Xnm, Ynm, BroadEps, Vbroad)
    WriteData2(len(Xnm),Xnm,RotFlux,FileName+'.rot')
    plot(Xnm, Ynm, "b-", Xnm, RotFlux, "r-")
    show()

plot(Wavelength, Flux, "b-", WaveInterp, FluxMed, "r-", Xnm, Ynm, "k-")
show()


Yerr = [0] * len(Xnm)
Edge = (Median3-1)/2
#print("Edge=",Edge)
diff4std=np.array(Yn) - np.array(FluxMed)
#list(Yn-Ynm)


for i in range(len(Xnm)):
    if (i > Edge-1) and (i < len(Xnm)-Edge):
        Yerr[i] = std(diff4std[(i-Edge):(i+Edge)])
for i in range(Edge):
    Yerr[i] = Yerr[Edge]
    Yerr[-1-i] = Yerr[-1-Edge]
        
WriteData3(len(Xn),Xn,Yn,Yerr,FileName+'.stddev')

exit()


for line in lines:
    FileName = line.rstrip('\n')
    CleanSpec(FileName)

#!/usr/bin/python
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    #!/usr/bin/env python
""" Fast algorithm for spectral analysis of unevenly sampled data

The Lomb-Scargle method performs spectral analysis on unevenly sampled
data and is known to be a powerful way to find, and test the
significance of, weak periodic signals. The method has previously been
thought to be 'slow', requiring of order 10(2)N(2) operations to analyze
N data points. We show that Fast Fourier Transforms (FFTs) can be used
in a novel way to make the computation of order 10(2)N log N. Despite
its use of the FFT, the algorithm is in no way equivalent to
conventional FFT periodogram analysis.

Keywords:
  DATA SAMPLING, FAST FOURIER TRANSFORMATIONS,
  SPECTRUM ANALYSIS, SIGNAL  PROCESSING

Example:
  > import numpy
  > import lomb
  > x = numpy.arange(10)
  > y = numpy.sin(x)
  > fx,fy, nout, jmax, prob = lomb.fasper(x,y, 6., 6.)

Reference:
  Press, W. H. & Rybicki, G. B. 1989
  ApJ vol. 338, p. 277-280.
  Fast algorithm for spectral analysis of unevenly sampled data
  bib code: 1989ApJ...338..277P

"""
from numpy import *
from numpy.fft import *
import sys
from math import ceil
from math import log10
from sys import stdin
import numpy
import matplotlib.pylab as mpl
from PyAstronomy.pyTiming import pyPeriod
from PyAstronomy.pyTiming import pyPDM

import matplotlib.pylab as plt


#import numpy as np
from scipy import signal
#import matplotlib.pyplot as plt

##########################################################################

def WriteSpec(nn,aa,bb,output_file_path):
    """
    Write two columns of data to an external ASCII text file
    """
    output_file = output_file_path.rstrip('\n')
    outfile = open(output_file,"w")
    for i in range (0, nn):
        outfile.write(' %15.7f \t %12.6f \n' %  (aa[i],bb[i]))
    outfile.close()

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

##########################################################################


def WriteSlidogram(num_spec,num_freq,DataList,output_file_path):
    """
    Write the matrix of data to an external ASCII text file
    """
    output_file = output_file_path.rstrip('\n')
    outfile = open(output_file,"w")
    for i in range (0, num_spec):
        for j in range (0, num_freq):
            outfile.write('%8.5e ' %  (DataList[i][j]))
        outfile.write('\n')
#        outfile.write(' %12.6f \t %12.6f \n' %  (aa[i],bb[i]))
    outfile.close()



#######################################

def enter_float():
    """
    Enter a floating point number and check its validity
    """
    number_float = None
    while not number_float:
        try:
            s =stdin.readline()
            number_float = float(s)
            if number_float == 0:
                break
        except ValueError:
            print("Invalid Number.  Enter float number. ")
    return number_float

########################################################################

def enter_int():
    """
    Enter an integer and check its validity
    """
    number_int = None
    while not number_int:
        try:
            s =stdin.readline()
            number_int = int(s)
        except ValueError:
            print("Invalid Number.  Enter integer number. ")
    return number_int

########################################################################


print("*****************************************")
print(" ")
print("myperiod.py infile outfile Ofac Fmax [Fmin ['W'/'D']]")
print("      Ofac: Oversampling factor")
print("      FMax: Maximum frequency")
print("      FMin: Minimum frequency")
print("      W[N]: Dewhitening the data using  the best sine frequency")
print("                                        and repeat calculations")
print("            If N (from 1 to 9) is given then N best frequencies")
print("                              will be removed one after another")
print("      D[T]: Calculation of  a dynamical  (sliding)  periodogram")
print("            If 'T' is given  then calculations will be based on")
print("                   a time window otherwise on an index interval")
print(" ")
print("*****************************************")
print(" ")
if (len(sys.argv) > 1):
  FileName = sys.argv[1]
else:
  #print " "
  print("\nEnter the input light curve filename: ")
  output_file_path = stdin.readline()
  FileName = output_file_path.rstrip('\n')

if (len(sys.argv) > 2):
    FileNameOutput = sys.argv[2]
else:
    #print(" ")
    print("\nEnter the output power spectrum filename: ")
    output_file_path = stdin.readline()
    FileNameOutput = output_file_path.rstrip('\n')

data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
Time=data0[0,:]
Flux=data0[1,:]
Npoints=len(Flux)

dTime = []
for i in range(Npoints-1):
    dTime.append(Time[i+1]-Time[i])
print('\nAverage time interval is',mean(dTime),'+/-',std(dTime))
dTmin = max(dTime)
dTime2 = set(dTime)
dTime2.remove(dTmin)
dTmin2 = max(dTime2)
print('Largest time interval is',dTmin)
print('Second largest time interval is',dTmin2)


# Build the TimeSeries instance
lc = pyPeriod.TimeSeries(Time, Flux)
HifacDef=lc.returnNyquist()

Ofac = 10
Fstep = 1 / (max(Time)-min(Time))
print("\nDefault frequency step is", Fstep)

if (len(sys.argv) > 3):
    try:
      Fstep = float(sys.argv[3])
      if Fstep<0:
          Ofac = int(abs(Fstep))
          print("Oversampling factor Ofac is", Ofac)
          Fst = False
      else:
          dF=Fstep
          print("Frequency step is", dF)
          Fst = True
      #Ofac = int(sys.argv[3])
      #Ofac = float(sys.argv[3])
      #print("Oversampling factor Ofac is", Ofac)
    except ValueError:
      print("Invalid Ofac. Enter Oversampling factor (negative integer) or Frequency step (positive number): ")
      #Ofac=enter_float()
      Fstep = enter_float()
      #Ofac=enter_int()
      if Fstep<0:
          Ofac = int(abs(Fstep))
          print("Oversampling factor Ofac is", Ofac)
          Fst = False
      else:
          dF=Fstep
          print("Frequency step is", dF)
          Fst = True
      #print("Oversampling factor Ofac is", Ofac)
else:
    #print(" ")
    print("\nDefault oversampling factor Ofac is", Ofac)
    Fst = False
# if (len(sys.argv) > 4):
#     try:
#       Hifac = float(sys.argv[4])
#     except ValueError:
#       print 'Invalid Hifac. Enter the highest frequency: '
#       Hifac=enter_float()
# else:
#     print " "
#     Hifac = 0.
#     print "Default Hifac is average Nyquist frequency"


########################################################################

if (len(sys.argv) > 4):
    try:
      Hifac = float(sys.argv[4])
      Fmax = Hifac
      print("Fmax is ", Hifac)
      Hifac = Hifac / HifacDef
    except ValueError:
      print("Invalid Fmax. Enter the highest frequency: ")
      Hifac=enter_float()
      Fmax = Hifac
      print("Hifac is ", Hifac)
      Hifac = Hifac / HifacDef
else:
    #print " "
    Hifac = 1.
    Fmax = HifacDef
    print("\nDefault Fmax is average Nyquist frequency", HifacDef)


if (len(sys.argv) > 5):
    try:
      Fmin = float(sys.argv[5])
      print("Fmin is ", Fmin)
      Hifac = Hifac / HifacDef
    except ValueError:
      print("Invalid Fmin. Enter the lowest frequency: ")
      Fmin=enter_float()
      print("Fmin is ", Fmin)
else:
    Fmin = 0.
    print("\nDefault Fmin is 0. It can be changed using command line.")

IsDeWhite = False
IsDyn = False
isDynTime = False
if (len(sys.argv) > 6):
    try:
      DeWhite = sys.argv[6]
      if DeWhite[0] == 'W' or DeWhite[0] == 'w':
          IsDeWhite = True
          print("LightCurve will be dewhitened and a new periodogram will be calculated")
          try:
              WN=int(DeWhite[1])
          except ValueError:
              WN=1
      elif DeWhite[0] == 'D' or DeWhite[0] == 'd':
          IsDyn = True          
          try:
              if DeWhite[1] == 'T' or DeWhite[1] == 't':
                  isDynTime = True
          except ValueError:
              isDynTime = False
      else:
          print("\nLightCurve will NOT be dewhitened")
    except ValueError:
      print("\nLightCurve will NOT be dewhitened")
else:
    print("\nLightCurve will NOT be dewhitened")

print("\n")

###############

#fs = 10e3
#N = 1e5
#amp = 2 * np.sqrt(2)
#noise_power = 0.01 * fs / 2
#time = np.arange(N) / float(fs)
#mod = 500*np.cos(2*np.pi*0.25*time)
#carrier = amp * np.sin(2*np.pi*3e3*time + mod)
#noise = rng.normal(scale=np.sqrt(noise_power),size=time.shape)
#noise *= np.exp(-time/5)
#x = carrier + noise


# f, t, Zxx = signal.stft(Flux, 1/Fstep, nperseg=86400)
# plt.pcolormesh(t, f, numpy.abs(Zxx), vmin=0, vmax=1, shading='gouraud')
# plt.title('STFT Magnitude')
# plt.ylabel('Frequency [Hz]')
# plt.xlabel('Time [sec]')
# plt.show()

# exit()

# fft = pyPeriod.Fourier(lc)
# fig, ax = fft.plot()
# print('Mean power level:', numpy.mean(fft.power))
# plt.show()

# NoutFFT=len(fft.freq)
# WriteSpec(NoutFFT,fft.freq,fft.power,FileNameOutput+'.fft.dat')

# exit()


if Fst:
    Freq = numpy.arange(Fmin,Fmax,dF)
    gls = pyPeriod.Gls(lc, fbeg=Fmin, fend=Fmax, ofac=Ofac, freq=Freq, norm='ZK', verbose=True)
else:
    gls = pyPeriod.Gls(lc, fbeg=Fmin, fend=Fmax, ofac=Ofac, hifac=Hifac, norm='ZK', verbose=True)
# Print helpful information to screen
#gls.info()

#print(gls.amp, gls.fbest)
#maxPower = gls.pmax
#print("GLS maximum power: ", maxPower)
#print("GLS statistics of maximum power peak: ", gls.stats(maxPower))


Nout=len(gls.freq)
WriteSpec(Nout,gls.freq,gls.power,FileNameOutput+'.GLombScargle.dat')
print('File'+FileNameOutput+'.GLombScargle.dat is created')
gls.plot(block=True)

if IsDeWhite:
    FluxWhite = Flux
    for i in range(WN):
        print("\nDewhitening...",i+1,"...\n")
        FluxWhite = FluxWhite - gls.sinmod(Time)
        lcW = pyPeriod.TimeSeries(Time, FluxWhite)
        if Fst:
            Freq = numpy.arange(Fmin,Fmax,dF)
            gls = pyPeriod.Gls(lcW, fbeg=Fmin, fend=Fmax, ofac=Ofac, freq=Freq, norm='ZK', verbose=True)
        else:
            gls = pyPeriod.Gls(lcW, fbeg=Fmin, fend=Fmax, ofac=Ofac, hifac=Hifac, norm='ZK', verbose=True)

        WriteSpec(Nout,gls.freq,gls.power,FileNameOutput+'.DeWhite'+str(i+1)+'.GLombScargle.dat')
        WriteSpec(len(Time),Time,FluxWhite,FileNameOutput+'.DeWhite'+str(i+1)+'.dat')
        gls.plot(block=True)
        print('File'+FileNameOutput+'.DeWhite'+str(i+1)+'.GLombScargle.dat is created')

if IsDyn:
    Slidogram = []
    TimeW = []
    FluxW = []
    print('\n\nYour time-series consists of',Npoints,'points covering','{:.5f}'.format(Time[-1]-Time[0]),'time units.')
    if isDynTime:
        print('Enter the window length (a time interval which will be used for calculations of single periodograms):')
        Nwin = enter_float()
        while Nwin<dTmin:
            print('Enter at least',dTmin,':')
            Nwin = enter_float()
        NcoverMin = ceil((Time[-1]-Time[0])/Nwin)
    else:
        print('Enter the window length (how many points will be used for calculations of single periodograms):')
        Nwin = enter_int()
        NcoverMin = ceil(Npoints/Nwin)
    print('Enter the number of windows per total time-series (at least',NcoverMin,'):')
    Ncover = enter_int()
    while Ncover<NcoverMin:
        print('Enter at least',NcoverMin,':')
        Ncover = enter_int()
    if isDynTime:
        Delta = (Time[-1]-Time[0] - Nwin)/(Ncover-1)
        print('Delta=',Delta)
        for Ncur in range(Ncover):
            print('Power spectrum ',Ncur+1)
            TimeWin=Time[find_nearest(Time,Delta*Ncur+Time[0]):find_nearest(Time,Nwin+Delta*Ncur+Time[0])]
            FluxWin=Flux[find_nearest(Time,Delta*Ncur+Time[0]):find_nearest(Time,Nwin+Delta*Ncur+Time[0])]
#            print(find_nearest(Time,Delta*Ncur+Time[0]),find_nearest(Time,Nwin+Delta*Ncur+Time[0]),find_nearest(Time,Nwin+Delta*Ncur+Time[0])-find_nearest(Time,Delta*Ncur+Time[0]))
            lcWin = pyPeriod.TimeSeries(TimeWin, FluxWin)            
            #print(average(TimeWin),average(FluxWin))
            TimeW.append(average(TimeWin))
            FluxW.append(average(FluxWin))
            if Fst:
                Freq = numpy.arange(Fmin,Fmax,dF)
                gls = pyPeriod.Gls(lcWin, fbeg=Fmin, fend=Fmax, ofac=Ofac, freq=Freq, norm='ZK')
            else:
                gls = pyPeriod.Gls(lcWin, fbeg=Fmin, fend=Fmax, ofac=Ofac, hifac=Hifac, norm='ZK')
            Slidogram.append(gls.power)
    else:
        Delta = int(round((Npoints - Nwin)/(Ncover-1)))
        for Ncur in range(Ncover):
            print('Power spectrum ',Ncur+1)
            if Ncur < Ncover/2:
                TimeWin=Time[Delta*Ncur:Nwin+Delta*Ncur-1]
                FluxWin=Flux[Delta*Ncur:Nwin+Delta*Ncur-1]
    #            print(Delta*Ncur+1,Nwin+Delta*Ncur,Nwin)
            else:
                TimeWin=Time[-(Ncover-Ncur-1)*Delta-Nwin:-1-(Ncover-Ncur-1)*Delta]
                FluxWin=Flux[-(Ncover-Ncur-1)*Delta-Nwin:-1-(Ncover-Ncur-1)*Delta]
    #            print(-(Ncover-Ncur-1)*Delta-Nwin+Npoints,-1-(Ncover-Ncur-1)*Delta+Npoints,Nwin)
            lcWin = pyPeriod.TimeSeries(TimeWin, FluxWin)
            #print(average(TimeWin),average(FluxWin))
            TimeW.append(average(TimeWin))
            FluxW.append(average(FluxWin))
            if Fst:
                Freq = numpy.arange(Fmin,Fmax,dF)
                gls = pyPeriod.Gls(lcWin, fbeg=Fmin, fend=Fmax, ofac=Ofac, freq=Freq, norm='ZK')
            else:
                gls = pyPeriod.Gls(lcWin, fbeg=Fmin, fend=Fmax, ofac=Ofac, hifac=Hifac, norm='ZK')
            Slidogram.append(gls.power)
    WriteSpec(Ncover,TimeW, FluxW,FileNameOutput+'.DynLightCurve.dat')
#    print(TimeW, FluxW)
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.plot(TimeW, FluxW, 'b.')
    plt.show()    
    
    plt.xlabel("Frequency")
    plt.ylabel("Power")
#    plt.pcolormesh(TimeW, FluxW, numpy.abs(Slidogram), vmin=0, vmax=1, shading='gouraud')
    for Ncur in range(Ncover):
        plt.plot(gls.freq, Slidogram[Ncur], 'b-')
    plt.show()
    WriteSlidogram(Ncover,len(Slidogram[Ncur]),Slidogram,FileNameOutput+'.Slidogram.dat')


# and plot power vs. frequency.
# plt.xlabel("Frequency")
# plt.ylabel("Power")
# plt.plot(gls.freq, gls.power, 'b.-')
# plt.show()


########################################
# ls = pyPeriod.LombScargle(lc, Ofac, Hifac)
# Nout=len(ls.freq)
# WriteSpec(Nout,ls.freq,ls.power,FileNameOutput+'.scargle.dat')
########################################

# fft = pyPeriod.Fourier(lc)
# fig, ax = fft.plot()
# print('Mean power level:', numpy.mean(fft.power))
# plt.show()

# NoutFFT=len(fft.freq)
# WriteSpec(NoutFFT,fft.freq,fft.power,FileNameOutput+'.fft.dat')

########################################

#S = pyPDM.Scanner(minVal=gls.freq[0], maxVal=gls.freq[Nout-1], dVal=gls.freq[1]-gls.freq[0], mode="frequency")
#P = pyPDM.PyPDM(Time, Flux)
#f1, t1 = P.pdmEquiBinCover(10, 3, S)
#WriteSpec(len(f1),f1,t1,FileNameOutput+'.pdm.dat')

# Show the result
#plt.figure(facecolor='white')
#plt.title("Result of PDM analysis")
#plt.xlabel("Frequency")
#plt.ylabel("Theta")
#plt.plot(f1, t1, 'bp-')
##plt.plot(f2, t2, 'gp-')
##plt.legend(["pdmEquiBinCover", "pdmEquiBin"])
#plt.show()



########################################################################

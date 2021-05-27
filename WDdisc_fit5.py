#!/usr/bin/python3
"""WDdisc_fit: ver. 20210118  (c) Vitaly Neustroev"""

import sys
import os
import os.path
#import numpy as np
from numpy import * 
from sys import stdin 
from numpy import trapz
import scipy.optimize as opt
from scipy.interpolate import splrep,splev
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.axes as pltaxes
from matplotlib.gridspec import GridSpec
import operator

#from multiprocessing import cpu_count
#https://www.techbeamers.com/python-multithreading-concepts/
#https://python-forum.io/Thread-Multithreading-in-a-loop
#https://stackoverflow.com/questions/38856172/simple-multithread-for-loop-in-python
#https://docs.python.org/2.7/library/multiprocessing.html#programming-guidelines
#https://realpython.com/intro-to-python-threading/
#https://towardsdatascience.com/10x-faster-parallel-python-without-python-multiprocessing-e5017c93cce1


##################################################################

def printHeader():
    """ Print the Header """
    print()
    print("*******************************************************************************")
    print("                WDdisc_fit: ver. 20210118  (c) Vitaly Neustroev")

def printUsage():
    """ Print the usage info """
    print("Usage:")
    print("WDdisc_fit.py --cut ObsSpectrum VelocityCut [He_lines: yes/no]")
    print("WDdisc_fit.py FitModel ObsSpectrum FileList_of_SpecModels [ResultDir]")
    print("WDdisc_fit.py FitModel ObsSpectrum FileList_of_SpecModels ResultDir Iterations")
    print()
    print("FitModels:")
    print("  1: WD only ")
    print("  2: WD + Power-Law ")
    print("  3: WD only. Normalized, only lines. The WD model spectra are normalized automatically.")
    print("  4: WD + Power-Law, automatically normalized, only lines.")
    print("  5: WD + a Fixed model (e.g. Slab). A model must be in the file 'FixedModel.dat' ")
    print("  Negative number: WD + Fixed index PL. The PL index is equal to this number.")
    print()    
    print("  Experimental models: ")
    print("    21: WD at the given distance + Power-Law ")
    print("        The distance will be asked and the WD spectrum will be scaled accordingly.")
    print("        FileList_of_SpecModels must have 3 columns: FileName_of_WDmodel  Temperature  logg*100.")
    print("    0: WD + only (Experimental!!!): ")
    print("    9: Normalized WD only (Experimental!!!): ")
    print("       The WD model spectra are normalized automatically.")
    print("FileList_of_SpecModels: List of Files with WD models.")
    print("                        For FitModel=21 it must have 3 columns: FileName  WDtemp  WDlogg*100.")     
    print()
    print("*******************************************************************************")
    print()


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

##################################################################

def WriteDataModel(nn,aa,bb,cc,dd,ee,ff,output_file_path):
    """
    Write two columns of data to an external ASCII text file
    """
    output_file = output_file_path.rstrip('\n')
    outfile = open(output_file,"w")
    outfile.write('#Wavelength  \t TotalModel  \t ObsSpec     \t WD          \t PowerLaw  \t Residuals \n')
    for i in range (0, nn):
        outfile.write(' %12.6f \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e  \n' %  
		(aa[i],bb[i],cc[i],dd[i],ee[i],ff[i]))
#        outfile.write(' %12.6f \t %12.6f \n' %  (aa[i],bb[i]))
    outfile.close()
########################################################################

def WriteData3(nn,aa,bb,cc,output_file_path):
    """
    Write three columns of data to an external ASCII text file
    """
    outfile = open(output_file_path,"w")
    for i in range (0, nn):
        outfile.write(' %s \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        #outfile.write(' %12.6f \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close() 

##########################################################################

def WriteData5(nn,aa,bb,cc,dd,ee,output_file_path):
    """
    Write five columns of data to an external ASCII text file
    """
    #output_file = output_file_path.rstrip('\n')
    outfile = open(output_file_path,"w")   
    for i in range (0, nn):
        outfile.write(' %s \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i],dd[i],ee[i]))
        #outfile.write(' %12.6f \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close() 

##########################################################################

def WriteData6(nn,aa,bb,cc,dd,ee,ff,output_file_path):
    """
    Write five columns of data to an external ASCII text file
    """
    #output_file = output_file_path.rstrip('\n')
    outfile = open(output_file_path,"w")   
    for i in range (0, nn):
        outfile.write(' %s \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i],dd[i],ee[i]),ff[i])
        #outfile.write(' %12.6f \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close() 

##########################################################################


def errfunc1(p, isErr, a1, a2, a3):
    N=len(a1)
    #print "p=",p
    Results = calc_Residual(N, isErr, a1, a2*abs(p), a3)
    return Results

#def calc_Residual1(N, isErr, a1, a2, a3):
    #""" calculate the residual between a1 and a2 (a3 is the error) """
    #if isErr:
        #Chi2 = sqrt(sum(((a1-a2)/a3)**2) / (N-1))
    #else:
        #Chi2 = sqrt(sum((a1-a2)**2) / (N-1))
    ##Chi2 = sqrt(sum((a1-a2)**2)) / (N-1)
    ##print "Chi2=",Chi2
    #return Chi2

##########################################################################

def errfunc2(p, isErr, x1, y1, y2, y3):
    #isErr,WaveSpec,Spectrum+SpecAddedNoise, SpecModelNew,SpecErr
    N=len(x1)
    powerlaw = list(x1)
    powerlaw = abs(y2*p[0])+abs(p[1])*x1**(-abs(p[2]))
    #for i in range(N):
        #powerlaw[i] = y2[i]*p[0]+p[1]*x1[i]**p[2]
    Results = calc_Residual(N, isErr, y1, powerlaw, y3)
    return Results

##########################################################################

def errfunc21(p, isErr, x1, y1, y2, y3, WDnorm):
    #isErr,WaveSpec,Spectrum+SpecAddedNoise, SpecModelNew,SpecErr
    N=len(x1)
    powerlaw = list(x1)
    powerlaw = abs(y2*WDnorm)+abs(p[0])*x1**(-abs(p[1]))
    #for i in range(N):
        #powerlaw[i] = y2[i]*p[0]+p[1]*x1[i]**p[2]
    Results = calc_Residual(N, isErr, y1, powerlaw, y3)
    return Results

##########################################################################

def errfunc3(p, isErr, x1, y1, x2, y2, y3):
    #isErr,WaveSpec,Spectrum+SpecAddedNoise,WaveModel,SpecModel,SpecErr
    N=len(x1)
    powerlaw = list(x2)
    powerlaw = abs(y2*p[0])+abs(p[1])*x2**(-abs(p[2]))
    
    WaveModelNorm,SpecModelNorm,Cont = WDcontNorm(x2,powerlaw,False)
    powerlawnew = interp(x1,WaveModelNorm,SpecModelNorm)
    
    #for i in range(N):
        #powerlaw[i] = y2[i]*p[0]+p[1]*x1[i]**p[2]
    Results = calc_Residual(N, isErr, y1, powerlawnew, y3)
    return Results

##########################################################################

def errfuncPLfixed(p, isErr, PL, x1, y1, y2, y3):
    #isErr,WaveSpec,Spectrum+SpecAddedNoise, SpecModelNew,SpecErr
    N=len(x1)
    powerlaw = list(x1)
    powerlaw = abs(y2*p[0])+abs(p[1])*x1**PL
    #for i in range(N):
        #powerlaw[i] = y2[i]*p[0]+p[1]*x1[i]**PL
    Results = calc_Residual(N, isErr, y1, powerlaw, y3)
    return Results

##########################################################################

def errfuncFixedModel(p, isErr, x1, y1, y2, y3):
    #isErr,WaveSpec,Spectrum+SpecAddedNoise, SpecModelNew,SpecErr
    N=len(x1)
    fixedmodel = list(x1)
    fixedmodel = y2*abs(p[0])+abs(p[1])*FixedModelSpec
    Results = calc_Residual(N, isErr, y1, fixedmodel, y3)
    return Results

##########################################################################


def calc_Residual(N, isErr, a1, a2, a3):
    """ calculate the residual between a1 and a2 (a3 is the error) """
    if isErr:
        Chi2 = sqrt(sum(((a1-a2)/a3)**2) / (N-1))
    else:
        Chi2 = sqrt(sum((a1-a2)**2) / (N-1))
    #print Chi2
    return Chi2

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

def SpecNorm(aa,bb,IsGraph=False):
    """
    The normalisation of a WD spectrum
    """
    idx_beg_1=find_nearest(aa,3857.)
    idx_end_1=find_nearest(aa,3863.)
    idx_beg_2=find_nearest(aa,3920.)
    idx_end_2=find_nearest(aa,3930.)
    idx_beg_3=find_nearest(aa,4015.)
    idx_end_3=find_nearest(aa,4040.)
    idx_beg_4=find_nearest(aa,4190.)
    idx_end_4=find_nearest(aa,4230.)
    idx_beg_5=find_nearest(aa,4530.)
    idx_end_5=find_nearest(aa,4620.)
    idx_beg_6=find_nearest(aa,5100.)
    idx_end_6=find_nearest(aa,6376.)
    idx_beg_7=find_nearest(aa,6750.)
    idx_end_7=find_nearest(aa,6779.5)
    wd_cont_wave = []
    wd_cont_flux = []
    wave = []
    flux = []    
    wd_cont_wave.extend(aa[idx_beg_1-1:idx_end_1+1])
    wd_cont_flux.extend(bb[idx_beg_1-1:idx_end_1+1])
    wd_cont_wave.extend(aa[idx_beg_2:idx_end_2+1])
    wd_cont_flux.extend(bb[idx_beg_2:idx_end_2+1])
    wd_cont_wave.extend(aa[idx_beg_3:idx_end_3+1])
    wd_cont_flux.extend(bb[idx_beg_3:idx_end_3+1])
    wd_cont_wave.extend(aa[idx_beg_4:idx_end_4+1])
    wd_cont_flux.extend(bb[idx_beg_4:idx_end_4+1])
    wd_cont_wave.extend(aa[idx_beg_5:idx_end_5+1])
    wd_cont_flux.extend(bb[idx_beg_5:idx_end_5+1])
    wd_cont_wave.extend(aa[idx_beg_6:idx_end_6+1])
    wd_cont_flux.extend(bb[idx_beg_6:idx_end_6+1])
    wd_cont_wave.extend(aa[idx_beg_7:idx_end_7+2])
    wd_cont_flux.extend(bb[idx_beg_7:idx_end_7+2])

    wave.extend(aa[idx_beg_1-1:idx_end_7+2])
    flux.extend(bb[idx_beg_1-1:idx_end_7+2])

    #spline = splrep(wd_cont_wave,wd_cont_flux,k=3)
    #continuum = splev(aa,spline)
    #normspec = bb / continuum
    #print continuum
    
    p15 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 15))
    continuum = p15(wave)
    #p20 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 20))
    #continuum20 = p20(aa)
    #p10 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 10))
    #continuum10 = p10(aa)
    
    normspec = flux / continuum
    
    
    if IsGraph:
        plt.figure(1)
        ##plt.plot(wd_cont_wave,wd_cont_flux,'b-',label='continuum')
        #plt.plot(aa,bb,'k-',label='spectrum')
        #plt.plot(aa,continuum15,'b-',label='p15')    
        #plt.plot(aa,continuum20,'g--',label='p20m')    
        #plt.plot(aa,continuum10,'r-',label='p10')  
        plt.plot(wave,flux,'r-',label='Spectrum')    
        plt.plot(wave,continuum,'b-',label='Continuum')    
        
        #plt.plot(wave,normspec,'b-',label='continuum')    
        plt.legend()
        #plt.show() # show the window
        #exit()
    return wave,normspec


##########################################################################

def WDcontNorm(aa,bb,IsGraph=False):
    """
    The normalisation of a WD spectrum
    """
    
    #H[0]=6562.80
    #dH[0]=80.
    #H[1]=4861.33
    #dH[1]=140.
    #H[2]=4340.47
    #dH[2]=120.
    #H[3]=4101.74
    #dH[3]=75.
    #H[4]=3970.07
    #H[5]=3889.05
    #H[6]=3835.40

    H =  [6562.80,4861.33,4340.47,4101.74,3970.07,3889.05,3835.40]
    dH = [80.,140.,120.,75.,50.,40.,30]
    wavenew = []
    contnew = []
    normspecnew = []
   

    ibeg1=find_nearest(aa,3811.)
    iend1=find_nearest(aa,3817.)
    ibeg2=find_nearest(aa,3857.)
    iend2=find_nearest(aa,3863.)
    ibeg3=find_nearest(aa,3920.)
    iend3=find_nearest(aa,3930.)
    ibeg4=find_nearest(aa,4010.)
    iend4=find_nearest(aa,4030.)
    wd_cont_wave = []
    wd_cont_flux = []
    wd_cont_wave.extend(aa[ibeg1:iend1])
    wd_cont_flux.extend(bb[ibeg1:iend1])
    wd_cont_wave.extend(aa[ibeg2:iend2])
    wd_cont_flux.extend(bb[ibeg2:iend2])
    wd_cont_wave.extend(aa[ibeg3:iend3])
    wd_cont_flux.extend(bb[ibeg3:iend3])
    wd_cont_wave.extend(aa[ibeg4:iend4])
    wd_cont_flux.extend(bb[ibeg4:iend4])
    wave = []
    flux = []        
    wave.extend(aa[ibeg1:iend4])
    flux.extend(bb[ibeg1:iend4])

    PolIdx=1
    if wd_cont_wave[0] < 3857.:
        PolIdx = 2
    p3 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, PolIdx))
    continuum = p3(wave)
    normspec = flux / continuum
    wavenew.extend(wave)
    contnew.extend(continuum)
    normspecnew.extend(normspec)

    ibeg1=find_nearest(aa,4023.)
    iend1=find_nearest(aa,4037.)
    ibeg2=find_nearest(aa,4190.)
    iend2=find_nearest(aa,4230.)
    wd_cont_wave = []
    wd_cont_flux = []
    wd_cont_wave.extend(aa[ibeg1:iend1])
    wd_cont_flux.extend(bb[ibeg1:iend1])
    wd_cont_wave.extend(aa[ibeg2:iend2])
    wd_cont_flux.extend(bb[ibeg2:iend2])
    wave = []
    flux = []        
    wave.extend(aa[iend4+1:iend2])
    flux.extend(bb[iend4+1:iend2])
    
    p2 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 2))
    continuum2 = p2(wave)
    normspec2 = flux / continuum2
    wavenew.extend(wave)    
    normspecnew.extend(normspec2)
    contnew.extend(continuum2)

    ibeg1=find_nearest(aa,4190.)
    iend1=find_nearest(aa,4230.)
    ibeg2=find_nearest(aa,4480.)
    iend2=find_nearest(aa,4500.)
    wd_cont_wave = []
    wd_cont_flux = []
    wd_cont_wave.extend(aa[ibeg1:iend1])
    wd_cont_flux.extend(bb[ibeg1:iend1])
    wd_cont_wave.extend(aa[ibeg2:iend2])
    wd_cont_flux.extend(bb[ibeg2:iend2])
    wave = []
    flux = []        
    wave.extend(aa[iend1+1:iend2])
    flux.extend(bb[iend1+1:iend2])
    
    p2 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 2))
    continuum2 = p2(wave)
    normspec2 = flux / continuum2
    contnew.extend(continuum2)
    wavenew.extend(wave)    
    normspecnew.extend(normspec2)

    ibeg1=find_nearest(aa,4640.)
    iend1=find_nearest(aa,4720.)
    ibeg2=find_nearest(aa,5000.)
    iend2=find_nearest(aa,5140.) #
    ibeg3=find_nearest(aa,5200.)
    iend3=find_nearest(aa,5250.)
    wd_cont_wave = []
    wd_cont_flux = []
    wd_cont_wave.extend(aa[ibeg1:iend1])
    wd_cont_flux.extend(bb[ibeg1:iend1])
    wd_cont_wave.extend(aa[ibeg2:iend2])
    wd_cont_flux.extend(bb[ibeg2:iend2])
    wd_cont_wave.extend(aa[ibeg3:iend3])
    wd_cont_flux.extend(bb[ibeg3:iend3])
    wave = []
    flux = []        
    wave.extend(aa[ibeg1:iend3])
    flux.extend(bb[ibeg1:iend3])
    
    p2 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 2))
    continuum2 = p2(wave)
    normspec2 = flux / continuum2
    contnew.extend(continuum2)
    wavenew.extend(wave)    
    normspecnew.extend(normspec2)

    wd_cont_wave = []
    wd_cont_flux = []
    wave = []
    flux = []        
    ibeg1=find_nearest(aa,6250.)
    iend1=find_nearest(aa,6440.)
    if ibeg1 != iend1:
        wd_cont_wave.extend(aa[ibeg1:iend1])
        wd_cont_flux.extend(bb[ibeg1:iend1])
    ibeg2=find_nearest(aa,6680.)
    iend2=find_nearest(aa,6779.5)
    if ibeg2 != iend1 or ibeg2 != iend2:
        wd_cont_wave.extend(aa[ibeg2:iend2])
        wd_cont_flux.extend(bb[ibeg2:iend2])
    wave.extend(aa[ibeg1:iend2])
    flux.extend(bb[ibeg1:iend2])
    
    if len(wd_cont_wave) > 2:
        p2 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 2))
        continuum2 = p2(wave)
        normspec2 = flux / continuum2
        contnew.extend(continuum2)
        wavenew.extend(wave)    
        normspecnew.extend(normspec2)    

    if IsGraph:
        plt.figure(1)
        ##plt.plot(wd_cont_wave,wd_cont_flux,'b-',label='continuum')
        plt.plot(aa,bb,'k-',label='spectrum')
        #plt.plot(aa,continuum15,'b-',label='p15')    
        #plt.plot(aa,continuum20,'g--',label='p20m')    
        #plt.plot(aa,continuum10,'r-',label='p10')  
        
        #plt.plot(wave,flux,'r-',label='Spectrum')    
        plt.plot(wavenew,contnew,'b-',label='Continuum')    
        plt.legend()
        plt.show() # show the window

        plt.figure(3)        
        plt.plot(wavenew,normspecnew,'b-',label='continuum')    
        plt.legend()
        plt.show() # show the window
        #exit()
        
    #return w,n
    return wavenew,normspecnew,contnew
    
    
    
    idx_beg_1=find_nearest(aa,3857.)
    idx_end_1=find_nearest(aa,3863.)
    idx_beg_2=find_nearest(aa,3920.)
    idx_end_2=find_nearest(aa,3930.)
    idx_beg_3=find_nearest(aa,4015.)
    idx_end_3=find_nearest(aa,4040.)
    idx_beg_4=find_nearest(aa,4190.)
    idx_end_4=find_nearest(aa,4230.)
    idx_beg_5=find_nearest(aa,4530.)
    idx_end_5=find_nearest(aa,4620.)
    idx_beg_6=find_nearest(aa,H[0]-220.)
    idx_end_6=find_nearest(aa,H[0]-140.)
    idx_beg_7=find_nearest(aa,H[0]+140.)
    idx_end_7=find_nearest(aa,6779.5)
    wd_cont_wave = []
    wd_cont_flux = []
    wave = []
    flux = []    
    wd_cont_wave.extend(aa[idx_beg_1-1:idx_end_1+1])
    wd_cont_flux.extend(bb[idx_beg_1-1:idx_end_1+1])
    wd_cont_wave.extend(aa[idx_beg_2:idx_end_2+1])
    wd_cont_flux.extend(bb[idx_beg_2:idx_end_2+1])
    wd_cont_wave.extend(aa[idx_beg_3:idx_end_3+1])
    wd_cont_flux.extend(bb[idx_beg_3:idx_end_3+1])
    wd_cont_wave.extend(aa[idx_beg_4:idx_end_4+1])
    wd_cont_flux.extend(bb[idx_beg_4:idx_end_4+1])
    wd_cont_wave.extend(aa[idx_beg_5:idx_end_5+1])
    wd_cont_flux.extend(bb[idx_beg_5:idx_end_5+1])
    wd_cont_wave.extend(aa[idx_beg_6:idx_end_6+1])
    wd_cont_flux.extend(bb[idx_beg_6:idx_end_6+1])
    wd_cont_wave.extend(aa[idx_beg_7:idx_end_7+2])
    wd_cont_flux.extend(bb[idx_beg_7:idx_end_7+2])

    wave.extend(aa[idx_beg_1-1:idx_end_7+2])
    flux.extend(bb[idx_beg_1-1:idx_end_7+2])

    #spline = splrep(wd_cont_wave,wd_cont_flux,k=3)
    #continuum = splev(aa,spline)
    #normspec = bb / continuum
    #print continuum
    
    p15 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 15))
    continuum = p15(wave)
    #p20 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 20))
    #continuum20 = p20(aa)
    #p10 = poly1d(polyfit(wd_cont_wave,wd_cont_flux, 10))
    #continuum10 = p10(aa)
    
    normspec = flux / continuum
    
    
    
    
    if IsGraph:
        plt.figure(1)
        ##plt.plot(wd_cont_wave,wd_cont_flux,'b-',label='continuum')
        #plt.plot(aa,bb,'k-',label='spectrum')
        #plt.plot(aa,continuum15,'b-',label='p15')    
        #plt.plot(aa,continuum20,'g--',label='p20m')    
        #plt.plot(aa,continuum10,'r-',label='p10')  
        plt.plot(wave,flux,'r-',label='Spectrum')    
        plt.plot(wave,continuum,'b-',label='Continuum')    
        
        #plt.plot(wave,normspec,'b-',label='continuum')    
        plt.legend()
        #plt.show() # show the window
        #exit()
    return wavenew,normspec

##########################################################################

def SpecLineCut(aa,bb,cc,IsErr=True):
    """
    Selection of spectral intervals for the fitting process
    """
    
    H =  [6562.80,4861.33,4340.47,4101.74,3970.07,3889.05,3835.40]
    dH = [80.,140.,120.,75.,50.,40.,30]
    
    wave = []
    flux = [] 
    fluxerr = []
    ibeg1=find_nearest(aa,3814.)
    iend1=find_nearest(aa,H[3]+dH[3])
    if ibeg1 != iend1:
        wave.extend(aa[ibeg1:iend1])
        flux.extend(bb[ibeg1:iend1])   
        fluxerr.extend(cc[ibeg1:iend1])   
    ibeg2=find_nearest(aa,H[2]-dH[2])
    iend2=find_nearest(aa,H[2]+dH[2])
    if ibeg2 != iend1 or ibeg2 != iend2:
        wave.extend(aa[ibeg2:iend2])
        flux.extend(bb[ibeg2:iend2])   
        fluxerr.extend(cc[ibeg2:iend2])   
    ibeg3=find_nearest(aa,H[1]-dH[1])
    iend3=find_nearest(aa,H[1]+dH[1])
    if ibeg3 != iend2 or ibeg3 != iend3:
        wave.extend(aa[ibeg3:iend3])
        flux.extend(bb[ibeg3:iend3])   
        fluxerr.extend(cc[ibeg3:iend3])   
    ibeg4=find_nearest(aa,H[0]-dH[0])
    iend4=find_nearest(aa,H[0]+dH[0])
    if ibeg4 != iend3 or ibeg4 != iend4:
        wave.extend(aa[ibeg4:iend4])
        flux.extend(bb[ibeg4:iend4])   
        fluxerr.extend(cc[ibeg4:iend4])   
        
    return wave,flux,fluxerr


##################################################################

def SpecEmissionLineCut(dV,aa,bb,cc,IsErr=True,IsHe=False):
    """
    Selection of spectral intervals for the fitting process
    """    
    c = 299792.4
    H =  [6562.80, 4861.33, 4340.47, 4101.74, 3970.07]
    #He = [4471.48,4921.93,5017.,5875.7,6678.]    
    He = [4471.48,4686.,5017.,5169.,5875.7,6678.]    
    if IsHe:
        H.extend(He)
        H.sort(reverse=True)
    
    wave = []
    flux = [] 
    fluxerr = []
    
    ibeg=0
    iend=len(aa)-1
    i=len(H)-1   
    
    if aa[0] <= H[i]+dV/c*H[i]:
        #ibeg=find_nearest(aa,H[i]+dV/c*H[i])
        ibeg=find_nearest(aa,H[i]+dV/c*H[i])+1
    else:
        i-=1
        while aa[0] >= H[i]-dV/c*H[i]:
            i-=1
        if aa[0] <= H[i+1]+dV/c*H[i+1]:
            #ibeg=find_nearest(aa,H[i+1]+dV/c*H[i+1])
            ibeg=find_nearest(aa,H[i+1]+dV/c*H[i+1])+1

    while i>0:
        i-=1
        iend=find_nearest(aa,H[i]-dV/c*H[i])
        if iend > ibeg:
            wave.extend(aa[ibeg:iend])
            flux.extend(bb[ibeg:iend])   
            fluxerr.extend(cc[ibeg:iend])
            #ibeg=find_nearest(aa,H[i]+dV/c*H[i])
            ibeg=find_nearest(aa,H[i]+dV/c*H[i])+1
        else: break
    if (aa[-1]>aa[iend] and aa[-1] > H[0]+dV/c*H[0]):
        iend=find_nearest(aa,6800)
        wave.extend(aa[ibeg:iend])
        flux.extend(bb[ibeg:iend])   
        fluxerr.extend(cc[ibeg:iend])

    return wave,flux,fluxerr


##################################################################

def FittingCycle(irun):
    NormFactor = []
    PowerFactor = []
    Chi2 = []
    Index = []
    ModelFile = []
    PowerModel = []
    Chi2Best = 1e9
    i = 0
    imodels = 0
    for line in lines:        
        imodels += 1
        #print('\033[K   Current model: ', imodels, '\r',)         # python2
        #print("{}\rCurrent model: ", imodels,end="")
        #print('\x1b[2K\r',)
        print('{}\rCurrent model: '.format(imodels),end="    ")   # python3
        #{0}\r        
        #print("\033[K")
        #sys.stdout.write("\033[K")        
        #print("Current model: ", imodels,end="")
        #sys.stdout.write("Current model: ", imodels)
        #print '\x1b[2K\r',
        sys.stdout.flush()
        
        sublines = line.split()
        #print(len(sublines))
        FileName = sublines[0]
        if len(sublines)>1:
            Teff = float(sublines[1])
#            print(Teff,'\r',)
        if len(sublines)>2:
            Logg = float(sublines[2])/100.
#            print(Logg,'\r',)
        
        #FileName = line.rstrip('\n')
        data0 = loadtxt(FileName, usecols=[0,1], unpack=True,skiprows=0)
        ModelFile.append(FileName)
        WaveModel = data0[0,:]
        SpecModel = data0[1,:]

        if FitModel == 3:
            #if isErr:
                #WaveSpec,Spectrum,SpecErr = SpecLineCut(array(WaveSpec),array(Spectrum),array(SpecErr),isErr)
            #else:
                #WaveSpec,Spectrum,SpecErr = SpecLineCut(array(WaveSpec),array(Spectrum),Spectrum,isErr)   
            WaveModelNorm,SpecModelNorm,Cont = WDcontNorm(WaveModel,SpecModel,False)
            SpecModelNew = interp(WaveSpec,WaveModelNorm,SpecModelNorm)
#### Temp    
            #plt.figure(4)        
            #plt.plot(WaveSpec,SpecModelNew,'b-',label='Model')    
            #plt.plot(WaveSpec,Spectrum,'r-',label='Star')    
            #plt.legend()
            #plt.show() # show the window
#### Temp
            
        elif FitModel == 9:
            WaveModelNorm,SpecModelNorm = SpecNorm(WaveModel,SpecModel)
            SpecModelNew = interp(WaveSpec,WaveModelNorm,SpecModelNorm)
        else:
            SpecModelNew = interp(WaveSpec,WaveModel,SpecModel)      
                   
        InitialGuessNorm = sum(Spectrum) / sum(SpecModelNew)
        #print(InitialGuessNorm)
        
        if FitModel <= 0.00:
            PowerIndex = InitialGuessPowerIndex = FitModel
            PowerNorm = InitialGuessPowerNorm = 0.1 * sum(Spectrum) / sum(WaveSpec**(FitModel))
        elif FitModel == 3 or FitModel == 4:
            PowerIndex = InitialGuessPowerIndex = -1.4
            PowerNorm = InitialGuessPowerNorm = 0.1         
        elif FitModel == 5:
            PowerIndex = InitialGuessPowerIndex = 0
            PowerNorm = InitialGuessPowerNorm = 0.1 * sum(Spectrum) / sum(FixedModelSpec)
        elif FitModel == 21:
            if len(sublines)>2:
                Rwd = (0.4015-0.03606*Logg)**2
                InitialGuessNorm = (Rwd/4.434e7/Dist)**2
            else:
                print('Log g is not given. The default value of 8.35 is used')
                InitialGuessNorm = 4.2654351716654616e-24
            PowerIndex = InitialGuessPowerIndex = -1.4
            PowerNorm = InitialGuessPowerNorm = 0.1 * sum(Spectrum) / sum(WaveSpec**(PowerIndex))
        else:
            PowerIndex = InitialGuessPowerIndex = -1.4
            PowerNorm = InitialGuessPowerNorm = 0.1 * sum(Spectrum) / sum(WaveSpec**(PowerIndex))
        
        if isErr:
            if irun>1:                
                SpecAddedNoise = random.normal(0., 1., len(Spectrum))*std(SpecErr)
                #random.shuffle(mylist)
            else:
                SpecAddedNoise = 0.
            #Norm_Factor, Chi2Min, Qiter, funcalls = opt.brent(errfunc, brack=(InitialGuess/10.,InitialGuess,InitialGuess*10.), args=(isErr,Spectrum+SpecAddedNoise, SpecModelNew,SpecErr),full_output=1)
            if FitModel == 1 or (FitModel == 3 or FitModel == 9):
                Norm_Factor, Chi2Min, Qiter, funcalls = opt.brent(errfunc1, brack=(InitialGuessNorm/10.,InitialGuessNorm,InitialGuessNorm*10.), args=(isErr,array(Spectrum)+SpecAddedNoise,SpecModelNew,SpecErr),full_output=1)
                PowerNorm = 0.0
                PowerIndex = 0.0
            elif FitModel == 2:
                p2 = [InitialGuessNorm, InitialGuessPowerNorm, InitialGuessPowerIndex] # Initial guess for the parameters
                Result = opt.minimize(errfunc2,p2,method='Nelder-Mead',args=(isErr,WaveSpec,Spectrum+SpecAddedNoise, SpecModelNew,SpecErr),
                                  options={'maxiter': 10000},tol=1e-4) #,bounds=bnds)
                #print Result
                Norm_Factor = Result.x[0]
                #print(PowerNorm)
                PowerNorm = Result.x[1]
                PowerIndex = Result.x[2]
                Chi2Min = Result.fun
            elif FitModel == 21:
                p2 = [InitialGuessPowerNorm, InitialGuessPowerIndex] # Initial guess for the parameters
                Result = opt.minimize(errfunc21,p2,method='Nelder-Mead',args=(isErr,WaveSpec,Spectrum+SpecAddedNoise,SpecModelNew,SpecErr,InitialGuessNorm),
                                  options={'maxiter': 10000},tol=1e-4) #,bounds=bnds)
                #print Result
                Norm_Factor = InitialGuessNorm
                PowerNorm = Result.x[0]
                PowerIndex = Result.x[1]
                Chi2Min = Result.fun
            elif FitModel == 4:
                p2 = [InitialGuessNorm, InitialGuessPowerNorm, InitialGuessPowerIndex] # Initial guess for the parameters
                Result = opt.minimize(errfunc3,p2,method='Nelder-Mead',args=(isErr,WaveSpec,array(Spectrum)+SpecAddedNoise, WaveModel,SpecModel,SpecErr),
                                  options={'maxiter': 10000},tol=1e-4) #,bounds=bnds)
                #print Result
                Norm_Factor = Result.x[0]
                PowerNorm = Result.x[1]
                PowerIndex = Result.x[2]
                Chi2Min = Result.fun
            #if FitModel == 9:
                #SpecNorm(WaveSpec,SpecModelNew)
            elif FitModel == 5:
                p0 = [InitialGuessNorm, InitialGuessPowerNorm] # Initial guess for the parameters
                Result = opt.minimize(errfuncFixedModel,p0,method='Nelder-Mead',args=(isErr,WaveSpec,Spectrum+SpecAddedNoise, SpecModelNew,SpecErr),options={'maxiter': 10000}) #,bounds=bnds)
                #print Result
                Norm_Factor = Result.x[0]
                PowerNorm = Result.x[1]
                PowerIndex = -999
                Chi2Min = Result.fun
            elif FitModel <= 0.00:
                p0 = [InitialGuessNorm, InitialGuessPowerNorm] # Initial guess for the parameters
                Result = opt.minimize(errfuncPLfixed,p0,method='Nelder-Mead',args=(isErr,FitModel,WaveSpec,Spectrum+SpecAddedNoise, SpecModelNew,SpecErr),options={'maxiter': 10000}) #,bounds=bnds)
                #print Result
                Norm_Factor = Result.x[0]
                PowerNorm = Result.x[1]
                PowerIndex = FitModel
                Chi2Min = Result.fun
            else:
                print("Wrong model. Stop executing!")
                exit()            
        else:
            if FitModel == 1 or (FitModel == 3 or FitModel == 9):
                Norm_Factor, Chi2Min, Qiter, funcalls = opt.brent(errfunc1, brack=(InitialGuessNorm/10.,InitialGuessNorm,InitialGuessNorm*10.), args=(isErr,Spectrum, SpecModelNew,Spectrum),full_output=1)
                PowerNorm = 0.0
                PowerIndex = 0.0
            elif FitModel == 2:
                p2 = [InitialGuessNorm, InitialGuessPowerNorm, InitialGuessPowerIndex] # Initial guess for the parameters
                Result = opt.minimize(errfunc2,p2,method='Nelder-Mead',args=(isErr,WaveSpec,Spectrum, SpecModelNew,Spectrum),options={'maxiter': 10000})
                #print Result
                Norm_Factor = Result.x[0]                
                PowerNorm = Result.x[1]                
                PowerIndex = Result.x[2]
                Chi2Min = Result.fun
            elif FitModel == 21:
                p2 = [InitialGuessPowerNorm, InitialGuessPowerIndex] # Initial guess for the parameters
                Result = opt.minimize(errfunc21,p2,method='Nelder-Mead',args=(isErr,WaveSpec,Spectrum,SpecModelNew,Spectrum,InitialGuessNorm),
                                  options={'maxiter': 10000},tol=1e-4) #,bounds=bnds)
                #print(Result)
                Norm_Factor = InitialGuessNorm
                PowerNorm = Result.x[0]
                PowerIndex = Result.x[1]
                Chi2Min = Result.fun                
            elif FitModel == 4:
                p2 = [InitialGuessNorm, InitialGuessPowerNorm, InitialGuessPowerIndex] # Initial guess for the parameters
                Result = opt.minimize(errfunc3,p2,method='Nelder-Mead',args=(isErr,WaveSpec,Spectrum, WaveModel,SpecModel,Spectrum),
                                  options={'maxiter': 10000},tol=1e-4) #,bounds=bnds)
                #print Result
                Norm_Factor = Result.x[0]
                PowerNorm = Result.x[1]
                PowerIndex = Result.x[2]
                #print InitialGuessPowerIndex,PowerIndex
                Chi2Min = Result.fun                
            elif FitModel == 5:
                p0 = [InitialGuessNorm, InitialGuessPowerNorm] # Initial guess for the parameters
                Result = opt.minimize(errfuncFixedModel,p0,method='Nelder-Mead',args=(isErr,WaveSpec,Spectrum, SpecModelNew,Spectrum),options={'maxiter': 10000}) #,bounds=bnds)
                #print Result
                Norm_Factor = Result.x[0]
                PowerNorm = Result.x[1]
                PowerIndex = -999
                Chi2Min = Result.fun
            elif FitModel <= 0.00:
                p0 = [InitialGuessNorm, InitialGuessPowerNorm] # Initial guess for the parameters
                Result = opt.minimize(errfuncPLfixed,p0,method='Nelder-Mead',args=(isErr,FitModel,WaveSpec,Spectrum, SpecModelNew,Spectrum),options={'maxiter': 10000})
                #print Result
                Norm_Factor = Result.x[0]
                PowerNorm = Result.x[1]
                PowerIndex = FitModel
                Chi2Min = Result.fun
            else:
                print("Wrong model. Stop executing!")
                exit()            
   
        NormFactor.append(abs(Norm_Factor))
        PowerFactor.append(abs(PowerNorm))
        Index.append(PowerIndex)
        Chi2.append(Chi2Min)    
        if Chi2Min < Chi2Best:
            Chi2Best = Chi2Min
            BestModel = i
        i+=1    
    ModelFileBest.append(ModelFile[BestModel])
    Chi2best.append(Chi2[BestModel])
    NormFactorBest.append(NormFactor[BestModel])
    PowerIndexBest.append(Index[BestModel])
    PowerFactorBest.append(PowerFactor[BestModel])
    return Chi2Best,BestModel,ModelFile,NormFactor,PowerFactor,Index,Chi2


##########################################################################


#https://python4mpia.github.io/fitting_data/least-squares-fitting.html
#http://www.astroml.org/book_figures/chapter4/fig_chi2_eval.html
#https://scipy.github.io/old-wiki/pages/Cookbook/Least_Squares_Circle.html

#########################################################################
#########################################################################
#########################################################################
#########################################################################
def SpecTesting(x,y):
    N=len(x)
    y2 = list(x)
    PL = input("Enter the Power-Law index: ")
    PL = float(PL)
    Norm = input("Enter the contribution of the Power-Law spectrum (%): ")
    Norm = float(Norm)
    if PL>0.0:
        PL = (-1) * PL
    y2 = x**PL
    sumy1=sum(y)
    sumy2=sum(y2)
    y2=y2/sumy2*sumy1
    y2=y2*Norm/100.+y
    #for i in range(N):
        #y[i] = y[i]+norm*x1[i]**PL
    return y2


#########################################################################

if __name__ == "__main__":

    #multiprocessing.cpu_count()
    printHeader()
    if (len(sys.argv) < 4):
        printUsage()
    
    if (len(sys.argv) > 1):
        try:
            FitModel=float(sys.argv[1])
        except:
            if sys.argv[1] == 'NL':
                FitModel = 91 # Normalization of the observed spectrum and the cutting Balmer lines
            elif sys.argv[1] == 'NC':
                FitModel = 92 # Normalization of the observed spectrum (another algorithm)
            elif sys.argv[1] == '--cut':
                HeLines = True
                if (len(sys.argv) > 2):
                    SpecFile = sys.argv[2]
                else:
                    print()
                    sys.stdout.write('Enter the filename of the observed spectrum: ')
                    sys.stdout.flush()
                    output_file_path = stdin.readline()
                    SpecFile = output_file_path.rstrip('\n')    
                if not os.path.exists(SpecFile): 
                    print("The file",SpecFile, "doesn't exist. Stop executing!")
                    exit()
                else:
                    data1 = loadtxt(SpecFile, unpack=True, skiprows=0)
                    WaveSpec = data1[0,:]
                    Spectrum = data1[1,:]
                    try:
                       SpecErr = data1[2,:]
                       isErr = 1
                       #print ("The Object spectrum has the 3rd column. Will assume that it includes error data.")
                    except IndexError:
                       isErr = 0 
                       SpecErr = data1[1,:]
                if (len(sys.argv) > 3):
                    VelCut = float(sys.argv[3])
                else:
                    print()
                    VelCut = input("Enter the range of velocities to be cut from the spectral lines (km/s): ")
                    VelCut = float(VelCut)
                if (len(sys.argv) > 4):
                    Helium=sys.argv[4]
                    if Helium.lower() == 'no':
                        HeLines = False                
                        print(Helium.lower())
                plt.figure(4)        
                plt.plot(WaveSpec,Spectrum,'b-',label='Original spectrum')                   
                WaveSpec,Spectrum,SpecErr = SpecEmissionLineCut(VelCut,WaveSpec,Spectrum,SpecErr,IsHe=HeLines)
                if isErr:
                    WriteData3(len(WaveSpec),WaveSpec,Spectrum,SpecErr,SpecFile+".cut")
                else: WriteData2(len(WaveSpec),WaveSpec,Spectrum,SpecFile+".cut")
                print()
                print ("The file",SpecFile+".cut","is created.")                      
                print()
                plt.plot(WaveSpec,Spectrum,'r-',label='Cut spectrum')    
                plt.legend()
                plt.show() # show the window                
                exit()            
            else:
                printUsage()
                exit()
    else:
        #printUsage()
        FitModel = input("Enter the ModelNumber parameter: ")
        FitModel = float(FitModel)
    
    if FitModel == 1:
        print("FitModel=",FitModel," --> WD only")
    elif FitModel == 2:
        print("FitModel=",FitModel," --> WD + Power-Law")
    elif FitModel == 21:
        print("FitModel=",FitModel," --> WD + Power-Law, Fixed WD normalisation")
        Dist = input("Enter the distance (pc): ")
        Dist = float(Dist)
    elif FitModel == 3:
        print("FitModel=",FitModel," --> Normalized, lines only")
    elif FitModel == 4:
        print("FitModel=",FitModel," --> WD + Power-Law, Normalized, only Lines")
    elif FitModel == 5:
        FixedModel = "FixedModel.dat"
        print("FitModel=",FitModel," --> WD + a Fixed Model")
        #print(os.getcwd())
        if not os.path.exists(FixedModel):
            print("The file",FixedModel,"doesn't exist. Stop executing!")
            exit()        
    elif FitModel == 9:
        print("FitModel=",FitModel," --> Normalized WD only (Experimental!!!)")
    elif FitModel == 91:
        print("Normalization of the observed spectrum and the cutting the Balmer lines")
    elif FitModel == 92:
        print("Normalization of the observed spectrum (another algorithm)")
    elif FitModel <= 0.0:
        print("FitModel is WD + Power-Law with the index ",FitModel)
    else: 
        print("Wrong model. Stop executing!")
        exit()            
        
    #print "FitModel=",FitModel
    
    if (len(sys.argv) > 2):
        SpecFile = sys.argv[2]
    else:
        print()
        sys.stdout.write('Enter the filename of the observed spectrum: ')
        sys.stdout.flush()
        output_file_path = stdin.readline()
        SpecFile = output_file_path.rstrip('\n')    
    if not os.path.exists(SpecFile): 
        print("The file",SpecFile, "doesn't exist. Stop executing!")
        exit()
    else:
        data1 = loadtxt(SpecFile, unpack=True, skiprows=0)
        WaveSpec = data1[0,:]
        Spectrum = data1[1,:]
        if FitModel == 5:
            data2 = loadtxt(FixedModel, unpack=True, skiprows=0)
            if data2[0,-1] < WaveSpec[-1] or data2[0,0] > WaveSpec[0]:
                print("The FixedModel spectrum is shorther than the observed one!")
                print("Wavelengths of the observed spectrum are in the range",WaveSpec[0],"-",WaveSpec[-1])
                print("while The FixedModel spectrum is from", data2[0,0],"to", data2[0,-1])
                print("Stop executing!")
                exit()                
            FixedModelSpec = interp(WaveSpec,data2[0,:],data2[1,:])
            #WaveSpec = data1[0,:]
            #Spectrum = data1[1,:]
        if FitModel == 9:
            if WaveSpec[0] < 3857. or WaveSpec[-1] > 6779.5:
                print ("Wavelengths of the observed spectrum must be in the range 3857.0-6779.5,")
                print ("but it is from", WaveSpec[0],"to", WaveSpec[-1])
                print ("Stop executing!")
                exit()
        try:
            SpecErr = data1[2,:]
            isErr = 1
            print ("The Object spectrum has the 3rd column. Will assume that it includes error data.")
        except IndexError:
            isErr = 0         
    
        if FitModel == 3 or FitModel == 4 or FitModel == 9:
            if isErr:
                WaveSpec,Spectrum,SpecErr = SpecLineCut(WaveSpec,Spectrum,SpecErr,isErr)
            else:
                WaveSpec,Spectrum,SpecErr = SpecLineCut(WaveSpec,Spectrum,Spectrum,isErr)   
    
        #isErr=False
        if FitModel == 91:        
            if isErr:
                WaveNorm,SpecNorm,Cont = WDcontNorm(WaveSpec,Spectrum,True)
                SpecErr = interp(WaveNorm,WaveSpec,SpecErr)
                SpecErr = SpecErr / Cont
                WaveSpec,Spectrum,SpecErr = SpecLineCut(array(WaveNorm),array(SpecNorm),array(SpecErr),isErr)            
                WriteData3(len(WaveSpec),WaveSpec,Spectrum,SpecErr,SpecFile+'.norm')
            else:
                WaveSpec,Spectrum,Cont = WDcontNorm(WaveSpec,Spectrum,True)
                WaveSpec,Spectrum,SpecErr = SpecLineCut(array(WaveSpec),array(Spectrum),Spectrum,isErr)   
                WriteData2(len(WaveSpec),WaveSpec,Spectrum,SpecFile+'.norm')
            print()
            print ("The file",SpecFile+".norm","is created.")
            exit()
    
        if FitModel == 92:        
            if WaveSpec[0] < 3857. or WaveSpec[-1] > 6779.5:
                print ("Wavelengths of the observed spectrum must be in the range 3857.0-6779.5,")
                print ("but it is from", WaveSpec[0],"to", WaveSpec[-1])
                print ("Stop executing!")
                exit()        
            if isErr:
                WaveNorm,SpecNorm,Cont = WDcontNorm(WaveSpec,Spectrum,True)
                SpecErr = interp(WaveNorm,WaveSpec,SpecErr)
                SpecErr = SpecErr / Cont
                WaveSpec,Spectrum,SpecErr = SpecLineCut(array(WaveNorm),array(SpecNorm),array(SpecErr),isErr)            
                WriteData3(len(WaveSpec),WaveSpec,Spectrum,SpecErr,SpecFile+'.norm')
            else:
                WaveSpec,Spectrum = SpecNorm(WaveSpec,Spectrum,True)
                WriteData2(len(WaveSpec),WaveSpec,Spectrum,SpecFile+'.norm')
    
            #WaveModelNorm,SpecModelNorm = SpecNorm(WaveModel,SpecModel)
            #SpecModelNew = interp(WaveSpec,WaveModelNorm,SpecModelNorm)
    #### Temp    
            plt.figure(4)        
            #plt.plot(WaveSpec,SpecModelNew,'b-',label='Model')    
            plt.plot(WaveSpec,Spectrum,'r-',label='Star')    
            plt.legend()
            plt.show() # show the window
    
            print()
            print("The file",SpecFile+".norm","is created.")
            exit()
    
    #Spectrum=SpecTesting(WaveSpec,Spectrum)
    
    if (len(sys.argv) > 3):
        FileList = sys.argv[3]
    else:   
        #print "Enter the filename of the list of model spectra: "
        print()    
        sys.stdout.write('Enter the filename of the list of model spectra: ')
        sys.stdout.flush()
        output_file_path = stdin.readline()
        FileList = output_file_path.rstrip('\n')    
    if not os.path.exists(FileList): 
        print ("The file",FileList, "doesn't exist. Stop executing!")
        exit()
    else:
        infile = open(FileList, "r")
        lines = infile.readlines()
        infile.close()
        #print "Number of models: ",len(lines)
    
    if (len(sys.argv) > 4):
        ResultDir = sys.argv[4]
    else:
        ResultDir = os.getcwd()
    
    if (len(sys.argv) > 5) and isErr:
        Runs = int(sys.argv[5])
        print()   
        print ("The comand line includes the 5th parameter which switches the program to the multi-run mode.")
        print ("The program will be executed",Runs,"times with the gaussian-noise-added input spectrum.")
        print()   
        print ("The best models will be written to the Dir ", ResultDir)
        #print(os.getcwd())
        #CPUcount= cpu_count()
        #print "There are", CPUcount, "CPUs"
    else:
        print()
        if ResultDir == os.getcwd():
            print ("Results will be written to the current dir ", ResultDir)
        else:
            print ("Results will be written to the Dir ", ResultDir)
        #print(os.getcwd())
        Runs=1
    
    #print "is Error?", isErr
    
    ModelFileBest = []
    NormFactorBest = []
    PowerFactorBest = []
    PowerIndexBest = []
    Chi2best = []
    
    print()   
    print("Number of models: ",len(lines))
    irun=1
    
    ################################ Fitting Cycle ###############################
    while irun<=Runs:        
        if irun == 2:
            print("****** Monte Carlo will NOT work for all models!!! ******")
            FileName = ModelFile[BestModel]
            data0 = loadtxt(FileName, usecols=[0,1], unpack=True,skiprows=0)
            WaveModel = data0[0,:]
            SpecModel = data0[1,:] * NormFactor[BestModel]
            SpecModel = SpecModel + PowerFactor[BestModel]*WaveModel**Index[BestModel]            
            SpecModelNew = interp(WaveSpec,WaveModel,SpecModel)
            Spectrum = SpecModelNew
        if Runs > 1:
            print ("Run=",irun,"                         ")
        
        Chi2Best,BestModel,ModelFile,NormFactor,PowerFactor,Index,Chi2 = FittingCycle(irun)
                
        #ModelFileBest.append(ModelFile[BestModel])
        #Chi2best.append(Chi2[BestModel])
        #NormFactorBest.append(NormFactor[BestModel])
        #PowerIndexBest.append(Index[BestModel])
        #PowerFactorBest.append(PowerFactor[BestModel])
        irun+=1        
    
    if Runs==1:
        print()
        print("Chi2Best=", Chi2Best)
        print("Results:",ModelFile[BestModel],NormFactor[BestModel],PowerFactor[BestModel],Index[BestModel])
        ResultFile=ResultDir+'/Res.'+SpecFile
        ResultBestModelFile=ResultDir+'/BestModel.'+SpecFile
        ResultBestADFile=ResultDir+'/AccDisk.'+SpecFile
        WriteData5(len(ModelFile),ModelFile,Chi2,NormFactor,PowerFactor,Index,ResultFile)
        
        FileName = ModelFile[BestModel]
        data0 = loadtxt(FileName, usecols=[0,1], unpack=True,skiprows=0)
        WaveModel = data0[0,:]
        SpecModel = data0[1,:] * NormFactor[BestModel]

        #print(os.getcwd())
        if FitModel != 5:
            SpecModel = SpecModel + PowerFactor[BestModel]*WaveModel**Index[BestModel]
        if FitModel == 3 or FitModel == 4:
            WaveModelNorm,SpecModelNorm,Cont = WDcontNorm(WaveModel,SpecModel,False)
            SpecModelNew = interp(WaveSpec,WaveModelNorm,SpecModelNorm)
            #WaveSpec = list(WaveSpec)
        elif FitModel == 9:
            WaveModelNorm,SpecModelNorm = SpecNorm(WaveModel,SpecModel,IsGraph=True)
            SpecModelNew = interp(WaveSpec,WaveModelNorm,SpecModelNorm)
        else:
            SpecModelNew = interp(WaveSpec,WaveModel,SpecModel)      
        
        if FitModel == 5:            
            SpecModelNew = SpecModelNew + PowerFactor[BestModel]*FixedModelSpec        
            WriteDataModel(len(WaveSpec),WaveSpec,SpecModelNew, Spectrum,
                interp(WaveSpec,WaveModel,data0[1,:] * NormFactor[BestModel]),
                PowerFactor[BestModel]*FixedModelSpec, Spectrum-SpecModelNew,ResultBestModelFile)
        else:
            WriteDataModel(len(WaveSpec),WaveSpec,SpecModelNew, Spectrum,
                interp(WaveSpec,WaveModel,data0[1,:] * NormFactor[BestModel]),
                PowerFactor[BestModel]*power(WaveSpec,Index[BestModel]), Spectrum-SpecModelNew,ResultBestModelFile)
            accretiondisk = Spectrum-interp(WaveSpec,WaveModel,data0[1,:] * NormFactor[BestModel])
            WriteData2(len(WaveSpec),WaveSpec,accretiondisk,ResultBestADFile)
        
    
    #######################################   Plotting code   #######################################
        fig = plt.figure(2,figsize=(12, 8))
        gs = gridspec.GridSpec(2, 1,
                           height_ratios=[4, 1]
                           )
        gs.update(left=0.08, right=0.95, hspace=0.0)
        figspec = plt.subplot(gs[0])
        figspec.set_title("$\chi^2$="+'%.2f' % Chi2Best+"   WD: "+ModelFile[BestModel],{'color': 'black', 'fontsize': 14})
        figres  = plt.subplot(gs[1])
        figspec.plot(WaveSpec,Spectrum, color="red", lw=1.5, label="Object Spectrum")
        figspec.plot(WaveSpec,SpecModelNew, color="blue", lw=1.5, label="Best Model Spectrum")
        if FitModel == 3:
            WaveModelNorm,SpecModelNorm,Cont = WDcontNorm(WaveModel,SpecModel,False)
            figspec.plot(WaveSpec,interp(WaveSpec,WaveModelNorm,array(SpecModelNorm) * NormFactor[BestModel]), color="gray",lw=1.0, linestyle='dashed',label="WD Spectrum")    
            #WDNorm = interp(WaveModelNorm,WaveModel,data0[1,:] * NormFactor[BestModel])/Cont
            #figspec.plot(WaveModelNorm,WDNorm,color="gray",lw=1.0, linestyle='dashed',label="WD Spectrum")
            WriteDataModel(len(WaveSpec),WaveSpec,SpecModelNew,Spectrum,
                           interp(WaveSpec,WaveModelNorm,array(SpecModelNorm) * NormFactor[BestModel]),
                               PowerFactor[BestModel]*power(WaveSpec,Index[BestModel]), Spectrum-SpecModelNew,ResultBestModelFile)
        elif FitModel == 4:       
            WDNorm = interp(WaveModelNorm,WaveModel,data0[1,:] * NormFactor[BestModel])/Cont
            #interp(WaveSpec,WaveModelNorm,array(SpecModelNorm) * NormFactor[BestModel]/Cont
            #WaveModelNorm,SpecModelNorm,Cont = WDcontNorm(WaveModel,data0[1,:],False)
            
            #SpecModel = data0[1,:] * NormFactor[BestModel] + PowerFactor[BestModel]*WaveModel**Index[BestModel]
            
            #tmp = NormFactor[BestModel]
            #q = SpecModelNorm * tmp
            #SpecModelNew = interp(WaveSpec,WaveModelNorm,array(SpecModelNorm) * NormFactor[BestModel])
            WDNorm = interp(WaveSpec,WaveModelNorm,WDNorm)
            figspec.plot(WaveSpec,WDNorm,color="gray",lw=1.0, linestyle='dashed',label="WD Spectrum")
            PowerLawNorm = interp(WaveModelNorm,WaveModel,PowerFactor[BestModel]*power(WaveModel,Index[BestModel]))/Cont
            PowerLawNorm = interp(WaveSpec,WaveModelNorm,PowerLawNorm)
            figspec.plot(WaveSpec, PowerLawNorm, color="green",lw=1.0, linestyle='dashed',label="Power-law, $\Gamma$="+'%.2e' % Index[BestModel])           
            WriteDataModel(len(WaveSpec),WaveSpec,SpecModelNew,Spectrum,
                   WDNorm,PowerLawNorm, Spectrum-SpecModelNew,ResultBestModelFile)
    
            
        elif FitModel == 9:
            WaveModelNorm,SpecModelNorm = SpecNorm(WaveModel,SpecModel,IsGraph=True)
            #SpecModelNew = interp(WaveSpec,WaveModelNorm,SpecModelNorm)        
            figspec.plot(WaveSpec,interp(WaveSpec,WaveModelNorm,SpecModelNorm * NormFactor[BestModel]), color="gray",lw=1.0, linestyle='dashed',label="WD Spectrum")
        elif FitModel == 5:
            figspec.plot(WaveSpec,interp(WaveSpec,WaveModel,data0[1,:] * NormFactor[BestModel]), color="gray",lw=1.0, linestyle='dashed',label="WD Spectrum")
            figspec.plot(WaveSpec, PowerFactor[BestModel]*FixedModelSpec, color="green",lw=1.0, linestyle='dashed',label="Fixed Model")   
        else:
            figspec.plot(WaveSpec,interp(WaveSpec,WaveModel,data0[1,:] * NormFactor[BestModel]), color="gray",lw=1.0, linestyle='dashed',label="WD Spectrum")
            figspec.plot(WaveSpec, PowerFactor[BestModel]*power(WaveSpec,Index[BestModel]), color="green",lw=1.0, linestyle='dashed',label="Power-law, $\Gamma$="+'%.3f' % Index[BestModel])   
        figspec.legend()
        Ymin = min(Spectrum[len(SpecModelNew)-1],SpecModelNew[len(SpecModelNew)-1],PowerFactor[BestModel]*WaveSpec[len(SpecModelNew)-1]**Index[BestModel])
        #figspec.text(WaveSpec[0],Ymin,"$\chi^2$="+'%.2f' % Chi2Best+"   WD: "+ModelFile[BestModel],{'color': 'black', 'fontsize': 14})   
        figspec.set_ylabel("Flux", size=14)
        figspec.tick_params(labelbottom=False)    
        
        figres.plot(WaveSpec, 100*(Spectrum-SpecModelNew)/Spectrum, "r|")
        figres.plot(WaveSpec, array(Spectrum)*0.0, "k:")
        figres.set_xlabel("Wavelength ($\mathrm{\AA}$)", size=14)
        figres.set_ylabel("Residuals ($\mathrm{\%}$)", size=14)
        plt.savefig(ResultDir+'/'+SpecFile+'.png', bbox_inches='tight')        
        plt.show()
    else:
        ResultFile=ResultDir+'/'+'Res.'+SpecFile
        WriteData5(len(ModelFileBest),ModelFileBest,Chi2best,NormFactorBest,
                   PowerFactorBest,PowerIndexBest,ResultFile)   
        
    print("This is the End!                        ")
    exit()

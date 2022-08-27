#!/usr/bin/env python3
import sys
from numpy import *
from pylab import *
from matplotlib import *
from sys import stdin

SDSS = ["u'","g'","r'","i'","z'"]
PS1 =  ["g", "r", "i", "z", "y"]
Bessell=["U","B","V","R","Rc","I","Ic","J","H","Ks","K"]
GALEX= ['FUV','NUV']
WISE= ['W1','W2','W3','W4']
AllBands = SDSS + PS1 + Bessell + GALEX + WISE

#########################################################################
def Flux(VegaSystem,ZP,Wave,WaveErr,Magnitude,MagErr):
    """
    Converting Magnitudes to Fluxes
    """
    if VegaSystem:
        #Magnitude = -2.5*math.log10(Flux / ZP)
        F = ZP*10**(-0.4*Magnitude)
        FErr = ZP/2.512**Magnitude*math.log(2.512)*MagErr
        #print("Vega System")
    else:
        #Magnitude = -5.*math.log10(ZP*1.e-8) -2.5*math.log10(Flux) - 42.408
        FluxJy=10**(Magnitude/(-2.5))*3631
        FluxJyErr=MagErr*10**(Magnitude/(-2.5))*3631*math.log(10)/2.5
        F=3e-5*FluxJy/Wave**2
        #Flux=10**(-(Magnitude-8.90)/2.5)
        FErr=MagErr*10**(Magnitude/(-2.5))*3631*3e-5*math.log(10)/2.5/Wave**2       
        #print("AB System")
    return F,FErr


#SDSS
#2*Col("b")*Col("FluxZero")*sinh(Col("Mag")/(-2.5/ln(10))-ln(Col("b")))
#abs((Col("FluxJy")*Col("MagErr")/(-2.5/ln(10)))/tanh((Col("Mag")/(-2.5/ln(10)))-ln(Col("b"))))
#3e-5*Col("FluxJy")/Col("Lambda")^2
#3e-5*Col("FluxJyErr")/Col("Lambda")^2

#Col(D)/2.512^Col(E)
#ZP/2.512**Magnitude*math.log(2.512)*MagErr


#########################################################################
def Convert(Band,Mag,MagErr):
    VegaSystem = None
    if Band in SDSS:
        i = SDSS.index(Band)
        #print("SDSS")        
        VegaSystem = False
        FilterZPs = [3521,4807,6253,7667,9113]
        Waves = [3521,4807,6253,7667,9113]
        WaveErrs = [320,700,680,740,700]
        ZP=FilterZPs[i]
        Wave=Waves[i]
        WaveErr=WaveErrs[i]

#FluxZero b
#3540	300	u	3767	1.4E-10
#4770	690	g'	3631	9E-11
#6222	690	r'	3631	1.2E-10
#7632	770	i'	3631	1.8E-10
#9049	690	z'	3565	7.4E-10
        
    elif Band in PS1:
        i = PS1.index(Band)
        #print("PS1")
        VegaSystem = False
        FilterZPs = [4810.,6156.,7504.,8669.,9613.]
        Waves = [4810.,6156.,7504.,8669.,9613.]
        WaveErrs = [525.,625.,600.,500.,320.]
        ZP=FilterZPs[i]
        Wave=Waves[i]
        WaveErr=WaveErrs[i]

    elif Band in GALEX:
        i = GALEX.index(Band)
        #print("PS1")
        VegaSystem = False
        FilterZPs = [1528.,2271.]
        Waves = [1528.,2271.]
        WaveErrs = [200.,500.]
        ZP=FilterZPs[i]
        Wave=Waves[i]
        WaveErr=WaveErrs[i]

    elif Band in Bessell:
        i = Bessell.index(Band)
        #print("Johnson-Cousins (Bessell)")
        VegaSystem = True
        FilterZPs = [417.5e-11,632.0e-11,363.1e-11,217.7e-11,217.7e-11,112.6e-11,112.6e-11,31.47e-11,11.38e-11,4.28e-11,4.09e-11]
        Waves = [3663,4380,5450,6410,6470,7865,7980,12200,16300,21590,21900]
        WaveErrs = [325,445,420,800,800,770,770,810,1255,1310,1620]
        ZP=FilterZPs[i]
        Wave=Waves[i]
        WaveErr=WaveErrs[i]

    elif Band in Bessell:
        i = Bessell.index(Band)
        #print("Johnson-Cousins (Bessell)")
        VegaSystem = True
        FilterZPs = [417.5e-11,632.0e-11,363.1e-11,217.7e-11,217.7e-11,112.6e-11,112.6e-11,31.47e-11,11.38e-11,4.28e-11,4.09e-11]
        Waves = [3663,4380,5450,6410,6470,7865,7980,12200,16300,21590,21900]
        WaveErrs = [325,445,420,800,800,770,770,810,1255,1310,1620]
        ZP=FilterZPs[i]
        Wave=Waves[i]
        WaveErr=WaveErrs[i]

    elif Band in WISE:
        i = WISE.index(Band)
        #print("Johnson-Cousins (Bessell)")
        VegaSystem = True
        FilterZPs = [8.1787E-12,2.415E-12,6.52e-14,5.09e-15]
        Waves = [33526.,46028.,115608.,220883.]
        WaveErrs = [3313.,5211.,27527.,20500.]
        ZP=FilterZPs[i]
        Wave=Waves[i]
        WaveErr=WaveErrs[i]

#4380	470	B	6.32E-9
#5450	425	V	3.631E-9
#6410	800	R	2.177E-9
#7980	745	I	1.126E-9
#12350	810	J	3.129E-10
#16620	1255	H	1.133E-10
#21590	1310	K	4.283E-11
#33526	--	W1	8.1787E-12
#46028	--	W2	2.415E-12


    #print(VegaSystem,ZP,Wave,WaveErr,Mag,MagErr)
    F,FErr=Flux(VegaSystem,ZP,Wave,WaveErr,Mag,MagErr)
    return Wave,WaveErr,F,FErr

#########################################################################
def ReadData(lines):
    Waves = []
    WaveErrs = []
    Fluxes = []
    FluxErrs = []
    for line in lines:        
        sublines = line.split()
        #print(len(sublines))
        Band = sublines[0]
        if len(sublines)>1:
            Mag = float(sublines[1])
#            print(Teff,'\r',)
        if len(sublines)>2:
            MagErr = float(sublines[2])
        else:
            MagErr = 0.0
        if Band in AllBands:
            Wave,WaveErr,Flux,FluxErr=Convert(Band,Mag,MagErr)
            Waves.append(Wave)
            WaveErrs.append(WaveErr)
            Fluxes.append(Flux)
            FluxErrs.append(FluxErr)
        else:
            print(Band,"is not recognised. Skipped.")
    return Waves,WaveErrs,Fluxes,FluxErrs

#########################################################################
def WriteData4(nn,aa,bb,cc,dd,output_file_path):
    """
    Write four columns of data to an external ASCII text file
    """
    outfile = open(output_file_path,"w")
    for i in range (0, nn):
        outfile.write(' %8.1f \t %8.1f \t %12.6e \t %12.6e \n' %  (aa[i],bb[i],cc[i],dd[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close()
#########################################################################

def print_header():
    print ("")
    print ("*******************************************************************************")
    print ("**                               sedplot.py                                  **")
    print ("** A little utility to plot a SED using multicolour photometric measurements **")
    print ("**                                                                           **")
    print ("**                            2022-August-28                                 **")
    print ("**                           Vitaly Neustroev                                **")
    print ("*******************************************************************************")
    print ("")

def usage():
    #print_header()
    "Usage function"
    print ("Usage: %s [options] FileName1 [FileName2] .. [FileName7]" % sys.argv[0])
    print (" ")
    print ("FileName is a SED filename (up to 7 files can be read).")
    print ("Each SED must give a band name in the first column and a magnitude in the second")
    print ("     (a magnitude error in the 3rd column is optional but will be used if given)")
    print ("Options: -hblws")
    print ("     -h: Help")
    print ("     -b: List of recognized photometric bands")
    print ("     -l: Plot in logariphmic scale [default: in linear]")
    print ("     -w: SEDs in fluxes will be written to text-files with names [FileName].sed consisted of 4 columns:") 
    print ("                                                                     Wavelength, WaveErr, Flux, FluxErr")
    print ("     -s: Plot of the sed-file [4 columns: Wavelength, WaveErr, Flux, FluxErr]")
    print ("")
    sys.exit(-1)


##########################################################################

print_header()
j=0
isLog = False
isWrite = False
isSED = False
FileList = []
pyplot.figure(figsize=(12, 8))
if len(sys.argv) == 1:
    usage()

for i in range(len(sys.argv)-1):
    CmdLinePar = sys.argv[i+1]
    if CmdLinePar[0] == '-' or CmdLinePar[0] == '+':
        if ('h' or 'H') in CmdLinePar[1:]:
            usage()
            exit()
        if ('l' or 'L') in CmdLinePar[1:]:
            isLog = True
        if ('w' or 'W') in CmdLinePar[1:]:
            isWrite = True
        if ('b' or 'B') in CmdLinePar[1:]:
            print('SDSS:')
            print(SDSS)
            print('Panstarrs-PS1:')
            print(PS1)
            print('Johnson-Cousins (Bessell):')
            print(Bessell)
            print('GALEX:')
            print(GALEX)
            print('WISE:')
            print(WISE)
            exit()
        if ('s' or 'S') in CmdLinePar[1:]:
            isSED = True
    else:
        FileName = CmdLinePar
        #print("Name: ",FileName)
        if not os.path.isfile(FileName):
            print("The File ",FileName," doesn't exist.")
            #exit(-1)               
        else:
            j+=1
            if isSED:
                data0 = loadtxt(FileName, usecols=[0,1,2,3], unpack=True)
                Waves=data0[0,:]
                WaveErrs=data0[1,:]
                Fluxes=data0[2,:]
                FluxErrs=data0[3,:]
            else:
                infile = open(FileName, "r")
                lines = infile.readlines()
                infile.close()
                Waves,WaveErrs,Fluxes,FluxErrs = ReadData(lines)
            if j==1:
                pyplot.errorbar(Waves,Fluxes,xerr=WaveErrs, yerr=FluxErrs,fmt='o',color='b',label=FileName)
            if j==2:                                                                              
                pyplot.errorbar(Waves,Fluxes,xerr=WaveErrs, yerr=FluxErrs,fmt='o',color='r',label=FileName)
            if j==3:                                                                              
                pyplot.errorbar(Waves,Fluxes,xerr=WaveErrs, yerr=FluxErrs,fmt='o',color='g',label=FileName)
            if j==4:                                                                              
                pyplot.errorbar(Waves,Fluxes,xerr=WaveErrs, yerr=FluxErrs,fmt='o',color='k',label=FileName)
            if j==5:                                                                              
                pyplot.errorbar(Waves,Fluxes,xerr=WaveErrs, yerr=FluxErrs,fmt='o',color='c',label=FileName)
            if j==6:                                                                              
                pyplot.errorbar(Waves,Fluxes,xerr=WaveErrs, yerr=FluxErrs,fmt='o',color='y',label=FileName)
            if j==7:                                                                              
                pyplot.errorbar(Waves,Fluxes,xerr=WaveErrs, yerr=FluxErrs,fmt='o',color='m',label=FileName)
            if isWrite and not isSED:
                WriteData4(len(Waves),Waves,WaveErrs,Fluxes,FluxErrs,FileName+'.sed')
#print("j=",j)
if j==0:
    print()
    #sys.stdout.write('Enter the filename: ')
    print('Enter the filename: ')
    output_file_path = stdin.readline()
    FileName = output_file_path.rstrip('\n')    
    if not os.path.exists(FileName): 
        print ("The file",FileName, "doesn't exist. Stop executing!")
        exit()
    else:
        j+=1
        infile = open(FileName, "r")
        lines = infile.readlines()
        infile.close()
        Waves,WaveErrs,Fluxes,FluxErrs = ReadData(lines)
        pyplot.errorbar(Waves,Fluxes,xerr=WaveErrs, yerr=FluxErrs,fmt='o',color='b',label=FileName)       
if isLog:
    yscale('log')
    xscale('log')
xlabel("Wavelength ($\mathrm{\AA}$)", size=14)
ylabel("Flux", size=14)
pyplot.legend()
show()

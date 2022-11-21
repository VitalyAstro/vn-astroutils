#!/usr/bin/python2
########################################################################
# program: hjd.py
# author: Vitaly Neustroev
# version: 0.1
# date: May 2, 2012
# description:  
#    
#  Convert geocentric Julian date to Heliocentric Julian date.
#  The file must have two columns: time(sec) & amplitude.
#              
########################################################################

import sys
from numpy import *
from sys import stdin
from astrolibpy import helio_jd

##########################################################################

def WriteHJD(nn,aa,bb,output_file_path):
    """
    Write two columns of data to an external ASCII text file
    """   
    output_file = output_file_path.rstrip('\n') 
    outfile = open(output_file,"w")  
    for i in range (0, nn):
        outfile.write(' %15.7f \t %15.7f \n' %  (aa[i],bb[i]))
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
            print 'Invalid Number.  Enter number. '      
    return number_float
    
########################################################################


print "****************************************************"
print " "
print "hjd.py inputfile outputfile RA(H.hh) Dec(D.dd) addJD"
print " "
print "****************************************************"
print " "
if (len(sys.argv) > 1):
  FileName = sys.argv[1]
else:
  print " "
  print "Enter the input light curve filename: "
  output_file_path = stdin.readline()
  FileName = output_file_path.rstrip('\n')

if (len(sys.argv) > 2):
    FileNameOutput = sys.argv[2]
else:
    print " "
    print "Enter the output light curve filename: "
    output_file_path = stdin.readline()
    FileNameOutput = output_file_path.rstrip('\n')

if (len(sys.argv) > 3):
    try:
      RA = float(sys.argv[3])
    except ValueError:
      print 'Invalid RA. Enter Right ascension of object for epoch 2000.0 (hours):'
      RA=enter_float()
else:
    print " "
    print "Enter Right ascension of object for epoch 2000.0 (hours):"
    RA=enter_float()

if (len(sys.argv) > 4):
    try:
      Dec = float(sys.argv[4])
    except ValueError:
      print 'Invalid Dec. Enter Declination of object for epoch 2000.0 (degrees):'
      Dec=enter_float()
else:
    print " "
    print "Enter Declination of object for epoch 2000.0 (degrees):"
    Dec=enter_float()

if (len(sys.argv) > 5):
    try:
      addJD = float(sys.argv[5])
    except ValueError:
      print 'Invalid addJD. Enter a constant to be added to Times in order to have full JDs:'
      addJD=enter_float()
else:
    print " "
    print "Enter a constant to be added to Times in order to have full JDs:"
    addJD=enter_float()

########################################################################

    
data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
Time=data0[0,:]
Flux=data0[1,:]
NumData = len(Time)

for j in range(NumData):
  hjd= helio_jd(Time[j]+addJD-2400000.,RA*15.,Dec)
  Time[j]=hjd+2400000.

WriteHJD(NumData,Time,Flux,FileNameOutput)



########################################################################


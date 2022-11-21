#!/usr/bin/python

#import numpy
#from numpy import *
import sys
import re
from sys import stdin
import numpy as np
#from pylab import *
#import pyfits
from astropy.io import fits as pyfits
from shutil import copy2
from os import rename 


##########################################################################
def print_header():
    print ("")
    print ("****************************************************************************************")
    print ("**                                  xs_reformat.py                                    **")
    print ("**   A little utility to rename FITS-files of spectra, obtained with ESO/X-shooter    **")
    print ("** or convert them into text-files, with wavelength correction (from nm to angstroms) **")
    print ("**                                    2022-Nov-21                                     **")
    print ("**                                 Vitaly Neustroev                                   **")
    print ("****************************************************************************************")
    print ("")


def usage():
    print_header()
    "Usage function"
    print ("Usage: %s [options] FileList" % sys.argv[0])
    print (" ")
    print ("FileList is a list of valid FITS files obtained with ESO/X-shooter")
    print ("Options:")
    print ("     -h: Help")
    print ("     -t: [Default] Convert FITS-files to text-files")
    print ("     -f: Copy of original FITS-files into current directory, giving them new names and")
    print ("         converting the wavelengths from nm to angstroms.")
    print ("     -c: Copy of original FITS-files into current directory, giving them new names and")
    print ("         keeping the original wavelengths unchanged.")
    print ("")
    sys.exit(-1)


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
    Write four columns of data to an external ASCII text file
    """
    outfile = open(output_file_path,"w")
    for i in range (0, nn):
        outfile.write(' %12.6f \t %12.6e \t %12.6e \t %12d \n' %  (aa[i],bb[i],cc[i],dd[i]))
        # outfile.write(' %12.6f \t %12.6f \t %12.6f \n' %  (aa[i],bb[i],cc[i]))
    outfile.close()
#########################################################################

def correctFITS(file_name,NewFileName):

    '''
    Read XShooter's fits files.
    :param file_name: name of the file in fit extension (string)
    :return: parameters
    '''

    hdulist = pyfits.open(file_name)
    header1 = hdulist[0].header
    fluxdata = hdulist[0].data
    
    header2 = hdulist[1].header
    errdata = hdulist[1].data

    header3 = hdulist[2].header
    qdata = hdulist[2].data
    
    header1['CRVAL1'] = header1['CRVAL1'] * 10.
    header1['CDELT1'] = header1['CDELT1'] * 10.
    header1['CUNIT1'] = 'angstroms'
    header2['CRVAL1'] = header2['CRVAL1'] * 10.
    header2['CDELT1'] = header2['CDELT1'] * 10.
    header2['CUNIT1'] = 'angstroms'
    header3['CRVAL1'] = header3['CRVAL1'] * 10.
    header3['CDELT1'] = header3['CDELT1'] * 10.
    header3['CUNIT1'] = 'angstroms'
    
    try:
        hdulist.writeto(NewFileName,checksum=True)
    except OSError:
        rename(NewFileName,NewFileName+'.bak')
        print('The existing file',NewFileName,'was renamed to',NewFileName+'.bak')
        try:
            hdulist.writeto(NewFileName,checksum=True)
        except OSError:
            print('Could not write the new file',NewFileName)

#########################################################################

def read_fits(file_name):

    '''
    Read XShooter's fits files.
    :param file_name: name of the file in fit extension (string)
    :return: parameters
    '''

    hdulist = pyfits.open(file_name)
    header1 = hdulist[0].header
    fluxdata = hdulist[0].data
    obj = header1['OBJECT']
    mjd = header1['MJD-OBS']
    obs_date = header1['DATE-OBS']
    arm = header1['ESO SEQ ARM']
    instrument = header1['INSTRUME']
    
    header2 = hdulist[1].header
    errdata = hdulist[1].data

    header3 = hdulist[2].header
    qdata = hdulist[2].data
    
    flux = np.copy(fluxdata)
    wave = np.ones(len(flux)) * header1['CRVAL1']
    wave = 10 * (wave + header1['CDELT1'] * np.arange(len(flux)))  # Convert to angstrom
    err = np.copy(errdata)
    qual = np.copy(qdata)
    
    return obj, instrument, arm, obs_date, mjd, wave, flux, err, qual




##############################



def ReadSpec(FileName,FITSformat):

    obj, instrument, arm, obs_date, mjd, wave, flux, err, qual = read_fits(FileName)
    mjdstr = "%.7f" %mjd   
    NewFileName = obj+ '_' + arm + '_MJD'+ mjdstr
    NewFileName = re.sub(r"\s+", '_', NewFileName)
    if FITSformat and isJustCopy:
        copy2(FileName,'./'+NewFileName+'.fits')
    elif FITSformat and not isJustCopy:
        correctFITS(FileName,'./'+NewFileName+'.fits')
    else:
        WriteData4(len(wave),wave,flux,err,qual,NewFileName+'.txt')
    #print("\nThe spectrum has %d points\n" % len(wave))
    
 
########################################################################################
#                               MAIN
########################################################################################

if len(sys.argv) == 1:
    usage()

#print len(sys.argv)
#print sys.argv
FileList = sys.argv[1]

FITSformat = False
if FileList == '-h':
    usage()
elif FileList == '-c':
    FITSformat = True
    isJustCopy = True
elif FileList == '-f':
    FITSformat = True
    isJustCopy = False
elif FileList == '-t':
    FITSformat = False
else:
    try:
        infile = open(FileList, "r")
        lines = infile.readlines()
        infile.close()      
    except NameError:
        print(" ")
        print("Enter the filename of the list of spectra: ")
        output_file_path = stdin.readline()
        FileList = output_file_path.rstrip('\n')
        try:
            infile = open(FileList, "r")
            lines = infile.readlines()
            infile.close()      
        except NameError:
            usage()
if (len(sys.argv) > 2):
    FileList = sys.argv[2]
else:
  print(" ")
  print("Enter the filename of the list of spectra: ")
  output_file_path = stdin.readline()
  FileList = output_file_path.rstrip('\n')
  try:
    infile = open(FileList, "r")
    lines = infile.readlines()
    infile.close()      
  except NameError:
      usage()


infile = open(FileList, "r")
lines = infile.readlines()
infile.close()

for line in lines:
    FileName = line.rstrip('\n')
    ReadSpec(FileName,FITSformat)

#!/usr/bin/python3

import sys
import re
from sys import stdin
import numpy as np
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
    print ("**                                    2022-Nov-22                                     **")
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
    print ("     -mjd: [Default] The file names of the converted files will be given in the MJD times.")
    print ("     -hjd: The file names of the converted files will be given in the HJD times (-2450000).")
    print ("     -bjd: The file names of the converted files will be given in the BJD times (-2450000).")
    print ("     -jd: The file names of the converted files will be given in the JD times (-2450000).")
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

def read_fits(file_name,isHJD,isBJD):

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
    exptime = header1['EXPTIME']
    ra  = header1['RA']
    dec = header1['DEC']
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
    
    jd = mjd + 2400000.5
    timeconv = jd - 2450000
    
    if isHJD or isBJD:
        from astropy import time, coordinates as coord, units as u    
        target = coord.SkyCoord(ra, dec, unit='deg', frame='icrs')
        paranal = coord.EarthLocation.of_site('paranal')
        time_corr = time.Time(jd+exptime/172800, format='jd', scale='utc', location=paranal)
        if isHJD:
            ltt_helio = time_corr.light_travel_time(target, 'heliocentric')            
            timeconv = (time_corr + ltt_helio).value - 2450000
        elif isBJD:
            ltt_bary =  time_corr.light_travel_time(target)
            timeconv = (time_corr + ltt_bary).value - 2450000
    
    return obj, instrument, arm, mjd, timeconv, wave, flux, err, qual

##############################



def ReadSpec(FileName,FITSformat,isHJD,isBJD,isJD):

    obj, instrument, arm, mjd, timeconv, wave, flux, err, qual = read_fits(FileName,isHJD,isBJD)
    if isHJD:
        #hjd = helio_jd(mjd+2400000.+exptime/172800,ra,dec)
        timestr = '_HJD'+ "%.6f" %timeconv
    elif isBJD:
        timestr = '_BJD'+ "%.6f" %timeconv
    elif isJD:
        timestr = '_JD'+ "%.6f"  %timeconv
    else:
        timestr = '_MJD'+ "%.6f" %mjd
    NewFileName = obj+ '_' + arm + timestr
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

isJD  = False
isHJD = False
isBJD = False
isFileList=False
isJustCopy = False
FITSformat = False
lines = []
s1 = 0

for i in range(len(sys.argv)-1):
    CmdLinePar = sys.argv[i+1]
    if (CmdLinePar[0] != '-'):
        if not isFileList:
            FileName = CmdLinePar
            if FileName[0] == '@':
                s1 = 1
                try:
                    infile = open(FileName[s1:], "r")
                    lines = infile.readlines()
                    infile.close()
                    #print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
                    isFileList=True                
                except:
                    print("Something wrong with the FileList ",FileName[s1:])
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
            else:
                FileName = CmdLinePar.rstrip('\n')
                lines.append(FileName)
    elif CmdLinePar[0] == '-':
        if CmdLinePar == '-h':
            usage()
            exit()
        if CmdLinePar == '-hjd':
            isHJD = True
            isBJD = False
            isJD = False
        if CmdLinePar == '-bjd':
            isHJD = False
            isBJD = True
            isJD = False
        if CmdLinePar == '-jd':
            isHJD = False
            isBJD = False
            isJD = True
        if CmdLinePar == '-c':
            FITSformat = True
            isJustCopy = True
        if CmdLinePar == '-f':
            FITSformat = True
            isJustCopy = False
        if CmdLinePar == '-t':
            FITSformat = False

if isHJD and not isBJD and not isJD:
    print('The file names will be given in the HJD times (-2450000).')
elif isBJD and not isHJD and not isJD:
    print('The file names will be given in the BJD times (-2450000).')
elif isJD and not isHJD and not isBJD:
    print('The file names will be given in the JD times  (-2450000).')
elif not isBJD and not isHJD and not isJD:
    print('The file names will be given in the MJD times.')
print()
            
for line in lines:
    FileName = line.rstrip('\n')
    ReadSpec(FileName,FITSformat,isHJD,isBJD,isJD)

#!/usr/bin/python
import matplotlib
import sys
import os
from sys import stdin
from numpy import *
import matplotlib.cm as cm
import matplotlib.pyplot as pyplot

def print_header():
    print ("")
    print ("***********************************************************************************")
    print ("**                                 vn-plot.py                                    **")
    print ("**   A little utility to plot simple line and scatter graphics from text files   **")
    print ("**                                                                               **")
    print ("**                               2025-January-23                                 **")
    print ("**                              Vitaly  Neustroev                                **")
    print ("***********************************************************************************")

def usage():
    #print_header()
    "Usage function"
    print ("Usage: %s [options] FileName1 [FileName2] .. [FileName7]" % sys.argv[0])
    print (" ")
    print ("FileName is a filename of a text file (up to 7 files can be read).")
    print ("     The file must have at least 2 columns")
    print ("     (the 3rd column can be used for error-bars, use -e option for it)")
    print ("Options: -hlse")
    print ("     -h: Help")
    print ("     -l: Plot in logariphmic scale [default: in linear]")
    print ("     -s: Scatter Plot [default: Line plot]")
    print ("     -e: Error bars [default: without error bars]")
    print ("")
    #sys.exit(-1)


##########################################################################

colors = ['k','b','r','g','c','y','m']
print_header()
j=0
isLog = False
isScatter = False
isErr = False
Scatter = ''
#FileList = []
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
        if ('s' or 'S') in CmdLinePar[1:]:
            isScatter = True
            Scatter = 'o'
        if ('e' or 'E') in CmdLinePar[1:]:
            isErr = True
    else:
        FileName = CmdLinePar
        isErrUsed = False
        if not os.path.isfile(FileName):
            print("The File ",FileName," doesn't exist.")
            #exit(-1)               
        else:
            j+=1
#            isErr = False
            try:
                data0 = loadtxt(FileName, unpack=True,skiprows=0)
                X0=data0[0,:]
                Y0=data0[1,:]
            except:
                print("Something wrong with the file ",FileName)
                print("Skiping it...")
                continue
            if isErr:
                try:
                   YErr=data0[2,:]
                   isErrUsed = True
#                   print ("The file has the 3rd column. We assume that ALL the files include error data.")
                except IndexError:    
                   print ("\nThe file doesn't have the 3rd column. ")
                   isErrUsed = False                  
            if isErrUsed:
                pyplot.errorbar(X0,Y0,yerr=YErr,fmt=Scatter,color=colors[j-1],label=FileName)
            else:
                pyplot.plot(X0, Y0,colors[j-1]+Scatter,label=FileName)

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
        try:
            data0 = loadtxt(FileName, unpack=True,skiprows=0)
            X0=data0[0,:]
            Y0=data0[1,:]
        except:
            print("Something wrong with the file ",FileName)
            exit()
        j+=1
        pyplot.plot(X0, Y0,colors[j-1]+Scatter,label=FileName)       
if isLog:
    pyplot.yscale('log')
    pyplot.xscale('log')

pyplot.xlabel("X", size=14)
pyplot.ylabel("Y", size=14)
pyplot.legend()
pyplot.show()


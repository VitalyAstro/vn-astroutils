#!/usr/bin/env python3
#import numpy, matplotlib
import sys
from numpy import *
from pylab import *
from matplotlib import *
#import matplotlib.axes3d as p3
import matplotlib.cm as cm
#import numpy
import matplotlib.pyplot as pyplot
#from matplotlib.mlab import griddata
#from enthought.mayavi import mlab
#from enthought.mayavi.mlab import *    !!!!!!!!!!!!!!!11


if (len(sys.argv) > 1):
    pyplot.figure()
    FileName = sys.argv[1]
    data0 = loadtxt(FileName, usecols=[0,1,2], unpack=True)
    X0=data0[0,:]
    Y0=data0[1,:]
    Y1=data0[2,:]
    pyplot.plot(X0, Y0,color='b',label="Column 2")
    pyplot.plot(X0, Y1,color='r',label="Column 3")
else:
    print("*****************************************")
    print(" ")
    print("plot3Col.py infile")
    print("     infile: a 3 column file XYY")
    print(" ")
    print("*****************************************")

    
#if (len(sys.argv) > 2):
    #FileName = sys.argv[2]
    #data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    #X0=data0[0,:]
    #Y0=data0[1,:]
    #pyplot.plot(X0, Y0,color='r',label=sys.argv[2])
##    legend.add(sys.argv[2])
#if (len(sys.argv) > 3):
    #FileName = sys.argv[3]
    #data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    #X0=data0[0,:]
    #Y0=data0[1,:]
    #pyplot.plot(X0, Y0,color='g',label=sys.argv[3])
##    legend.add((sys.argv[3]))
#if (len(sys.argv) > 4):
    #FileName = sys.argv[4]
    #data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    #X0=data0[0,:]
    #Y0=data0[1,:]
    #pyplot.plot(X0, Y0,color='bl',label=sys.argv[4])
#    legend.add((sys.argv[4]))

xlabel("X")
ylabel("Y")
pyplot.legend()
show()

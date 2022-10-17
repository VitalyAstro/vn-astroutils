#!/usr/bin/python3
#import numpy, matplotlib
import sys
from numpy import *
#from pylab import *
from matplotlib import *
#import matplotlib.axes3d as p3
import matplotlib.cm as cm
#import matplotlib.pyplot as plt
#from matplotlib.mlab import griddata
#from enthought.mayavi import mlab
#from enthought.mayavi.mlab import *    !!!!!!!!!!!!!!!11


if (len(sys.argv) > 1):
    pyplot.figure()
    FileName = sys.argv[1]
    data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    X0=data0[0,:]
    Y0=data0[1,:]
    pyplot.plot(X0, Y0,color='b',label=sys.argv[1])
if (len(sys.argv) > 2):
    FileName = sys.argv[2]
    data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    X0=data0[0,:]
    Y0=data0[1,:]
    pyplot.plot(X0, Y0,color='r',label=sys.argv[2])
#    legend.add(sys.argv[2])
if (len(sys.argv) > 3):
    FileName = sys.argv[3]
    data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    X0=data0[0,:]
    Y0=data0[1,:]
    pyplot.plot(X0, Y0,color='g',label=sys.argv[3])
#    legend.add((sys.argv[3]))
if (len(sys.argv) > 4):
    FileName = sys.argv[4]
    data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    X0=data0[0,:]
    Y0=data0[1,:]
    pyplot.plot(X0, Y0,color='k',label=sys.argv[4])
#    legend.add((sys.argv[4]))
if (len(sys.argv) > 5):
    FileName = sys.argv[5]
    data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    X0=data0[0,:]
    Y0=data0[1,:]
    pyplot.plot(X0, Y0,color='m',label=sys.argv[5])
#    legend.add((sys.argv[4]))
if (len(sys.argv) > 6):
    FileName = sys.argv[6]
    data0 = loadtxt(FileName, usecols=[0,1], unpack=True)
    X0=data0[0,:]
    Y0=data0[1,:]
    pyplot.plot(X0, Y0,color='c',label=sys.argv[6])
#    legend.add((sys.argv[4]))


# Legend the plot
#title("Phase diagram")
xlabel("X")
ylabel("Y")
pyplot.legend()
#leg1 = sys.argv[1]
#    legend((sys.argv[1]))
#legend((sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))
#legend()
show()

# define grid.
#xi = np.linspace(-0.10,1.40,100)
#yi = np.linspace(-0.75,0.75,100)

# grid the data.
#mi = griddata(X,Y,M,xi,yi)
#ti = griddata(X,Y,T,xi,yi)

# contour the gridded data, plotting dots at the randomly spaced data points.
#contour(xi,yi,ti,12,linewidths=0.5,colors='k')
#contourf(xi,yi,zi,12,cmap=cm.jet)
#surf(xi,yi,ti)
#mesh(xi, yi, ti)
#colorbar() # draw colorbar

#xlim(-0.10,1.40)
#ylim(-0.75,0.75)
#zlim(-0.75,0.75)

#!!!!!!!!!!!!!!!!!!!!!
#points3d(X, Y, Z, T, colormap="copper", scale_factor=0.03,vmin=3000)
#points3d(X, Y, Z, M, colormap="copper", scale_factor=0.03)

#Draws lines between points.
#plot3d(X, Y, Z, T, colormap="copper",representation="points",vmin=1000)

#pcolor(xi, yi, ti, vmin=3000)
# plot data points.
#scatter(X,Y,marker='o',c='b',s=5)
#scatter(X,Y,marker='o')


#plt.show()

#plot(X,Y,'.')
#pcolor(X, Y, T)

#plt.subplots_adjust(hspace=0.5)
#plt.subplot(121)
#hexbin(X,Y,T, cmap=cm.jet)

#!!!!!!!!!!!!!!!!!!!!!
#axes(z_axis_visibility=True)
#axes(z_axis_visibility=True,extent=[-0.10,1.40,-0.75,0.75,-0.75,0.75])
#show()

#print "Press ^C!"
#while True:
#	pass

#scatter(X,Y,marker='o',c='b',s=5)
#plot(X,Y,marker='o')
#show()


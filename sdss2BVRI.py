#!/usr/bin/env python
""" Transformations between SDSS/PS1/VST and UBVRI magnitudes

"""
import sys
from sys import stdin

########################################################################
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
            print('Invalid Number.  Enter number. ')
    return number_float
    
########################################################################

def print_header():
    print ("")
    print ("***********************************************************************************************")
    print ("**                                     sdss2BVRI.py                                          **")
    print ("**        A utility to transform SDSS/PS1/VST magnitudes to UBVRI magnitudes                 **")
    print ("**                                     2022-Oct-07                                           **")
#2022-Oct-07: python2 --> python3; format of results.
    print ("**                                   Vitaly Neustroev                                        **")
    print ("***********************************************************************************************")
    print ("")

########################################################################

print_header()
#print('Transformations between SDSS/PS1/VST and UBVRI magnitudes')
#print()
print('Enter the u magnitude (or a dummy number): ')
u=enter_float()
print('Enter the g magnitude: ')
g=enter_float()
print('Enter the r magnitude: ')
r=enter_float()
print('Enter the i magnitude: ')
i=enter_float()
print('Enter the z magnitude: ')
z=enter_float()

#print("Entered numbers:")
#print('SDSS u=',"{:7.3f}".format(u))
#print('SDSS u=',"{:7.3f}".format(g))
#print('SDSS u=',"{:7.3f}".format(r))
#print('SDSS u=',"{:7.3f}".format(i))
#print('SDSS u=',"{:7.3f}".format(z))

B1 = u - 0.8116*(u - g) + 0.1313
B2 = g + 0.3130*(g - r) + 0.2271

V1 = g - 0.2906*(u - g) + 0.0885
V2 = g - 0.5784*(g - r) - 0.0038

R1 = r - 0.1837*(g - r) - 0.0971
R2 = r - 0.2936*(r - i) - 0.1439

I1 = r - 1.2444*(r - i) - 0.3820
I2 = i - 0.3780*(i - z) - 0.3974

U1 = 0.78*(u-g) - 0.88 + (B1+B2)/2.
U2 = 0.52*(u-g) + 0.53 * (g-r) - 0.82 + (B1+B2)/2.

print ("")
print('Johnson magnitudes, converted from SDSS magnitudes:')
print('U magnitudes: U1=',"{:7.3f}".format(U1),'\tU2=',"{:7.3f}".format(U2))
print('B magnitudes: B1=',"{:7.3f}".format(B1),'\tB2=',"{:7.3f}".format(B2))
print('V magnitudes: V1=',"{:7.3f}".format(V1),'\tV2=',"{:7.3f}".format(V2))
print('R magnitudes: R1=',"{:7.3f}".format(R1),'\tR2=',"{:7.3f}".format(R2))
print('I magnitudes: I1=',"{:7.3f}".format(I1),'\tI2=',"{:7.3f}".format(I2))

uo = u + 0.01*(u - g) - 0.03
go = g + 0.05*(g - r) - 0.06
ro = r + 0.03*(g - r) - 0.035
io = i - 0.015
zo = z - 0.04*(i - z) + 0.05

B1 = uo - 0.8116*(uo - go) + 0.1313
B2 = go + 0.3130*(go - ro) + 0.2271

V1 = go - 0.2906*(uo - go) + 0.0885
V2 = go - 0.5784*(go - ro) - 0.0038

R1 = ro - 0.1837*(go - ro) - 0.0971
R2 = ro - 0.2936*(ro - io) - 0.1439

I1 = ro - 1.2444*(ro - io) - 0.3820
I2 = io - 0.3780*(io - zo) - 0.3974

U1 = 0.78*(uo - go) - 0.88 + (B1 + B2)/2.
U2 = 0.52*(uo - go) + 0.53 * (go - ro) - 0.82 + (B1 + B2)/2.

print ("")
print('Johnson and SDSS magnitudes, converted from VST magnitudes:')
print('SDSS u=',"{:7.3f}".format(uo),'   U magnitudes: U1=',"{:7.3f}".format(U1),'\tU2=',"{:7.3f}".format(U2))
print('SDSS g=',"{:7.3f}".format(go),'   B magnitudes: B1=',"{:7.3f}".format(B1),'\tB2=',"{:7.3f}".format(B2))
print('SDSS r=',"{:7.3f}".format(ro),'   V magnitudes: V1=',"{:7.3f}".format(V1),'\tV2=',"{:7.3f}".format(V2))
print('SDSS i=',"{:7.3f}".format(io),'   R magnitudes: R1=',"{:7.3f}".format(R1),'\tR2=',"{:7.3f}".format(R2))
print('SDSS z=',"{:7.3f}".format(zo),'   I magnitudes: I1=',"{:7.3f}".format(I1),'\tI2=',"{:7.3f}".format(I2))


B1 = g + 0.561*(g - r) + 0.194
B2 = g + 0.540*(g - r)  + 0.016*(g - r)**2 + 0.199

V1 = g - 0.508*(g - r) - 0.017
V2 = r + 0.492*(g - r) - 0.017

R1 = r - 0.166*(g - r) - 0.142
R2 = r - 0.275*(r - i) - 0.166

I1 = i - 0.167*(g - r) - 0.376
I2 = i - 0.214*(r - i) - 0.416

print ("")
print('Johnson magnitudes, converted from PANSTARRS magnitudes:')
print('B magnitudes: B1=',"{:7.3f}".format(B1),'\tB2=',"{:7.3f}".format(B2))
print('V magnitudes: V1=',"{:7.3f}".format(V1),'\tV2=',"{:7.3f}".format(V2))
print('R magnitudes: R1=',"{:7.3f}".format(R1),'\tR2=',"{:7.3f}".format(R2))
print('I magnitudes: I1=',"{:7.3f}".format(I1),'\tI2=',"{:7.3f}".format(I2))


########################################################################

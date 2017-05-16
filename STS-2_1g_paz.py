#!/bin/env python
''' this is for a sts2 1g.'''
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi as pi


# generation 1
start = 0.
stop = 500.
step = 1.
n = np.arange(start,stop+step,step)
fn = n*2
wn=2*pi*fn
wmix = -2*pi*29.8
print (wmix)

#zeros
z0=-318.6 + 401.2j
z1=-318.6 - 401.2j
z2=-15.15
#poles
p0=-7.454e3-7.142e3j
p1=-7.454e3+7.142e3j
p2=-417.1
p3=-100.9 + 401.9j
p4=-100.9 - 401.9j
p5=-15.99
p6=-0.037-0.037j
p7=-0.037+0.037j

#calc gfit
multiplier = ((wn*1j)**2)
numerator = 8.6177e15*( (1j*wn - z0)*(1j*wn - z1)* (1j*wn - z2))
denominator = (1j*wn -p0)* (1j*wn -p1) * (1j*wn -p2) * (1j*wn -p3) * (1j*wn -p4) * (1j*wn -p5) * (1j*wn -p6) * (1j*wn -p7) 
gfit= multiplier*(numerator/denominator*(1j*wn-wmix))

#plot
plt.figure()
#plt.plot(fn, np.abs(gfit))
#plt.semilogx(fn, abs(gfit*ifit))
plt.semilogx(fn, abs(gfit))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
#plt.xlim([0.2,100])
#plt.ylim([400,1600])
plt.grid()
plt.show()

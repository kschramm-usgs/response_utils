#!/bin/env python
''' based on yet another ringler matlab script.  this will take
    poles and zeros and try to find the pole of the calibration 
    coil 
'''
import numpy as np

# create some vectors
n = np.arange(0.0001,500,0.0001)
fn = n*0.2
wn = 2*np.pi*fn

# define some poles from STS2.5 documentation
z1=-2*pi*1.5;
z2=0;
z3=0;
p1=-2*pi*1.55;
p2=-2*pi*35+j*2*pi*22;
p3=-2*pi*35-j*2*pi*22;
p4=-2*pi*35+j*2*pi*109;
p5=-2*pi*35-j*2*pi*109;
p6=-0.037+j*0.037;
p7=-0.037-j*0.037;

# calculate the gfit:
gfit=1500*((j*wn-z2)*(j*wn-z3))*(p1*(abs(p2)^2)*(abs(p4)^2)/z1)*((j*wn-z1)/((j*wn-p1)*(j*wn-p2)*(j*wn-p3)*(j*wn-p4)*(j*wn-p5)*(j*wn-p6)*(j*wn-p7)));

#define the the inverse filter:
px=2*pi*180;
hx=0.5;
py=2*pi*100;

#:wq
calculate the fit of the inverse filter
ifit=(px^2+2*hx*px*i*wn-wn.^2).*(i*wn+py)/(abs(py)*px^2);

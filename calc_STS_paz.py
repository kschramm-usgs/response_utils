#!/bin/env python
''' this is from Adam Ringler's matlab code on the n drive '''

import numpy as np
import matplotlib.pyplot as plt
from numpy import pi as pi


# set up some variables:
start=1
stop=500
step=1
n=np.arange(start,stop+step,step)
fn=n*2;
wn=2*pi*fn;


#convert poles and zeros
z1=2*pi*1.5;
p1=2*pi*1.55;
p2=2*pi*35+1j*2*pi*22;
p3=2*pi*35-1j*2*pi*22;
p4=2*pi*35+1j*2*pi*109;
p5=2*pi*35-1j*2*pi*109;

#calculate the GF
gfit=1500*(p1*(abs(p2)**2)*(abs(p4)**2)/z1)*((1j*wn+z1)/((1j*wn+p1)*(1j*wn+p2)*(1j*wn+p3)*(1j*wn+p4)*(1j*wn+p5)));

px=2*pi*180;
hx=0.5;
py=2*pi*100;

#calculate the IF
ifit=(px**2+2*hx*px*1j*wn-wn**2)*(py+1j*wn)/(py*px**2);


#plot
#h, f = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

plt.figure()
#plt.subplot(121)
plt.semilogx(fn, abs(gfit))
plt.semilogx(fn, abs(gfit*ifit))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.xlim([0.2,100])
plt.ylim([400,1600])
plt.grid()


#figure(1)
#clf
#p1=semilogx(fn,abs(gfit),'Color','r','LineWidth',2);
#hold on
#p2=semilogx(fn,abs(gfit.*ifit),'Color','b','LineWidth',2);
#legend([p1 p2],'A_{fit_n}','A_{G_n}','FontSize',16);
#xlabel('Frequency (Hz)','FontSize',16);
#ylabel('Amplitude (V/m/s)','FontSize',16);
#set(gca,'FontSize',16);
#title('STS-2.5 Amplitude Response','FontSize',16);
#grid on;
#xlim([.2 100]);
#orient Landscape 
#print('-djpeg','STS25Amp.jpg');
#
plt.figure()
#clf
plt.semilogx(fn,np.angle(gfit)*180/np.pi)
#,'Color','r','LineWidth',2);
#hold on
plt.semilogx(fn,np.angle(gfit*ifit)*180/np.pi)
plt.grid()
#'Color','b','LineWidth',2);
#legend([p1 p2],'A_{fit_n}','A_{G_n}','FontSize',16);
#plt.xlabel('Frequency (Hz)','FontSize',16);
#plt.ylabel('Phase (degrees)','FontSize',16);
#set(gca,'FontSize',16);
#title('STS-2.5 Phase Response','FontSize',16);
#grid on;
plt.xlim([.2,100]);
#orient Landscape 
#print('-djpeg','STS25Phase.jpg');
plt.show()

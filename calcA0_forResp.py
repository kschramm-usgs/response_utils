#! /bin/env python
import numpy as np
import matplotlib.pyplot as plt

from obspy.signal.invsim import paz_2_amplitude_value_of_freq_resp, paz_to_freq_resp
# this is from Austin Holland's seisutils
from response import Response
print("script to calculate response normalization")


#for the STS-6
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.3068, -22.3068, -501.105, -501.105, -501.105, -501.105, -501.105, -501.105, -121.58+1j*647.231, -121.58-1j*647.231, -351.858]
#zer=[0, 0, -21.9110,-21.9911, -351.858, -605.071, -521.5+1j*960.699, -521.5-1j*960.699]
#instName='STS-6 '

#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.3068, -22.3068, -501.105, -501.105, -501.105, -501.105, -501.105, -501.105, -121.58+1j*647.231, -121.58-1j*647.231, -351.858]
# removing 2 of the 6 repeated poles:
pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.3068, -22.3068, -501.105, -501.105, -501.105, -501.105, -121.58+1j*647.231, -121.58-1j*647.231, -351.858]
zer=[0, 0, -21.9110,-21.9911, -351.858]
instName='STS-6 no coil, 4 repeats'
#
#for the STS-6B
#pol=[-0.01234-1j*0.01234, -0.01234+1j*0.01234, -255, -15.6, -97.3+1j*401, -97.3-1j*401, -520, -375, -13300, -10500-1j*10100, -10500+1j*10100]
#zer=[0 , 0, -5910+1j*3410, -5910-1j*3410, -555, -295, -10.8, -684-1j*176, -684+1j*176 ] 
#instName='STS-6B'

#zer=[0,0, 15.707963267948966, 15.707963267948966, 592.0645514955324, 592.0645514955324]
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, 22.399555620095224, 22.399555620095224, 510.006151383767, 510.006151383767, 510.006151383767, 510.006151383767]
#instName='STS-6 newer resp'



print("calculating info for a "+instName)
# use Austin's response class to build the resp from the poles and zeros
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
norm_freq=0.02
#norm_freq=1.0
n,f = inst.check_normalization(freq=norm_freq, nfft=2**24,t_sample=0.001)
scale_fac= 1.0/n
print ('The A0 norm factor is: '+str(scale_fac)+' for f='+str(norm_freq))

#check the value
inst.a0=1.0/n
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
print('This should be close to 1: '+str(1.0/n))


h, f = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

plt.figure()
plt.subplot(121)
plt.semilogx(f, abs(h))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
#plt.xlim([0.0001,10000])


plt.subplot(122)
# take negative of imaginary part
phase = np.unwrap(np.arctan2(-h.imag, h.real))
plt.semilogx(f, phase)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [radian]')
# title, centered above both subplots
plt.suptitle('Frequency Response of '+inst.desc +' Seismometer')
plt.xlim([0.0001,10000])
# make more room in between subplots for the ylabel of right plot
plt.subplots_adjust(wspace=0.3)
plt.show()

plt.show()


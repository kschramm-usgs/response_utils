#!/bin/env python
''' plot up some different poles and zeros to see how they match'''
import numpy as np
import matplotlib.pyplot as plt

from obspy.signal.invsim import paz_to_freq_resp
# this is from Austin Holland's seisutils
from response import Response
 
norm_freq = 1.0

#define the poles and zeros.
#tr-120 NRL paz
pol=[-32.55, -142, -364 + 404j, -364 - 404j]
zer=[-31.63,   -350]
instName='TR120 NRL high-f paz'
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n
h1, f1 = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

#tr-120 one high-f cal signal
pol=[-28.15385, -146.36053, -313.52396 - 429.04933j, -313.52396 + 429.04933j]
zer=[-27.63939,   -467.15426]
instName='TR120 one high-f paz'
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n
h2, f2 = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

#tr120  30 min of high f cal signal
pol=[-23.67549, -147.92203, -312.53386 - 436.54774j, -312.53386 + 436.54774j]
zer=[ -23.24459,   -489.502]
instName='TR120 30 min high-f paz'
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n
h3, f3 = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

#tr120  45 min of high f cal signal
pol=[-20.68972, -148.72745, -314.83534 - 438.03747j, -314.83534 + 438.03747j]
zer=[ -20.3227, -493.60654]
instName='TR120 45 min high-f paz'
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n
h4, f4 = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

#tr-120 2 min high-f cal signal
pol=[-35.90014, -145.49746, -314.38332 - 425.38735j, -314.38332 + 425.38735j]
zer=[-35.41853,   -458.82071]
instName='TR120 2 min high-f paz'
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n
h5, f5 = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

#tr120  60 min of high f cal signal
pol=[-18.4309, -150.02235, -315.06088 - 444.01042j, -315.06088 + 444.01042j]
zer=[-18.1278,   -512.28696]
instName='TR120 60 min high-f paz'
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n
h6, f6 = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

#tr120 no zoom of high f cal signal
pol=[-26.61705, -146.2397, -313.59657 - 428.50791j, -313.59657 + 428.50791j]
zer=[-26.0876,   -465.27966]
instName='TR120 no zoom high-f paz'
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n
h7, f7 = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)


# plot up the information
plt.figure()
plt.subplot(121)
plt.semilogx(f1, abs(h1), color='r',label='NRL response' )
plt.semilogx(f2, abs(h2), label='  Only cal window')
plt.semilogx(f7, abs(h7), label='NoZoom cal window')
plt.semilogx(f5, abs(h5), label=' 2 min cal window')
plt.semilogx(f3, abs(h3), label='30 min cal window')
plt.semilogx(f4, abs(h4), label='45 min cal window')
plt.semilogx(f6, abs(h6), label='60 min cal window')
plt.grid()
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude [units]')
plt.legend()
plt.xlim([0.2,40])

plt.subplot(122)
# take negative of imaginary part
phase1 = np.unwrap(np.arctan2(-h1.imag, h1.real))
phase2 = np.unwrap(np.arctan2(-h2.imag, h2.real))
phase3 = np.unwrap(np.arctan2(-h3.imag, h3.real))
phase4 = np.unwrap(np.arctan2(-h4.imag, h4.real))
plt.semilogx(f1, phase1, color='red')
plt.semilogx(f2, phase2)
plt.semilogx(f3, phase3)
plt.semilogx(f4, phase4)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [radian]')
plt.grid()
# title, centered above both subplots
plt.suptitle('Frequency Response different windows of high-f cal')
plt.xlim([0.2,40])
# make more room in between subplots for the ylabel of right plot
plt.subplots_adjust(wspace=0.3)
plt.show()

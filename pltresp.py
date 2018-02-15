#!/usr/bin/env python

import numpy as np
from obspy.signal.invsim import paz_2_amplitude_value_of_freq_resp, paz_to_freq_resp
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

instName='AAK 10'
# use Austin's response class to build the resp from the poles and zeros
# for ii aak
# this is the resp we have in /APPS/metadata

# this is the IDA start resp:
paz1 = {'zeros':[
                  0.00000E+00+  0.00000E+00j,
                  0.00000E+00+  0.00000E+00j,
                 -2.41120E+00+  0.00000E+00j,
                 -2.81068E+01+  0.00000E+00j,
                 -7.37047E+01+ -6.85162E+01j,
                 -7.37047E+01+  6.85162E+01j],
        'poles':[
                 -5.91568E-03+  5.86737E-03j,
                 -5.91568E-03+ -5.86737E-03j,
                 -2.45700E+00+  0.00000E+00j,
                 -1.54921E+01+ -6.37734E+01j,
                 -1.54921E+01+  6.37734E+01j,
                 -5.96513E+01+  0.00000E+00j,
                 -8.28083E+01+  0.00000E+00j,
                 -1.67590E+03+ -1.59951E+03j,
                 -1.67590E+03+  1.59951E+03j,
                 -2.11676E+03+  0.00000E+00j,
                 -4.06000E+01+  0.00000E+00j],
        'gain': 1}
paz1['gain'] = 1./paz_2_amplitude_value_of_freq_resp(paz1, .1)

# this is the fit from IDA
paz2={'zeros':[
                0.00000E+00+  0.00000E+00j,
                0.00000E+00+  0.00000E+00j,
               -3.45984E+00+  0.00000E+00j,
               -2.76605E+01+  0.00000E+00j,
               -7.37047E+01+ -6.85162E+01j,
               -7.37047E+01+  6.85162E+01j],
      'poles':[
               -5.91568E-03+  5.86737E-03j, 
               -5.91568E-03+ -5.86737E-03j,
               -3.64739E+00+  0.00000E+00j,
               -2.32375E+01+ -7.10784E+01j,
               -2.32375E+01+  7.10784E+01j,
               -3.46549E+01+  0.00000E+00j,
               -4.56665E+01+  0.00000E+00j,
               -1.67590E+03+ -1.59951E+03j,
               -1.67590E+03+  1.59951E+03j,
               -2.11676E+03+  0.00000E+00j, 
               -3.40583E+01+  0.00000E+00j],
      'gain': 1}  
paz2['gain'] = 1./paz_2_amplitude_value_of_freq_resp(paz2, .1)

# this is from sensor test suite, starting with AAK resp 
paz3 = {'zeros':[
                0.00000E+00+  0.00000E+00j,
                0.00000E+00+  0.00000E+00j,
               -3.45984E+00+  0.00000E+00j,
               -3.34238E+01+  0.00000E+00j,
               -4.65364E+02+ -3.18375E+01j,
               -4.65364E+02+  3.18375E+01j],
        'poles':[
               -5.91568E-03+  5.86737E-03j, 
               -5.91568E-03+ -5.86737E-03j,
               -3.64739E+00+  0.00000E+00j,
               -4.00447E+01+ -1.11554E+02j,
               -4.00447E+01+  1.11554E+02j,
               -5.16208E+01+  0.00000E+00j,
               -6.129300+01+  0.00000E+00j,
               -1.67590E+03+ -1.59951E+03j,
               -1.67590E+03+  1.59951E+03j,
               -2.11676E+03+  0.00000E+00j, 
               -3.22515E+01+  0.00000E+00j],
       'gain':1}
paz3['gain'] = 1./paz_2_amplitude_value_of_freq_resp(paz3, .1)

h1,f1 = paz_to_freq_resp(paz1['poles'], paz1['zeros'], paz1['gain'], 1./200., 2**18, freq=True)
h2,f2 = paz_to_freq_resp(paz2['poles'], paz2['zeros'], paz2['gain'], 1./200., 2**18, freq=True)
h3,f3 = paz_to_freq_resp(paz3['poles'], paz3['zeros'], paz3['gain'], 1./200., 2**18, freq=True)
#h4,f4 = paz_to_freq_resp(paz4['poles'], paz4['zeros'], paz4['gain'], 1./200., 2**18, freq=True)
#plotting....
plt.figure(1, figsize=(8,8))
plt.subplots_adjust(hspace=0.001)
plt.subplot(211)

print('plotting amplitude')
plt.semilogx(f1, 20*np.log10(h1), label='IDA starting resp', linewidth=6.0, alpha=0.35)
plt.semilogx(f2, 20*np.log10(h2), label='IDA method resp')
#plt.semilogx(f4, 20*np.log10(h4), label='STS 2.5 no coil', linewidth=6.0, alpha=0.35)
plt.semilogx(f3, 20*np.log10(h3), label='test suite resp')
plt.legend()
plt.ylabel('Amplitude')
#plt.xlim([100,1000])
#plt.xlim([0.0001,10000])
#plt.xlim([0.001,1])
#plt.xlim([40,800])
#plt.xlim([0.2,50])
plt.xticks()
plt.minorticks_on()
plt.grid(which='both', linestyle=':')
plt.subplot(212)
# take negative of imaginary part
print('plotting phase')
phase = np.unwrap(np.arctan2(-h1.imag, h1.real))
plt.semilogx(f1, phase, label='IDA starting resp', linewidth=6.0, alpha=0.35)
phase = np.unwrap(np.arctan2(-h2.imag, h2.real))
plt.semilogx(f2, phase, label='IDA method resp')
#phase = np.unwrap(np.arctan2(-h4.imag, h4.real))
#plt.semilogx(f4, phase, label='STS 2.5 no coil', linewidth=6.0, alpha=0.35)
phase = np.unwrap(np.arctan2(-h3.imag, h3.real))
plt.semilogx(f3, phase, label='test suite resp')
plt.legend()
#plt.xlabel('Period [S]')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase (radian)')
plt.xticks()
plt.minorticks_on()
plt.grid(which='both', linestyle=':')
# title, centered above both subplots
plt.suptitle('Response of '+instName +' Seismometer')
#plt.xlim([40,800])
#plt.xlim([100,1000])
#plt.xlim([0.001,1])
#plt.xlim([0.0001,10000])
#plt.xlim([0.2,50])
# make more room in between subplots for the ylabel of right plot
plt.subplots_adjust(wspace=0.3)
print('writing file')
figName=instName.replace(' ','_') + 'compareResp.jpg'
plt.savefig(figName, format='JPEG',dpi=400)
plt.show()

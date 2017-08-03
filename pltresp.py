#!/usr/bin/env python

import numpy as np
from obspy.signal.invsim import paz_2_amplitude_value_of_freq_resp, paz_to_freq_resp
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

instName='CTAO 00'
# use Austin's response class to build the resp from the poles and zeros
# for ctao
# this is the resp we have in /APPS/metadata

paz1 = {'zeros':[0.000000E+00+ 0.000000E+00j, 0.000000E+00+ 0.000000E+00j],'poles': 
        [-1.234000E-02+ 1.234000E-02j, -1.234000E-02+-1.234000E-02j, -3.918000E+01+ 4.912000E+01j, -3.918000E+01+-4.912000E+01j],'gain': 1}

paz1['gain'] = 1./paz_2_amplitude_value_of_freq_resp(paz1, .1)


# this is from /r02/scratch/CALS/NEW/IU_CTAO_2017_045/STS1VBBE3/ZP_CTAO_00_BHZ_2017_046_01_25_20_0000-14400
paz2={'zeros': [0.000000E+00+ 0.000000E+00j, 0.000000E+00+ 0.000000E+00j, -4.147725E-02+ 0.000000E+00j, -4.147725E-02+ 0.000000E+00j],
        'poles': [-1.703608E-02+ 1.234000E-02j,-1.703608E-02+-1.234000E-02j, -1.778991e-02+ 0.000000e-00j,-5.605308e-02+ 0.000000e-00j,
            -3.918000E+01+ 4.912000E+01j, -3.918000E+01+-4.912000E+01j],'gain': 1}  
paz2['gain'] = 1./paz_2_amplitude_value_of_freq_resp(paz2, .1)
# this is from sensor test suite
paz3 ={'zeros': [0.000000E+00+ 0.000000E+00j, 0.000000E+00+ 0.000000E+00j], 'poles':
        [-1.263000E-02+ 1.265000E-02j, -1.263000E-02+-1.265000E-02j, 
            -3.620107E+01+ 6.850121E+01j, -3.620107E+01+-6.850121E+01j],'gain':1}
paz3['gain'] = 1./paz_2_amplitude_value_of_freq_resp(paz3, .1)

h1,f1 = paz_to_freq_resp(paz1['poles'], paz1['zeros'], paz1['gain'], 1./200., 2**18, freq=True)
h2,f2 = paz_to_freq_resp(paz2['poles'], paz2['zeros'], paz2['gain'], 1./200., 2**18, freq=True)
h3,f3 = paz_to_freq_resp(paz3['poles'], paz3['zeros'], paz3['gain'], 1./200., 2**18, freq=True)
#plotting....
plt.figure(1, figsize=(8,8))
plt.subplots_adjust(hspace=0.001)
plt.subplot(211)

print('plotting amplitude')
plt.semilogx(1./f1, 20*np.log10(h1), label='metadata resp')
plt.semilogx(1./f2, 20*np.log10(h2), label='old method resp')
plt.semilogx(1./f3, 20*np.log10(h3), label='test suite resp')
plt.legend()
#plt.xlabel('Period [S]')

plt.ylabel('Amplitude')
plt.xlim([0.0001,10000])
#plt.xlim([40,800])
#plt.xlim([0.2,40])
plt.xticks([])
plt.subplot(212)
# take negative of imaginary part
print('plotting phase')
phase = np.unwrap(np.arctan2(-h1.imag, h1.real))
plt.semilogx(1/f1, phase, label='metadata resp')
phase = np.unwrap(np.arctan2(-h2.imag, h2.real))
plt.semilogx(1/f2, phase, label='old method resp')
phase = np.unwrap(np.arctan2(-h3.imag, h3.real))
plt.semilogx(1/f3, phase, label='test suite resp')
plt.legend()
plt.xlabel('Period (s)')
#plt.xlabel('Period [S]')
plt.ylabel('Phase (radian)')
# title, centered above both subplots
plt.suptitle('Response of '+instName +' Seismometer')
#plt.xlim([40,800])
plt.xlim([0.0001,10000])
# make more room in between subplots for the ylabel of right plot
plt.subplots_adjust(wspace=0.3)
print('writing file')
figName=instName.replace(' ','_') + 'compareResp.jpg'
plt.savefig(figName, format='JPEG',dpi=400)

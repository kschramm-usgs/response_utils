#! /bin/env python
import numpy as np
import matplotlib.pyplot as plt

from obspy.signal.invsim import paz_to_freq_resp
# this is from Austin Holland's seisutils
from response import Response
print("script to calculate response normalization")



# use Austin's response class to build the resp from the poles and zeros
# for ctao
# this is the resp we have in /APPS/metadata
zer1=[
         0.000000E+00+ 0.000000E+00j,
         0.000000E+00+ 0.000000E+00j,
        -9.424780E+00+ 0.000000E+00j,
        -6.283190E+02+ 0.000000E+00j,
        -5.654870E+02+ 9.794520E+02j,
        -5.654870E+02+-9.794520E+02j]
pol1=[
        -3.731580E-02+-3.670000E-02j,
        -3.731580E-02+ 3.670000E-02j,
        -9.738940E+00+ 0.000000E+00j,
        -2.199110E+02+ 1.382300E+02j,
        -2.199110E+02+-1.382300E+02j,
        -2.199110E+02+ 6.848670E+02j,
        -2.199110E+02+-6.848670E+02j]

# this is from the CTBTO
pol2=[
        -3.700000E-02+ 3.700000E-02j,
        -3.700000E-02+-3.700000E-02j,
        -9.738940E+00+ 0.000000E+00j,
        -2.199110E+02+ 1.382300E+02j,
        -2.199110E+02+-1.382300E+02j,
        -2.199110E+02+ 6.848670E+02j,
        -2.199110E+02+-6.848670E+02j]
zer2=[
        -9.424780E+00+ 0.000000E+00j,
        -6.283190E+02+ 0.000000E+00j,
        -5.654870E+02+ 9.794520E+02j,
        -5.654870E+02+-9.794520E+02j,
         0.000000E+00+ 0.000000E+00j,
         0.000000E+00+ 0.000000E+00j]

# this is from sensor test suite
zer3=[
         0.000000E+00+ 0.000000E+00j,
         0.000000E+00+ 0.000000E+00j]
pol3=[
        -1.263000E-02+ 1.265000E-02j,
        -1.263000E-02+-1.265000E-02j,
        -3.620107E+01+ 6.850121E+01j,
        -3.620107E+01+-6.850121E+01j]


instName='SFJD_10'
print("calculating info for a "+instName)

inst1 = Response(desc=instName,units='Radians')
inst1.zeros = zer1
inst1.poles = pol1
#norm_freq=0.02
norm_freq=1.0
n1,f1 = inst1.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n1
print ('The A0 norm factor is: '+str(scale_fac)+' for f='+str(norm_freq))
#check the value
inst1.a0=1.0/n1
A01=inst1.a0
n,f = inst1.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
print('This should be close to 1: '+str(1.0/n))
h1, f1 = paz_to_freq_resp(inst1.poles, inst1.zeros, scale_fac, 0.001, 2**26, freq=True)
print(h1)

# and now for the second resp....
inst2 = Response(desc=instName,units='Radians')
inst2.zeros = zer2
inst2.poles = pol2
#norm_freq=0.02
norm_freq=1.0
n2,f2 = inst2.check_normalization(freq=norm_freq, nfft=2**24,t_sample=0.001)
scale_fac= 1.0/n2
print ('The A0 norm factor is: '+str(scale_fac)+' for f='+str(norm_freq))
#check the value
inst2.a0=1.0/n2
A02=inst2.a0
n,f = inst2.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
print('This should be close to 1: '+str(1.0/n))
h2, f2 = paz_to_freq_resp(inst2.poles, inst2.zeros, scale_fac, 0.001, 2**26, freq=True)
print(h2)

# and now for the third resp....
inst = Response(desc=instName,units='Radians')
inst.zeros = zer3
inst.poles = pol3
#norm_freq=0.02
norm_freq=1.0
n3,f3 = inst.check_normalization(freq=norm_freq, nfft=2**24,t_sample=0.001)
scale_fac= 1.0/n3
print ('The A0 norm factor is: '+str(scale_fac)+' for f='+str(norm_freq))
#check the value
inst.a0=1.0/n3
A03=inst2.a0
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
print('This should be close to 1: '+str(1.0/n))
h3, f3 = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)
print(h3)

#plotting....
plt.figure(1,figsize=(8,8))
plt.subplot(121)
plt.ioff()
#plt.loglog(f, abs(h))
print('plotting amplitude')
#plt.semilogx(1/f1,abs(h1)/A01, label='metadata resp')
#plt.semilogx(1/f2,abs(h2)/A02, label='old method resp')
#plt.semilogx(1/f3,abs(h3)/A03, label='test suite resp')
plt.semilogx(1/f1,20*np.log10(abs(h1)/A01), label='metadata resp')
plt.semilogx(1/f2,20*np.log10(abs(h2)/A02), label='CTBTO resp')
#plt.semilogx(1/f3,20*np.log10(abs(h3)/A03), label='test suite resp')
plt.legend()
plt.xlabel('Period [S]')
#plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.xlim([0.0001,10000])
#plt.xlim([40,800])
#plt.xlim([0.2,40])

plt.subplot(122)
# take negative of imaginary part
print('plotting phase')
phase = np.unwrap(np.arctan2(h1.imag, h1.real))
plt.semilogx(1/f1, phase, label='metadata resp')
phase = np.unwrap(np.arctan2(h2.imag, h2.real))
plt.semilogx(1/f2, phase, label='CTBTO resp')
#phase = np.unwrap(np.arctan2(h3.imag, h3.real))
#plt.semilogx(1/f3, phase, label='test suite resp')
plt.legend()
#plt.xlabel('Frequency [Hz]')
plt.xlabel('Period [S]')
plt.ylabel('Phase [radian]')
# title, centered above both subplots
plt.suptitle('Response of '+instName +' Seismometer')
#plt.xlim([40,800])
plt.xlim([0.0001,10000])
# make more room in between subplots for the ylabel of right plot
plt.subplots_adjust(wspace=0.3)
print('writing file')
figName=instName + 'compareResp.pdf'
plt.draw()
plt.savefig(figName,dpi=400)
#plt.show()


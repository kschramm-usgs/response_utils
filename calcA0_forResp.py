#! /bin/env python
import numpy as np
import matplotlib.pyplot as plt

from obspy.signal.invsim import paz_to_freq_resp
# this is from Austin Holland's seisutils
from response import Response

# for the STS 5/360 
#pol=[-0.01234-1j*0.01234, -0.01234+1j*0.01234, -39.18-1j*49.12, -39.18+1j*49.12]
#zer=[0 , 0]
#instName='STS-5/360'

# for the TR-360
#pol=[-0.01234+1j*0.01234, -0.01234-1j*0.01234, -196-1j*231, -196+1j*231, -732-1j*1450, -732+1j*1450, -173]
#zer=[0, 0, -108, -161]
#instName='TR-360'

#for the STS-1t5
#pol=[ -39.18-1j*49.12, -39.18+1j*49.12]
#zer=[-0.02427184, -0.02427184]
#instName='STS-1 Aaron'

#for the STS-1
##pol=[-0.01234+1j*0.01234, -0.01234-1j*0.01234, -39.18-1j*49.12, -39.18+1j*49.12]
#zer=[0, 0]
#instName='STS-1'

#for the STS-6
pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.3068, -22.3068, -501.105, -501.105, -501.105, -501.105, -501.105, -501.105, -121.58+1j*647.231, -121.58-1j*647.231, -351.858]
zer=[0, 0, -21.991, -21.991, -351.858, -605.071, -521.5+1j*960.699, -521.5-1j*960.699]
instName='STS-6'

#for the STS-6B
#pol=[-0.01234-1j*0.01234, -0.01234+1j*0.01234, -255, -15.6, -97.3+1j*401, -97.3-1j*401, -520, -375, -13300, -10500-1j*10100, -10500+1j*10100]
#zer=[0 , 0, -5910+1j*3410, -5910-1j*3410, -555, -295, -10.8, -684-1j*176, -684+1j*176 ] 
#instName='STS-6B'

# use Austin's response class to build the resp from the poles and zeros
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
norm_freq=0.02
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
scale_fac= 1.0/n
print ('The A0 norm factor is: '+str(scale_fac)+' for f='+str(norm_freq))


#poles = [-4.440 + 4.440j, -4.440 - 4.440j, -1.083 + 0.0j]
#zeros = [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j]
#scale_fac = 0.4

h, f = paz_to_freq_resp(inst.poles, inst.zeros, scale_fac, 0.001, 2**26, freq=True)

plt.figure()
plt.subplot(121)
plt.loglog(f, abs(h))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')

plt.subplot(122)
# take negative of imaginary part
phase = np.unwrap(np.arctan2(-h.imag, h.real))
plt.semilogx(f, phase)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [radian]')
# title, centered above both subplots
plt.suptitle('Frequency Response of '+inst.desc +' Seismometer')
# make more room in between subplots for the ylabel of right plot
plt.subplots_adjust(wspace=0.3)
plt.show()


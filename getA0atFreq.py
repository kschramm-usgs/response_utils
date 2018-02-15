#! /bin/env python
import numpy as np
import matplotlib.pyplot as plt

from obspy.signal.invsim import paz_2_amplitude_value_of_freq_resp, paz_to_freq_resp
# this is from Austin Holland's seisutils
from response import Response
print("script to calculate response normalization")


# for the STS 5/360 
#pol=[-0.01234-1j*0.01234, -0.01234+1j*0.01234, -39.18-1j*49.12, -39.18+1j*49.12]
#zer=[0 , 0]
#instName='STS-5/360'

# for the TR-360
#pol=[-0.01234+1j*0.01234, -0.01234-1j*0.01234, -196-1j*231, -196+1j*231, -732-1j*1450, -732+1j*1450, -173]
#zer=[0, 0, -108,-161]
#instName='TR-360'
#with the pole close to 161 removed.
#pol=[-0.01234+1j*0.01234, -0.01234-1j*0.01234, -196-1j*231, -196+1j*231, -732-1j*1450, -732+1j*1450]
## with the zero at 161 removed
#zer=[0, 0, -108]
#instName='TR-360 no -173 zero '

#for the STS-1t5
#pol=[ -39.18-1j*49.12, -39.18+1j*49.12]
#zer=[-0.02427184, -0.02427184]
#instName='STS-1 Aaron'

#for the STS-1
##pol=[-0.01234+1j*0.01234, -0.01234-1j*0.01234, -39.18-1j*49.12, -39.18+1j*49.12]
#zer=[0, 0]
#instName='STS-1'

# -6.050710e+02+ 0.000000e+00j,
#-5.215040e+02+ -9.613270e+02j,
#-5.215040e+02+  9.613270e+02j,
#for the sts 2.5 no coil
#zer=[
#  0.000000e+00+ 0.000000e+00j,
#  0.000000e+00+ 0.000000e+00j,
# -1.570800e+01+ 0.000000e+00j,
# -1.570800e+01+ 0.000000e+00j,
# -9.738940e+02+ 0.000000e+00j]
#
#pol=[
# -3.702000e-02+ 3.702000e-02j,
# -3.702000e-02+-3.702000e-02j,
# -1.604100e+01+ 0.000000e+00j,
# -1.604100e+01+ 0.000000e+00j,
# -3.273540e+02+-7.414160e+01j,
# -3.273540e+02+ 7.414160e+01j,
# -9.738940e+02+ 0.000000e+00j]
###instName='STS-2.5 nrl no coil'

#for the STS-6
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.3068, -22.3068, -501.105, -501.105, -501.105, -501.105, -501.105, -501.105, -121.58+1j*647.231, -121.58-1j*647.231, -351.858]
#zer=[0, 0, -21.991, -21.991, -351.858, -605.071, -521.5+1j*960.699, -521.5-1j*960.699]
#instName='STS-6 NRL'

#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.43, -22.43, -493.23, -501.105, -501.105, -501.105, -501.105, -501.105, -121.58+1j*647.231, -121.58-1j*647.231, -351.858]
#zer=[0, 0, -21.991, -21.991, -351.858, -605.071, -521.5+1j*960.699, -521.5-1j*960.699]
#instName='STS-6 New'

#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.3068, -22.3068, -501.105, -501.105, -501.105, -501.105, -501.105, -501.105, -121.58+1j*647.231, -121.58-1j*647.231, 351.858]
#zer=[0, 0, -21.9110,-21.9911, -351.858, -605.071, -521.5+1j*960.699, -521.5-1j*960.699]
#instName='STS-6 '

#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.3068, -22.3068, -501.105, -501.105, -501.105, -501.105, -501.105, -501.105, -121.58+1j*647.231, -121.58-1j*647.231]
## removing 2 of the 6 repeated poles:
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.3068, -22.3068, -501.105, -501.105, -501.105, -501.105, -501,105,-501.105 ]
#zer=[0, 0, -21.991, -21.991, -521.5+1j*960.699, -521.5-1j*960.699]
#zer=[0, 0, -21.9110,-21.9911, -351.858]
#instName='STS-6 no coil'
#
#for the STS-6B
#pol=[-0.01234-1j*0.01234, -0.01234+1j*0.01234, -255, -15.6, -97.3+1j*401, -97.3-1j*401, -520, -375, -13300, -10500-1j*10100, -10500+1j*10100]
#zer=[0 , 0, -5910+1j*3410, -5910-1j*3410, -555, -295, -10.8, -684-1j*176, -684+1j*176 ] 
#instName='STS-6B'

#zer=[0,0, -21.99, -21.99, -296.03, -296.03]
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206, -22.40, -22.40, -510.006, -510.006, -510.006, -510.006]
#instName='STS-6 150804'

#zer=[0,0, -21.99115, -21.99115, -605.0707450, -260.7522+960.699033j, -260.7522-960.699033j,355.0]
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206,-22.242476,-22.242476,-516.100841,-516.100841,516.100841,-516.100841,-516.100841,-516.100841,-104.3088+684.8620j,-104.3088-684.8620j,-355.0]
#zer=[0,0, -21.99115, -21.99115, 355.0]
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206,-22.242476,-22.242476,-516.100841,-516.100841,516.100841,-516.100841,-516.100841,-516.100841,-104.3088+684.8620j,-104.3088-684.8620j,-355.0]
#zer=[0,0, -21.99115, -21.99115, -605.0707450, -260.7522+960.699033j, -260.7522-960.699033j,355.0]
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206,-22.242476,-22.242476,-516.100841,-516.100841,-516.100841,-516.100841,-104.3088+684.8620j,-104.3088-684.8620j,-355.0]
#zer=[0,0, -21.99115, -21.99115, -260.7522+960.699033j, -260.7522-960.699033j,355.0]
#pol=[-0.0123+1j*0.01206, -0.0123-1j*0.01206,-22.242476,-22.242476,-516.100841,-516.100841,-516.100841,-516.100841,-104.3088+684.8620j,-104.3088-684.8620j,-355.0]
#instName='STS-6 176211'



# for a TR-120
#instName='LB Trillium 120'
#zer=[
#    +0.000000e+00    +0.000000e+00j,
#    +0.000000e+00    +0.000000e+00j,
#    +0.000000e+00    +0.000000e+00j,
#    -1.060000e+02    +0.000000e+00j,
#    -1.580000e+02    +0.000000e+00j
#]
#pol=[
#    -3.859000e-02    +3.649000e-02j,
#    -3.859000e-02    -3.649000e-02j,
#    -1.580000e+02    +1.930000e+02j,
#    -1.580000e+02    -1.930000e+02j,
#    -6.390000e+02    +1.418000e+03j,
##    -6.390000e+02    -1.418000e+03j,
#    -1.900000e+02    +0.000000e+00j
#]
#A0 = 8.531186e+17

#for TR-120 g2 no cal coil
##zer = [
# 0.000000E+00+ 0.000000E+00j,
# 0.000000E+00+ 0.000000E+00j,
#-3.163000E+01+ 0.000000E+00j,
#-3.500000E+02+ 0.000000E+00j]
#pol=[
#-3.661400E-02+ 3.705900E-02j,
#-3.661400E-02+-3.705900E-02j,
#-3.255000E+01+ 0.000000E+00j,
#-1.590000E+02+ 0.000000E+00j,
#-3.640000E+02+ 4.040000E+02j,
#-3.640000E+02+-4.040000E+02j,
#-1.260000E+03+ 0.000000E+00j,
#-4.900000E+03+ 5.200000E+03j,
#-4.900000E+03+-5.200000E+03j,
#-7.100000E+03+ 1.700000E+03j,
#-7.100000E+03+-1.700000E+03j]
#instName='TR120 2g no coil'
# for a T120
#instName='T120P'
#zer=[
#     0.000000e+00+ 0.000000e+00j,
#     0.000000e+00+ 0.000000e+00j,
#    -9.000000e+01+ 0.000000e+00j,
#    -1.607000e+02+ 0.000000e+00j,
#    -3.108000e+03+ 0.000000e+00j]
#pol=[
#    -3.852000e-02+ 3.658000e-02j,
#    -3.852000e-02+-3.658000e-02j,
#    -1.780000e+02+ 0.000000e+00j,
#    -1.350000e+02+ 1.600000e+02j,
#    -1.350000e+02+-1.600000e+02j,
#    -6.710000e+02+ 1.154000e+03j,
#    -6.710000e+02+-1.154000e+03j]

#gs-13
#instName='GS-13'
#zer=[0.00+0.000j, 0.000-0.000j]
##pol=[-4.443000+4.443000j,-4.443000-4.44300j]
#pol=[-4.449800+4.628900j,-4.449800-4.62890j]

#resp info fit to cals on sts-6 (WVT/00)
#instName='WVT STS-6'
#zer = [
# 0.000000E+00+ 0.000000E+00j,
# 0.000000E+00+ 0.000000E+00j,
#-0.983300E+01+ 0.000000E+00j,
#-0.983300E+01+ 0.000000E+00j,
#-2.832250E+02+ 0.000000E+00j,
#-6.050710E+02+ 0.000000E+00j,
#-5.215040E+02+-9.606990E+02j,
#-5.215040E+02+ 9.606990E+02j]
#pol= [
#-1.253000E-02+ 1.254000E-02j,
#-1.253000E-02+-1.254000E-02j,
#-0.990855E+01+ 0.000000E+00j,
#-0.990855E+01+ 0.000000E+00j,
#-6.086859E+02+ 0.000000E+00j,
#-6.086859E+02+ 0.000000E+00j,
#-6.086859E+02+ 0.000000E+00j,
#-6.086859E+02+ 0.000000E+00j,
#-6.086859E+02+ 0.000000E+00j,
#-6.086859E+02+ 0.000000E+00j,
#-1.215800E+02+ 6.472310E+02j,
#-1.215800E+02+-6.472310E+02j,
#-3.599300E+02+ 0.000000E+00j]
#

#instName='SFJD 10 for ctbto'
#zer=[
#    0.000000E+00+ 0.000000E+00j,
#    0.000000E+00+ 0.000000E+00j,
#   -9.424780E+00+ 0.000000E+00j,
#   -6.283190E+02+ 0.000000E+00j,
#   -5.654870E+02+ 9.794520E+02j,
#   -5.654870E+02+-9.794520E+02j]
#pol=[
#   -3.731580E-02+-3.670000E-02j,
#   -3.731580E-02+ 3.670000E-02j,
#   -9.738940E+00+ 0.000000E+00j,
#   -2.199110E+02+ 1.382300E+02j,
#   -2.199110E+02+-1.382300E+02j,
#   -2.199110E+02+ 6.848670E+02j,
#   -2.199110E+02+-6.848670E+02j]

instName='ANMO 10 for CTBTO cals'
zer=[
  0.000000E+00+ 0.000000E+00j,
  0.000000E+00+ 0.000000E+00j,
 -3.163000E+01+ 0.000000E+00j,
 -1.600000E+02+ 0.000000E+00j,
 -3.500000E+02+ 0.000000E+00j,
 -3.177000E+03+ 0.000000E+00j]
pol=[
 -3.680560E-02+ 3.626490E-02j,
 -3.680560E-02+-3.626490E-02j,
 -3.255000E+01+ 0.000000E+00j,
 -1.420000E+02+ 0.000000E+00j,
 -3.640000E+02+ 4.040000E+02j,
 -3.640000E+02+-4.040000E+02j,
 -1.260000E+03+ 0.000000E+00j,
 -4.900000E+03+-5.200000E+03j,
 -4.900000E+03+ 5.200000E+03j,
 -7.100000E+03+-1.700000E+03j,
 -7.100000E+03+ 1.700000E+03j]





print("calculating info for a "+instName)
# use Austin's response class to build the resp from the poles and zeros
inst = Response(desc=instName,units='Radians')
inst.zeros = zer
inst.poles = pol
#norm_freq=0.05
#norm_freq=0.02
norm_freq=1.0
n,f = inst.check_normalization(freq=norm_freq, nfft=2**26,t_sample=0.001)
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
plt.grid(which='both',linestyle=':')
#plt.xlim([0.2,40])


plt.subplot(122)
# take negative of imaginary part
phase = np.unwrap(np.arctan2(h.imag, h.real))
plt.semilogx(f, phase)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [radian]')
# title, centered above both subplots
plt.suptitle('Frequency Response of '+inst.desc +' Seismometer')
#plt.xlim([0.0001,10000])
# make more room in between subplots for the ylabel of right plot
plt.subplots_adjust(wspace=0.3)
plt.grid(which='both',linestyle=':')

plt.show()


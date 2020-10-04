import SEOBNRE
from pycbc.waveform import get_td_waveform
import numpy as np
import lal
import lalsimulation as lalsim
import pylab


longAscNodes = 0
eccentricity = 0 
meanPerAno = 0
approx=lalsim.SEOBNRv1
nonGRdict = lal.CreateDict()
m1 = 30
m2 = 30
s1 = [0,0,0]#[0.4,-0.2,0.43]
s2 = [0,0,0]#[-0.1,0.8,0]
dist = 400.
#iota = np.pi*0.4
iota = 0.
phiRef = 0
deltaT = 1./4096/4
f_ref = 0
f_low = 30


for apx in ['SEOBNRE']:
    hp, hc = get_td_waveform(approximant=apx,
                                 mass1=m1,
                                 mass2=m2,
                                 delta_t=deltaT,
                                 eccentricity = 0,
                                 distance = dist,
                                 coa_phase = phiRef,
                                 f_lower=f_low,
                                 longAscNodes = longAscNodes)

    pylab.plot(hp.sample_times, hp, label='hplus:'+apx)



hpt1, hct1 = lalsim.SimInspiralChooseTDWaveform(m1 * lal.MSUN_SI, m2 * lal.MSUN_SI,\
                                                          s1[0], s1[1], s1[2],\
                                                          s2[0], s2[1], s2[2],\
                                                          dist * 1e6 * lal.PC_SI, iota, phiRef,\
                                                          longAscNodes, eccentricity, meanPerAno,
                                                          deltaT, f_low, f_ref,\
                                                          nonGRdict, approx)

t1 = np.arange(hpt1.data.length, dtype=float) * hpt1.deltaT
t1 = t1 + hpt1.epoch
pylab.plot(t1, hpt1.data.data, label='hplus:SEOBNRv4,lal')

pylab.ylabel('Strain')
pylab.xlabel('Time (s)')
pylab.legend()
pylab.show()
import SEOBNRE
from pycbc.waveform import get_td_waveform
import numpy as np
import lal
import lalsimulation as lalsim
import pylab


longAscNodes = 0
eccentricity = 0 
meanPerAno = 0
approx=lalsim.SEOBNRv4
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
f_low = 40


for coa in [1,2,3]:
    hp, hc = get_td_waveform(approximant='SEOBNRv4',
                                 mass1=m1,
                                 mass2=m2,
                                 delta_t=deltaT,
                                 eccentricity = 0,
                                 distance = dist,
                                 coa_phase = 0,
                                 f_lower=f_low,
                                 long_asc_nodes = coa)

    pylab.plot(hp.sample_times, hp, label='hplus:longAscNodes='+str(coa))

pylab.ylabel('Strain')
pylab.xlabel('Time (s)')
pylab.legend()
pylab.show()
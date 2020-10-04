from pycbc.waveform import td_approximants, fd_approximants
import SEOBNRE
# List of td approximants that are available


import pylab
from pycbc.waveform import get_td_waveform

for apx in ['SEOBNRE']:
    hp, hc = get_td_waveform(approximant=apx,
                                 mass1=30,
                                 mass2=30,
                                 delta_t=1.0/4096,
                                 eccentricity = 0,
                                 distance = 400,
                                 f_lower=40,
                                 coa_phase = 0,
                                 longAscNodes = 0)

    pylab.plot(hp.sample_times, hp, label=apx)

import numpy as np 
comp = np.loadtxt("./comp.txt")

pylab.plot(comp[:,0],comp[:,1],label='comp')
pylab.ylabel('Strain')
pylab.xlabel('Time (s)')
pylab.legend()
pylab.show()

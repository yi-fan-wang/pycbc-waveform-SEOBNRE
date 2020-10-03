from pycbc.waveform import td_approximants, fd_approximants
import SEOBNRE
# List of td approximants that are available


import pylab
from pycbc.waveform import get_td_waveform

for apx in ['SEOBNRv4','SEOBNRE']:
    hp, hc = get_td_waveform(approximant=apx,
                                 mass1=30,
                                 mass2=30,
                                 delta_t=1.0/4096,
                                 eccentricity = 0,
                                 distance = 400,
                                 coa_phase = 0,
                                 f_lower=40)

    pylab.plot(hc.sample_times, hc, label='hplus:'+apx)

pylab.ylabel('Strain')
pylab.xlabel('Time (s)')
pylab.legend()
pylab.show()
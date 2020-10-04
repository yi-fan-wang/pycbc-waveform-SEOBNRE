from pycbc.waveform import td_approximants, fd_approximants
import SEOBNRE
# List of td approximants that are available


import pylab
from pycbc.waveform import get_td_waveform

for apx in ['SEOBNRv4']:
  for coa in [0,1]:
    hp, hc = get_td_waveform(approximant=apx,
                                 mass1=30,
                                 mass2=30,
                                 delta_t=1.0/4096,
                                 eccentricity = 0,
                                 distance = 400,
                                 coa_phase = coa,
                                 f_lower=20)

    pylab.plot(hp.sample_times, hp, label='hplus:'+apx+'coa_phase'+str(coa))

pylab.ylabel('Strain')
pylab.xlabel('Time (s)')
pylab.legend()
pylab.show()
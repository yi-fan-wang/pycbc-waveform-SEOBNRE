import pylab
from pycbc.waveform import get_td_waveform

for apx in ['SEOBNRv2', 'SEOBNRv4']:
    hp, hc = get_td_waveform(approximant=apx,
                                 mass1=30,
                                 mass2=30,
                                 spin1z=0,
                                 delta_t=1.0/4096,
                                 f_lower=20)

    pylab.plot(hp.sample_times, hp, label=apx)

pylab.ylabel('Strain')
pylab.xlabel('Time (s)')
pylab.legend()
pylab.show()
from ctypes import *
import numpy
from pycbc.types import TimeSeries,FrequencySeries
import os
def SEOBNRE_td(**kwargs):

    lib = CDLL('$LD_LIBRARY_PATH/libSEOBNRE.so')

    waveform_generation = lib.genwaveform
    waveform_generation.argtypes = [POINTER(c_double), POINTER(c_double),\
    POINTER(c_double),POINTER(c_int),\
    c_double, c_double,c_double, c_double,c_double, c_double,c_double, c_double,c_double, c_double, c_double]

    params = {'coa_phase': 0.0,'delta_t': None,'mass1': None,'mass2': None,\
    'spin1z': 0.0,'spin2z': 0.0,'f_lower': None, 'eccentricity': 0.0,'distance': 1.0,'inclination': 0.0,\
    'long_asc_nodes': 0.0}
    for value in params:
        if value in kwargs:
            params[value] = kwargs[value]

    # TODO: avoid this ugly method
    # Initilization
    arraysize = 100000

    hplus = numpy.zeros(arraysize)
    hcross = numpy.zeros(arraysize)
    truesize = numpy.array([0])
    t0 = numpy.array([1.])

    # Generate the waveform
    waveform_generation(hplus.ctypes.data_as(POINTER(c_double)),hcross.ctypes.data_as(POINTER(c_double)),
        t0.ctypes.data_as(POINTER(c_double)),truesize.ctypes.data_as(POINTER(c_int)),
        params['coa_phase'],params['delta_t'],params['mass1'],params['mass2'],params['spin1z'],
        params['spin2z'],params['f_lower'],params['eccentricity'],params['distance'],params['inclination'],params['long_asc_nodes'])

    hplus = hplus[:truesize[0]]
    hcross = hcross[:truesize[0]]
    time = numpy.arange(t0[0],t0[0]+params['delta_t']*truesize[0], params['delta_t'])

    # Build the TimeSeries format

    hp = TimeSeries(hplus[:], dtype=hplus.dtype,delta_t=params['delta_t'], epoch=time[0])
    hc = TimeSeries(hcross[:],dtype=hcross.dtype,delta_t=params['delta_t'], epoch=time[0])

    return hp,hc

def SEOBNRE_fd(**kwargs):
    hp, hc = SEOBNRE_td(**kwargs)

    return (hp.to_frequencyseries(),
            hc.to_frequencyseries())


import numpy

def mkricker(amplitude, width, t0):

    def ricker(t):
        tp = t - t0
        s = (numpy.pi*width*tp)**2
        return amplitude * (1.0 - 2.0*s) * numpy.exp(-s)

    return ricker
        

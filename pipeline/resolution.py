import numpy as np


def window(k,p,R):
    """
    Deconvolution kernel.
    k: velocity wavenumber
    p: pixel width (in velocity, km/s)
    R: resolution
    """
    resolutionCorrection= np.exp(-0.5*k**2*R**2)
    tophat =  np.sinc(0.5*k*p/np.pi) # THIS PI NEEDS TO BE HERE#np.sin(k*p/2)/(k*p/2)
    deconvolution_kernel = resolutionCorrection*tophat
    return deconvolution_kernel

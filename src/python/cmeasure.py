import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.fftpack as sfft
import scipy.linalg as slag
import json

import genscript.parameter as par
import genscript.myarray as mar
import genscript.mpiutil as mpi

import misc.helper as helper
import covmat as covm

import cyth.covm as cyth_cov








def autocorr(d, auto_type='FFT', zero_padding=True):
    ''' ->> conduct the auto-correlation of the map <<- ''' 

    if auto_type=='scipy':
        cor=sig.correlate2d(d, d) 
	print 'cor shape:', cor.shape

    if auto_type=='FFT':
        # ->> zero padding <<- #
        if zero_padding==True:
            dx, dy=d.shape
            dd=np.zeros((dx*2, dy*2))
            dd[:dx,:dy]=d
        else:
            dd=d

        dk=np.fft.rfft2(dd)
        cor=np.fft.irfft2(dk*np.conjugate(dk))
        #cor=np.fft.fftshift(cor)

    return cor[:d.shape[0], :d.shape[1]]



def pk_fft_2d(d, dmap_res, zero_padding=False):
    ''' ->> measure 2D power spectrum with FFT <<- '''

    if zero_padding==True:
        # ->> zero padding <<- #
        dx, dy=d.shape
        dd=np.zeros((dx*2, dy*2))
        dd[:dx,:dy]=d
    else:
        dd=d
    
    dk=np.fft.fft2(dd)
    pk=np.fft.fftshift(np.absolute(dk)**2.)

    #return pk*np.prod(dmap_res)/np.prod(d.shape)*(2.*np.pi)**2.
    return pk*np.prod(dmap_res)**2./(2.*np.pi)**2.





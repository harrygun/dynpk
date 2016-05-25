import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.fftpack as sfft

import misc.helper as helper







def autocorr(d, auto_type='FFT'):
    ''' ->> conduct the auto-correlation of the map <<- ''' 

    if auto_type=='scipy':
        cor=sig.correlate2d(d, d) 

    if auto_type=='FFT':
       dk=np.fft.rfft2(d)
       #cor=(dk*np.conjugate(dk).astype(np.float64)
       cor=np.fft.irfft2(dk*np.conjugate(dk))

    #print 'auto correlation:', cor.shape, type(cor[0,0])
    print cor.max(), cor.min()


    return cor





def quad_pk(d, pk_fid):
    ''' ->> construct fiducial estimator, given the fiducial pk <<- '''

    


    return

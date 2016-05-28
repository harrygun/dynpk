import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.fftpack as sfft

import genscript.parameter as par

import misc.helper as helper





''' ------------------------------------------------------------ 
              ->>  Quadratic Estimator <<-
    ------------------------------------------------------------
'''
defaultQuadEestValueDict={
    'map_dimension':   2100,
    }


defaultQuadEestNameDict={
    'map_dimension':   'm_dim',
    }



class QuadEstimator(par.Parameters):


    def __init__(self, paramfname=None, section=None, def_var=True, 
                 prog_control=None, **pardict):

        return




def wpixel_fourier(kl, x_a, dx, pixel_type='rect'):
    # ->> 
    w=np.sinc(kl)

    return


def quad_init(d):
    # ->> initialization of quadratic estimator <<- #
    # ->> 

    return



def quad_pk(epar, d, fcm):
    ''' ->> construct fiducial estimator, given the fiducial pk <<- 
            epar:  parameters for quadratic estimator 
            fcm: fidicual covariance matrix  
    '''
    # ->>  




    return







''' ->> other routines <<- '''
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





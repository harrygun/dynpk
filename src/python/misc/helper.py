''' ->> some routines to help to start, like random mask generator <<- '''
import numpy as np
import scipy.fftpack as sfft
import pylab as pl
import random as rnd




rbinary_list=lambda n: [rnd.randint(0,1) for b in range(1,n+1)]


def bin_mask(size, mask_type='default'):
    ''' ->> random generator of binary mask <<- '''

    if mask_type=='default':
        m=np.array(rbinary_list(np.prod(size))).reshape(size)
    else:
        raise Exception

    return m




def gen_mask(mapsize, npt=10, mask_type='default'):
    ''' ->> random mask generator <<- '''

    s=np.min(mapsize)
    m=np.ones(mapsize)
     
    #->> 
    x, y=np.mgrid[:mapsize[0],:mapsize[1]]
    rmin, rmax = s*0.05, s*0.3

    rl=rmin+np.random.rand(npt)*(rmax-rmin)
    cl_x=mapsize[0]*np.random.rand(npt)
    cl_y=mapsize[1]*np.random.rand(npt)

    #print rl, cl_x, cl_y

    for i in range(npt):
        m*=1.-np.exp(-((x-cl_x[i])**2.+(y-cl_y[i])**2.)/rl[i]**2.)

    if mask_type=='binary':
         _frac_=0.3
         m[np.where(m<_frac_)]=0
         m[np.where(m>=_frac_)]=1

    return m



'''-----------------------------------------------------------
       ->>     simple estimation of the map     <<- 
   -----------------------------------------------------------
'''

def map_stat_est(d):
    # ->> some simple estimation of the map <<- #
    return np.array([np.mean(d), np.std(d)])



def random_noise_generator(d, snr):
    ''' ->> given the original data map <<- 
        ->> snr:   S/N ratio
    '''
    #->> statistics
    mean, std=map_stat_est(d)

    # ->> 
    nvar=std/snr




    return










''' ------------------------------------------------------------ 
             ->>       Some other routines          <<-
    ------------------------------------------------------------
'''

def klist_fft(real_size, kdim):

    k_min=np.ones(real_size.shape)*2.*np.pi/np.array(real_size)

    print 'klist_fft:', real_size.shape, kdim, len(kdim), k_min

    return np.array([sfft.fftfreq(kdim[i], d=1./float(kdim[i]))*k_min[i]\
                    for i in range(len(kdim)) ] )

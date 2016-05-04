''' ->> some routines to help to start, like random mask generator <<- '''
import numpy as np
import pylab as pl
import random as rnd




rbinary_list=lambda n: [rnd.randint(0,1) for b in range(1,n+1)]



def gen_mask(mapsize, pix_frac=10, mask_type='default'):
    '''->> random mask generator <<- '''

    crit_pix=5
    mapsize=np.array(mapsize)
    pixsize=np.ceil(mapsize/pix_frac).astype(int)


    if mask_type=='default':
        mask=np.array(rbinary_list(np.prod(pixsize))).reshape(pixsize)
    else:
        raise Exception

    #->> 
    m=np.zeros(mapsize) 
    m[np.where(mask )] = 1

    return m



def covm_sig_guess(d, method='default'):
    ''' ->> Guess the signal covariance matrix <<- '''

    if method=='default':
        # ->>
        pass


    return

''' ->> some routines to help to start, like random mask generator <<- '''
import numpy as np
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




def gen_mask(mapsize, mask_type='default'):
    ''' ->> random mask generator <<- '''
     
    nn=20

    for i in range(nn):
        #->> 
         


    return m




def covm_sig_guess(d, method='default'):
    ''' ->> Guess the signal covariance matrix <<- '''

    if method=='default':
        # ->>
        pass


    return

''' ->> some routines to help to start, like random mask generator <<- '''
import numpy as np
import pylab as pl
import random as rnd




rbinary_list=lambda n: [rnd.randint(0,1) for b in range(1,n+1)]




def gen_mask(size):
    '''->> random mask generator <<- '''

    m=rbinary_list(np.prod(size)).reshape(size)

    return m



def covm_sig_guess(d, method='default'):
    ''' ->> Guess the signal covariance matrix <<- '''

    if method=='default':
        # ->>
        pass


    return

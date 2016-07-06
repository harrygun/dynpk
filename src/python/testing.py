import os
import pylab as pl
import numpy as np
import matplotlib.colors as colors

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import genscript.myarray as mar

import misc.helper as helper
import quadestimator as qde
import cmeasure as cms





def do_test(p, d):

    if p.do_testing=='False':
        return

    # ->> 
    dmean=d.mean()
    print 'mean: ', dmean
    d=(d-dmean)[:200,:200]


    #test_type='fft'  #'mask'
    test_type='auto_corr'  

    if test_type=='mask':
        #size=(200,400)
        size=d.shape  #(400,900)

        mask_type='default'  #'binary'
        #mask_type='binary'

        m=helper.gen_mask(size, npt=8, mask_type=mask_type)

        cb=pl.imshow(m*d)
        pl.colorbar(cb)
        pl.show()

    if test_type=='fft':
        #->> 
        #dk= 
	pass

    if test_type=='auto_corr':
        cor=qde.autocorr(d, auto_type='FFT')
        #cor=qde.autocorr(d, auto_type='scipy')

	#cb=pl.imshow(cor, norm=colors.LogNorm() )
	#cb=pl.imshow(cor, vmin=-50, vmax=50)

	cb=pl.imshow(cor)
	pl.show()



    p.finalize()
    quit()



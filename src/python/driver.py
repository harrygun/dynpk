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







def do_test(p, d):

    if p.do_testing=='False':
        return

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
        cor=qde.autocorr(d)

	#cb=pl.imshow(cor, norm=colors.LogNorm() )
	cb=pl.imshow(cor, vmin=-50, vmax=50)
	#pl.hist(cor.flatten(), (-50, 50)) 

	pl.show()



    p.finalize()
    quit()






param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'a_init': 1e-2,
    'smooth_R': 0,
    'smooth_type': 'Gauss', 
    'smooth_R_list_type':  'linear'
    }

prog_control={
    'do_testing': True, 
    #-------------------------------#
    #-------------------------------#
    }





if __name__=='__main__':

    # ->> initialization <<- #
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(**init_dict)

    root='../../workspace/'


    # ->> data importing <<- #
    fn=root+'data/DynSpec/Rickett_53560dspec.npy'
    d=np.load(fn)
    print 'data shape:', d.shape


    #->> do some testing <<- #
    dmean=d.mean()
    dd=d-dmean
    print 'mean: ', dmean

    do_test(p, dd)


    # ->> construction KL modes <<- #
    # ->> (1). guess the covariance matrix <<- # 
    mask=helper.gen_mask(d.shape, npt=10, mask_type='default')


    # ->> (2). construct <<- # 
     




    #->> write files <<- #





    
    # ->> The End <<- #
    p.finalize()




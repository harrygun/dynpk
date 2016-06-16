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






param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'cosmology_parameter_fname': 'parameters/cosparameter.cfg',
    'cosmology_parameter_sec': 'Cosmology_Parameters',
    'qestmator_sec':       'Quadratic_Estimator',
    }

prog_control={
    'do_testing': False, 
    #-------------------------------#
    #-------------------------------#
    }





if __name__=='__main__':

    # ->> initialization <<- #
    sec='General'
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(section=sec, **init_dict)

    root='../../workspace/'

    # ->> data importing <<- #
    #fn=root+'data/DynSpec/Rickett_53560dspec.npy'
    fn=root+'data/sims/ds1.npy'

    d=np.load(fn)
    #print 'data shape:', d.shape

    #->> do some testing <<- #
    do_test(p, d)


    ''' ->> now start to work <<- '''

    #->> define class <<- #
    qe_dict={}
    dmap=d[:200,:200]

    #->>parafname='same as parameter file' 
    qe=qde.QuadestPara(paramfname=p.paramfname, section=p.qestmator_sec,
            prog_control=p, dmap=dmap, **qe_dict)

    print '\n->> QuadestPara parameters:\n', qe.paramdict
    print qe.plist.shape

    # ->> construction KL modes <<- #
    # ->> (1). guess the covariance matrix <<- # 
    # ->> mask=helper.gen_mask(d.shape, npt=10, mask_type='default')


    # ->> (2). construct the quadratic estimator <<- # 
    # a). 2D axis list 
    n_it=0
    qe.quadest_iteration(self, pk_fid, n_it)


    # b).  

    # ->> initialization <<- #

    # 


    #->> write files <<- #





    
    # ->> The End <<- #
    p.finalize()




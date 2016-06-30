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

    print 'mpi.rank=', mpi.rank


    '''----------------------------------------------
                ->>    now we start    <<- 
       ----------------------------------------------'''
    #->> data initialization <<- #
    qe_dict={'calculate_dcov': False, 
             'fname_dcov':     root+'result/dcov.npz',
	    }
    dmap=d[:100,:100]

    ''' #->> Initialzing quadratic estimator class <<- #
        #->> parafname='same as parameter file' 
        #->> calculating dcov
    '''
    qe=qde.QuadestPara(paramfname=p.paramfname, section=p.qestmator_sec,
            prog_control=p, dmap=dmap, **qe_dict)

    print '\n->> QuadestPara parameters:\n', qe.paramdict


    # ->> start estimator evaluation, n_it=0 without iternation <<- #
    n_it=0

    # ->> set fiducial power spectrum <<- #
    pk_fid=qe.fid_pk_first_guess(guess_option='simplest')
    print '->> fiducial pk initialization done.'

    # ->> initialize noise covmat <<- #
    noise_level='noiseless'
    qe.covn_vec_init(noise_level=noise_level)
    print '->> noise level initialized.'

    # ->> run estimator <<- #
    print '->> start the quadratic estimator.'
    Qi=qe.quadest_iteration(pk_fid, n_it)



    #->> write files <<- #
    fn_qi=root+'result/Qi.npz'

    if mpi.rank0:
        np.savez(fn_qi, Qi=Qi)
        mpi.barrier()
     


    
    # ->> The End <<- #
    p.finalize()




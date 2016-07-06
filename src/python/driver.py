import os
import pylab as pl
import numpy as np
import scipy as sp
import matplotlib.colors as colors

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import genscript.myarray as mar
import genscript.myplot as mpl

import misc.helper as helper
import quadestimator as qde
import cmeasure as cms

import testing as test








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
    fn=root+'data/sims/dsU1.npy'

    d=np.load(fn)

    #->> do some testing <<- #
    test.do_test(p, d)

    print 'mpi.rank=', mpi.rank


    '''----------------------------------------------
                ->>    now we start    <<- 
       ----------------------------------------------'''

    # ->> lower data resolution <<- #
    zoom_factor=0.5
    _dmap_=d[:100,:100]-np.mean(d[:100,:100])
    dmap=sp.ndimage.interpolation.zoom(_dmap_, zoom_factor)



    #->> data initialization <<- #
    qe_dict={'calculate_dcov':   False, 
             'fname_dcov':       root+'result/dcov_fft_50x50.npz',
	     'map_zoom_factor':  zoom_factor,
	     'get_bp_type':      'FFT',
	    }

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
    pk_fid=qe.fid_pk_first_guess()
    print '->> fiducial pk initialization done.'

    # ->> initialize noise covmat <<- #
    noise_level='noiseless'
    qe.covn_vec_init(noise_level=noise_level)
    print '->> noise level initialized.'

    # ->> run estimator <<- #
    print '->> start the quadratic estimator.'
    Qi, Fij=qe.quadest_iteration(pk_fid, n_it)

    print mpi.rank, '->> quadratic estimator done.'

    #->> write files <<- #
    fn_qi=root+'result/Qi_fft.npz'

    mpi.barrier()
    if mpi.rank0:
        np.savez(fn_qi, Qi=Qi, Fij=Fij)
     

    # ->> The End <<- #
    p.finalize()
    mpi.finalize()




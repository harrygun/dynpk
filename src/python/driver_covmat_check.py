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





def auto_corr_test(p, d):

    if p.do_testing=='False':
        return

    cor=qde.autocorr(d, auto_type='FFT')
    cb=pl.imshow(cor)
    pl.show()

    quit()





param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'cosmology_parameter_fname': 'parameters/cosparameter.cfg',
    'cosmology_parameter_sec': 'Cosmology_Parameters',
    'qestmator_sec':       'Quadratic_Estimator',
    }

prog_control={
    'do_testing': True, 
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
    fn=root+'data/sims/dsU1.npy'

    d=np.load(fn)
    #print 'data shape:', d.shape

    #->> do some testing <<- #
    auto_corr_test(p, d)

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







    #->> write files <<- #
    fn_qi=root+'result/Qi.npz'

    if mpi.rank0:
        np.savez(fn_qi, Qi=Qi)
        mpi.barrier()
     


    
    # ->> The End <<- #
    p.finalize()




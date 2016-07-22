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


import cyth.covm as cyth_cov



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
    fn_in=root+'result/dcov_out.dat'

    if mpi.rank0:
        data=np.fromfile(fn_in).reshape(2500,50,50)
    else: data=0

    dcov=mpi.bcast(data, root=0)

    # ->> 
    npt=2500
    mdim_t, mdim_f = 50, 50
    dcov_full=np.zeros((2500,2500))

    mpirange=mpi.mpirange(2500)
    for i in mpirange:

        fn_out=root+'result/dcov_full_fft/dcov_full_'+str(i)+'.dat'
	print 'rank=', mpi.rank, fn_out

        cyth_cov.get_full_dcov(dcov, dcov_full, i, npt, mdim_t, mdim_f)
        dcov_full.tofile(fn_out)


    # ->> The End <<- #
    p.finalize()
    mpi.finalize()




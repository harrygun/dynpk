''' ->> for generating the 1D/2D data sending to Carly <<- '''
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
import cyth.covm as cyth_cov






def fisher(dcov, cov):




    return



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
    fn=root+'data/sims/dsU1.npy'
    #fn=root+'data/sims/ds1.npy'

    d=np.load(fn)
    #print 'data shape:', d.shape

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
             'fname_dcov':       root+'result/dcov.npz',
	     'map_zoom_factor':  zoom_factor,
	     'get_bp_type':      'FFT',
	    }


    ''' #->> Initialzing quadratic estimator class <<- #
        #->> parafname='same as parameter file' 
        #->> calculating dcov
    '''
    skip_init=True  # ->> skip initialize bandpowre and dcov <<- #
    qe=qde.QuadestPara(paramfname=p.paramfname, section=p.qestmator_sec,
            prog_control=p, dmap=dmap, skip_init=True, **qe_dict)

    print '\n->> QuadestPara parameters:\n', qe.paramdict

    # ->> initialize FFT band power <<- #
    bp_dict={'dmap_shape':   qe.dmap.shape, } 
    qe.band_power_init(**bp_dict)


    # ->> initialize dcov <<- #
    fname_dcov_fft=root+'result/dcov_fft_50x50.npz'
    qe.dcov_init(fname_dcov_fft)

    # ->> 
    print 'klist', qe.klist.shape, qe.kt_list.shape, qe.kf_list.shape
    print qe.kf_list[25]

    print 'index where kf=0:', np.where(qe.klist[1]==0)
    print 'or in 2D:', np.where(qe.klist.reshape(2,50,50)[1]==0)


    print 'Let me just select k_f=0, and output dcov, cov, and icov.'
    print 'Also, for dcov in 1D, I have to set delta_f=0 among all pixels.'

    print 'dcov:', qe.dcov.shape


    #
    fno_dcov=root+'result/1d/dcov.dat'
    fno_cov=root+'result/1d/cov.dat'
    ddcov=qe.dcov.reshape(50,50,100,100)[:,25,:50,0]


    ff=np.fromfile(root+'result/cov.dat').reshape(50,50,50,50)
    ccov=ff[:,0,:,0]



    # ->> output <<- #
    ddcov.tofile(fno_dcov)
    ccov.tofile(fno_cov)
     








    # ->> The End <<- #
    p.finalize()




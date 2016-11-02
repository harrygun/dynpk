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
import quadestimator_r1D as qde1d
import cmeasure as cms
import cyth.covm as cyth_cov









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

    print 'mpi.rank=', mpi.rank

    # ->> data importing <<- #
    root='../../workspace/'

    fn=root+'data/sims/dsU1.npy'
    d=np.load(fn)



    '''----------------------------------------------
                ->>    now we start    <<- 
       ----------------------------------------------'''
    # ->> lower data resolution <<- #
    #zoom_factor=0.5
    zoom_factor=1  #0.5
    #n_pts=500
    n_pts=50
    ncol=10
    _dmap_=d[:n_pts,ncol]-np.mean(d[:n_pts,ncol])
    dmap=sp.ndimage.interpolation.zoom(_dmap_, zoom_factor)



    #->> 
    folder=root+'result/cr1d/npt_'+str(n_pts)
    print 'folder existence:', os.path.isdir(folder)

    if not os.path.isdir(folder):
        cmd="mkdir "+folder
	os.system(cmd)

    dout_name=folder+'/dmap.dat'
    pk_name=folder+'/bpk.dat'
    dcov_name=folder+'/dcov_py.dat'






    #->> data initialization <<- #
    qe_dict={'calculate_dcov':   True,
             'fname_dcov':       root+'result/cr1d/npt_50/dcov_r1d_fft_24.npz',
	     'map_zoom_factor':  zoom_factor,
	     'get_bp_type':      'FFT',
	    }


    ''' #->> Initialzing quadratic estimator class <<- #
        #->> parafname='same as parameter file' 
        #->> calculating dcov
    '''
    skip_init=False # ->> skip initialize bandpowre and dcov <<- #
    qe=qde1d.Quadest1dPara(paramfname=p.paramfname, section=p.qestmator_sec,
            prog_control=p, dmap=dmap, skip_init=skip_init, **qe_dict)

    print '\n->> QuadestPara parameters:\n', qe.paramdict


    # ->> 
    # ->> start estimator evaluation, n_it=0 without iternation <<- #
    n_it=0

    # ->> set fiducial power spectrum <<- #
    _pk_=qe.fid_pk_first_guess()
    kdim=len(_pk_)
    pk_fid=_pk_[kdim/2+1:]


    print 'n_bp:', len(pk_fid), qe.dcov.shape, pk_fid.shape, dmap.shape
    #print qe.klist, qe.klist_low, qe.klist_up

    bpk=np.concatenate((qe.klist, qe.klist_low, qe.klist_up, pk_fid), axis=0)

    print 'output bpk.shape:', bpk.shape, pk_fid.shape, 
    print kdim, qe.klist[kdim/2+1:].shape, qe.klist.shape, qe.klist_low.shape, \
          qe.klist_up.shape

    #->> write files <<- #
    bpk.tofile(pk_name)
    dmap.tofile(dout_name)


    #print 'dcov:\n', qe.dcov[124]
    qe.dcov.tofile(dcov_name)

    
    # ->> The End <<- #
    p.finalize()




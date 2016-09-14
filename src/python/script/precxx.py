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





def auto_corr_test(p, d):

    #if p.do_testing=='False':
    #    return

    cor=cms.autocorr(d, auto_type='FFT')
    cb=pl.imshow(cor) #, norm=colors.LogNorm())
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

    print 'mpi.rank=', mpi.rank

    # ->> data importing <<- #
    root='../../workspace/'

    fn=root+'data/sims/dsU1.npy'
    d=np.load(fn)


    #->> 
    dout_name=root+'result/cr1d/dmap.dat'
    pk_name=root+'result/cr1d/bpk.dat'


    '''----------------------------------------------
                ->>    now we start    <<- 
       ----------------------------------------------'''
    # ->> lower data resolution <<- #
    #zoom_factor=0.5
    zoom_factor=1  #0.5
    n_pts=500
    ncol=10
    _dmap_=d[:n_pts,ncol]-np.mean(d[:n_pts,ncol])
    dmap=sp.ndimage.interpolation.zoom(_dmap_, zoom_factor)




    #->> data initialization <<- #
    qe_dict={'calculate_dcov':   False,
             'fname_dcov':       root+'result/r1d/dcov_r1d_fft_24.npz',
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

    bpk=np.concatenate((qe.klist[kdim/2+1:], qe.klist_low[kdim/2+1:], 
        qe.klist_up[kdim/2+1:], pk_fid), axis=0)

    #->> write files <<- #
    bpk.tofile(pk_name)
    dmap.tofile(dout_name)



    
    # ->> The End <<- #
    p.finalize()




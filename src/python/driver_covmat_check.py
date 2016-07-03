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
    #->> data initialization <<- #
    qe_dict={'calculate_dcov': False, 
             'fname_dcov':     root+'result/dcov.npz',
	    }

    # ->> lower data resolution <<- #
    zoom_factor=0.5
    _dmap_=d[:100,:100]-np.mean(d[:100,:100])
    dmap=sp.ndimage.interpolation.zoom(_dmap_, zoom_factor)

    raise Exception('update map resolution here, & redefine band pk etc.') 

    # ->> generate band power <<- #
    print 'map resolution', qe.dmap_res
    rsize=qe.dmap_res*np.array(dmap.shape)
    kdim=np.array(dmap.shape)  # ->> assuming full FFT instead of rfft <<- #
    klist=helper.klist_fft(rsize, kdim)

    print klist.shape, 'klist range:', min(klist[0]), max(klist[0]), min(klist[1]), max(klist[1])
    print 'Now, I got the full FFT klist; next lower the resolution, to get P(k).'
    print klist.shape

    
    raise Exception('update everything before initialize qe')


    ''' #->> Initialzing quadratic estimator class <<- #
        #->> parafname='same as parameter file' 
        #->> calculating dcov
    '''
    qe=qde.QuadestPara(paramfname=p.paramfname, section=p.qestmator_sec,
            prog_control=p, dmap=dmap, **qe_dict)

    print '\n->> QuadestPara parameters:\n', qe.paramdict
    print 'dcov shape:', qe.dcov.shape


    # ->> covariance matrix testing <<- #
    print '->> now test dcov and covariace matrix <<- #'
    print '->> Basically, assuming the correct power spectrum, I`d like to check whether dcov would produce the correct covariance matrix. <<- '


    # ->> measure pk and correlation function <<- #
    pk=cms.pk_fft_2d(dmap)
    cor=cms.autocorr(dmap, auto_type='FFT')

    print 'pk.shape, cor.shape:', pk.shape, cor.shape
    print 'pk min/max,  cor min/max: ', pk.min(), pk.max(), cor.min(), cor.max()

    _show_=False
    if _show_:
        nplt, ncol = 2, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                              gap_size=0.5,return_figure=True)

        cb1=ax[0].imshow(cor) 
        cb2=ax[1].imshow(pk, norm=colors.LogNorm()) 

        pl.show()





    #->> write files <<- #


    
    # ->> The End <<- #
    p.finalize()



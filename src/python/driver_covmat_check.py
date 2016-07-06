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

    # ->> measure pk and correlation function from FFT <<- #
    pk2d=cms.pk_fft_2d(dmap, qe.dmap_res)
    cor_fft=cms.autocorr(dmap, auto_type='FFT', zero_padding=True)

    print 'pk.shape, cor.shape:', pk2d.shape, cor_fft.shape
    print 'pk min/max,  cor_fft min/max: ', pk2d.min(), pk2d.max(), cor_fft.min(), cor_fft.max()


    # ->> covariance matrix testing <<- #
    print 'Assuming the correct power spectrum, to check if dcov would reproduce the correct covariance matrix. <<- '

    # ->> now check the covariance matrix <<- #
    corf=np.zeros((qe.m_dim[0], qe.m_dim[1]))
    cyth_cov.get_correlation(corf,  qe.dcov,  pk2d.flatten(), qe.npt)

    print 'cov comp:', corf.min(), corf.max(), cor_fft.min(), cor_fft.max()

    fname_cov_comp=root+'result/cov_comparison_50x50.npz'
    np.savez(fname_cov_comp, corf=corf, cor_fft=cor_fft)

    _show_=True
    if _show_:
        _t=np.arange(qe.m_dim[0])*qe.dmap_res[0]
        _f=np.arange(qe.m_dim[1])*qe.dmap_res[1]

	t, f = mar.meshgrid(_t, _f)

        nplt, ncol = 3, 3
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                              gap_size=0.5,return_figure=True)

        cb1=ax[0].pcolormesh(t, f, cor_fft, shading='gouraud') 
	cb1.set_edgecolor('face')
        pl.colorbar(cb1)


	cb2=ax[1].pcolormesh(t, f, corf, shading='gouraud')
	cb2.set_edgecolor('face')
        pl.colorbar(cb2)

        cb3=ax[2].imshow(pk2d, norm=colors.LogNorm()) 
	
	for i in range(2):
	    ax[i].set_xlim([0, _t[15]])
	    ax[i].set_ylim([0, _f[15]])

	    #ax[i].set_xlim([0, _t[-1]])
	    #ax[i].set_ylim([0, _f[-1]])

        pl.show()





    #->> write files <<- #


    
    # ->> The End <<- #
    p.finalize()




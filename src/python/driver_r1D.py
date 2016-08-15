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


    '''----------------------------------------------
                ->>    now we start    <<- 
       ----------------------------------------------'''
    # ->> lower data resolution <<- #
    zoom_factor=0.5
    _dmap_=d[:100,10]-np.mean(d[:100,10])
    dmap=sp.ndimage.interpolation.zoom(_dmap_, zoom_factor)

    if False:
        dk=np.fft.fft(dmap)
        pk=np.abs(np.fft.fftshift(dk))**2.

        pl.plot(pk)
	pl.show()


    #->> data initialization <<- #
    qe_dict={'calculate_dcov':   True,
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

    # ->> 
    qe.dcov.tofile(root+'result/r1d/dcov.dat')
    pk_fid.tofile(root+'result/r1d/plist.dat')
    dmap.tofile(root+'result/r1d/dmap.dat')


    # ->> check the covariance matrix <<- #
    print 'dcov', qe.dcov.shape
    pk_rec=np.zeros(pk_fid.shape)
    cyth_cov.get_correlation_r1d(pk_rec,  qe.dcov,  pk_fid)

    pk_rec.tofile(root+'result/r1d/pk_rec.dat')
    quit()

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
    fn_qi=root+'result/r1d/Qi_fft.npz'

    mpi.barrier()
    if mpi.rank0:
        np.savez(fn_qi, Qi=Qi, Fij=Fij)
     

    quit()



    # ->> initialize FFT band power <<- #
    #bp_dict={'dmap_shape':   qe.dmap.shape, } 
    #qe.band_power_init(**bp_dict)

    #print 'klist len:', len(qe.klist), 'dt:', qe.dt, 'npt:', qe.npt


    # ->> initialize dcov <<- #
    #fname_dcov_fft=root+'result/r1d/dcov_r1d_fft_49.npz'
    #qe.dcov_init(fname_dcov_fft)



    # ->> measure pk and correlation function from FFT <<- #
    pk=cms.pk_fft(dmap, qe.dmap_res)
    dk=np.fft.fft(dmap)
    pk=np.abs(np.fft.fftshift(dk))**2.


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


    _output_bp_dcov_=True
    if _output_bp_dcov_:
        fn_out_dcov=root+'result/dcov_out.dat'
        fn_out_plist=root+'result/plist.dat'

        (qe.dcov[:,:qe.m_dim[0],:qe.m_dim[1]]).tofile(fn_out_dcov) 
        (pk2d.flatten()).tofile(fn_out_plist) 


    print 'corf shape:', corf.shape, cor_fft.shape


    _show_=True
    if _show_:
        _t=np.arange(qe.m_dim[0])*qe.dmap_res[0]
        _f=np.arange(qe.m_dim[1])*qe.dmap_res[1]

	t, f = mar.meshgrid(_t, _f)

        nplt, ncol = 2, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                              gap_size=0.5,return_figure=True)

        cb1=ax[0].pcolormesh(t, f, cor_fft, shading='gouraud') 
	cb1.set_edgecolor('face')
        #pl.colorbar(cb1)


	cb2=ax[1].pcolormesh(t, f, corf, shading='gouraud')
	cb2.set_edgecolor('face')
        #pl.colorbar(cb2)

        #cb3=ax[2].imshow(pk2d, norm=colors.LogNorm()) 
        #pl.colorbar(cb3)
	
	for i in range(2):
	    ax[i].set_xlim([0, _t[15]])
	    ax[i].set_ylim([0, _f[15]])

	    #ax[i].set_xlim([0, _t[-1]])
	    #ax[i].set_ylim([0, _f[-1]])

        pl.show()





    #->> write files <<- #


    
    # ->> The End <<- #
    p.finalize()




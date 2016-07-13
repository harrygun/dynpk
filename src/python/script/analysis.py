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









param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'cosmology_parameter_fname': 'parameters/cosparameter.cfg',
    'cosmology_parameter_sec': 'Cosmology_Parameters',
    'qestmator_sec':       'Quadratic_Estimator',
    }

prog_control={
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

    fn_qi=root+'result/Qi_fft.npz'

    # ->> import dcov <<- #

    if False:

        fn_dcov=root+'result/dcov_out.dat'
        dcov=np.fromfile(fn_dcov).reshape(2500,50,50)

        print dcov.shape

        nplt, ncol = 4, 2
        fig,ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5.,\
                              gap_size=0.5,return_figure=True)

        ax[0].imshow(dcov[0])
        ax[1].imshow(dcov[1])
        ax[2].imshow(dcov[200])
        ax[3].imshow(dcov[-1])

	pl.show()


    if True:
        fn_dcov=root+'result/plist.dat'
        plist=np.fromfile(fn_dcov).reshape(50,50)

	print plist.shape

	cb=pl.imshow(plist, norm=colors.LogNorm())
	#pl.colorbar(cb)
	pl.show()


    # ->> The End <<- #
    p.finalize()
    mpi.finalize()




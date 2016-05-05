import os
import pylab as pl
import numpy as np

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import genscript.myarray as mar

import misc.helper as hp







def test(p):

    #size=(200,400)
    size=(400,900)

    #mask_type='default'  #'binary'
    mask_type='binary'

    m=hp.gen_mask(size, mask_type=mask_type)

    cb=pl.imshow(m)
    pl.colorbar(cb)
    pl.show()

    return





param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'a_init': 1e-2,
    'smooth_R': 0,
    'smooth_type': 'Gauss', 
    'smooth_R_list_type':  'linear'
    }

prog_control={
    'do_testing': False, 
    #-------------------------------#
    #-------------------------------#
    }





if __name__=='__main__':

    # ->> initialization <<- #
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(**init_dict)

    root='../../workspace/'


    # ->> data importing <<- #
    fn=root+'data/DynSpec/Rickett_53560dspec.npy'
    d=np.load(fn)
    print 'data shape:', d.shape


    #->> do some testing <<- #
    test(p)



    # ->> construction KL modes <<- #





    #->> write files <<- #





    
    # ->> The End <<- #
    p.finalize()




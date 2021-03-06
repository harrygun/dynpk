import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.fftpack as sfft
import scipy.linalg as slag
import json

from genscript.extendclass import *
import genscript.parameter as par
import genscript.myarray as mar
import genscript.mpiutil as mpi

import misc.helper as helper
import covmat as covm

import cyth.covm as cyth_cov
import cmeasure as cms






def one_dim_band_power_init(bp_init_type, fname=None, **pdict):

    if bp_init_type=='from_file':
        raise Exception('from_file band_power_init() NOT supported yet.')

    # ->>  <<- #
    if bp_init_type=='internal_log':
        raise Exception('NOT SUPPORTED NOW YET, NEED TO MODIFY BEFORE USE.')
        
        kt_min, kt_max, kt_num = pdict['kt_list_para']
        kf_min, kf_max, kf_num = pdict['kf_list_para']

        _kt_l=np.linspace(kt_min, kt_max, kt_num+1)
        _kf_l=np.linspace(kf_min, kf_max, kf_num+1)

        # ->> middle point <<- #
	kt_list=10.**((_kt_l[:-1]+_kt_l[1:])/2.)
	kf_list=10.**((_kf_l[:-1]+_kf_l[1:])/2.)

        klist=mar.meshgrid(kt_list, kf_list) 
        Dk_list=mar.meshgrid(10.**_kt_l[1:]-10.**_kt_l[:-1], \
                             10.**_kf_l[1:]-10.**_kf_l[:-1])

        klist_low=mar.meshgrid(10.**_kt_l[:-1], 10.**_kf_l[:-1])
        klist_up =mar.meshgrid(10.**_kt_l[1:],  10.**_kf_l[1:])

        sk=klist.shape
        #sdk=Dk_list.shape

	print 'Dk_list shape:', Dk_list.shape

        return klist.reshape(2,sk[1]*sk[2]), klist_low.reshape(2,sk[1]*sk[2]), \
               klist_up.reshape(2,sk[1]*sk[2]), kt_list, kf_list, \
               Dk_list.reshape(2,sk[1]*sk[2])


    if bp_init_type=='FFT':
        try:
            dmap_res=pdict['dmap_res']
	    dmap_shape=pdict['dmap_shape']
	except:
	    raise Exception('FFT band power initialization error')

        rbsize=dmap_res*np.array(dmap_shape)
        # ->> assuming full FFT instead of rfft <<- #
        kdim=np.array(dmap_shape)
        k_min=2.*np.pi/np.array(rbsize)

	print 'kdim', kdim

        klist=helper.klist_fft(rbsize, kdim)[0,int(kdim[0]/2)+1:]

        # ->> boundary <<- #
        klist_low=klist-k_min/2.
        klist_up =klist+k_min/2.

	#print 'klist_low', klist_low
	#print 'klist_up', klist_up

	#klist_low[0]=0.


        return klist, klist_low, klist_up







def quade_iter(dmap, dcov, covn_vec, pfid, klist, npt, npix, m_dim, nit=0, do_mpi=False):

    covf=np.zeros((npix, npix))

    Qi=np.zeros(npt)
    Fij=np.zeros((npt, npt))

    # ->> first run <<- #
    cyth_cov.quad_estimator_r1d_wrapper(dmap, covf, dcov, covn_vec, pfid, \
                                    Qi, Fij, npt, npix, m_dim, do_mpi=do_mpi)

    # ->> iteration <<- #
    for it in range(nit):
        cyth_cov.quad_estimator_r1d_wrapper(dmap, covf, dcov, covn_vec, Qi, \
                                        Qi, Fij, npt, npix, m_dim, do_mpi=do_mpi)

    else:
        raise Exception('ONLLY SUPPORT CYTHON VERSION.')

    return Qi, Fij





''' ------------------------------------------------------------ 
              ->>  Quadratic Estimator <<-
    ------------------------------------------------------------
'''
defaultQuadest1dParaValueDict={
    'map_dimension':   [100],
    'get_bandpower_list_type':    'from_file', 
    'bandpower_list_fname':       'x.dat', 
    'kt_list_para':               [-2, 2, 10],
    'map_resolution':             [1],
    'do_mpi':                     False,
    'calculate_dcov':             True,
    'fname_dcov':                 'y.dat',
    'map_zoom_factor':            0.5,
    'fiducial_pk_type':          'FFT_P(k)',
    }


defaultQuadest1dParaNameDict={
    'map_dimension':             'm_dim',
    'get_bandpower_list_type':   'get_bp_type', 
    'bandpower_list_fname':      'bp_list_fname',
    'kt_list_para':              'kt_list_para',
    'map_resolution':            'dmap_res',
    'do_mpi':                    'do_mpi',
    'calculate_dcov':            'calculate_dcov',
    'fname_dcov':                'fname_dcov',
    'map_zoom_factor':           'map_zoom_factor',
    'fiducial_pk_type':          'fid_pk_type',
    }


class Quadest1dPara(par.Parameters):

    def __init__(self, paramfname=None, section=None, def_var=True, 
                prog_control=None, dmap=None, skip_init=False, **pardict):

        super(Quadest1dPara,self).__init__(defaultQuadest1dParaValueDict, 
	                names_dict=defaultQuadest1dParaNameDict, paramfname=paramfname, 
		        section=section, def_var=def_var, **pardict)

	# ->> set the map data <<- #
	self.dmap=dmap
	#->> # of pixels <<- #
        self.npix=np.prod(dmap.shape)

        # ->> scale map resolution with zoom_factor <<- #
        _m_res_=[self.dmap_res[i]/self.map_zoom_factor \
                 for i in range(len(self.dmap_res)) ]

        if isinstance(dmap, np.ndarray):
	    # ->> update m_dim & dmap_res <<- #
	    update={'m_dim':      list(dmap.shape), 'dmap_res':  _m_res_ }
            self.update_params(**update)
            #print 'dmap is np.ndarray', self.m_dim, type(self.m_dim)

	self.dt=np.array(self.dmap_res)  # in unit of ...

        if skip_init==True:
	    return

	# ->> initialization band-power list <<- #
        self.band_power_init()

        # ->> initialization of quadratic estimator <<- #
        self.dcov_init(self.fname_dcov)

	#->> self.covn initialization <<- #
        #self.covn_vec_init()

	print 'initialization is done.'

        return


    def band_power_init(self, **paradict):
        # ->> # of band powers <<- #

	if self.get_bp_type=='FFT':
            #self.npt=np.prod(self.dmap.shape)

            _kdim=np.prod(self.dmap.shape)
            self.npt=int(_kdim/2)-1
	    paradict={'dmap_shape':   self.dmap.shape, }

        # ->> might need to modify this part later, DON'T think so though. <<-# 
        elif 'internal' in (self.get_bp_type.split('_')):  #=='internal_log':
	    self.npt=self.kt_list_para[-1]
	    raise Exception('NOT USED NOW.')

	else:
	    raise Exception()

        self.klist, self.klist_low, self.klist_up=one_dim_band_power_init(\
	                             self.get_bp_type, fname=self.bp_list_fname, \
			             **(myDict(self.paramdict)+myDict(paradict)) )

        return 


    def dcov_init(self, fn_dcov):

        if self.calculate_dcov==True:

            if self.do_mpi==True:
                _dcov_=covm.dcov_r1d(self.klist_low, self.klist_up, self.dt, \
                             self.npt, self.m_dim, speedup=True, do_mpi=self.do_mpi)
                self.dcov=mpi.gather_unify(_dcov_, root=0)
	    else:
                self.dcov=covm.dcov_r1d(self.klist_low, self.klist_up, self.dt, \
                             self.npt, self.m_dim, speedup=True, do_mpi=False)

            # ->> save data <<- #
            if mpi.rank0:
                np.savez(fn_dcov, dcov=self.dcov)
                mpi.barrier()

        else:
	    # ->> import from files <<- #
	    self.dcov=np.load(fn_dcov)['dcov']

        return


    def covn_vec_init(self, noise_level='noiseless'):
        self.covn_vec=covm.covn_vec(self.npix, noise_level=noise_level)
        return 



    def fid_pk_first_guess(self, guess_option=None):

        if guess_option==None:
	    guess_option=self.fid_pk_type

        if guess_option=='simplest':
            pk_fid=np.ones(self.npt)

        elif guess_option=='FFT_P(k)':

	    if self.get_bp_type!='FFT':  
	        raise Exception('Inconsistent fiducial pk setting.')


	    pk_fid=cms.pk_fft_1d(self.dmap, self.dmap_res).flatten() 

        return pk_fid




    def quadest_iteration(self, pk_fid, n_iteration):
        # ->> run iterative quadratic estimator <<- #
        return quade_iter(self.dmap, self.dcov, self.covn_vec, pk_fid, \
                          self.klist, self.npt, self.npix, self.m_dim, \
                          nit=n_iteration, do_mpi=self.do_mpi)

    def get_derived(self):
        pass




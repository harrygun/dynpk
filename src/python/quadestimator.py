import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.fftpack as sfft
import scipy.linalg as slag
import json

import genscript.parameter as par
import genscript.myarray as mar
import genscript.mpiutil as mpi

import misc.helper as helper
import covmat as covm

import cyth.covm as cyth_cov






def band_power_init(bp_init_type, fname=None, **pdict):

    if bp_init_type=='from_file':
        raise Exception('from_file band_power_init() NOT supported yet.')

    # ->>  <<- #
    if bp_init_type=='internal_log':
        
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
	except:
	    raise Exception('FFT band power initialization error')

        rsize=dmap_res*np.array(dmap.shape)

        # ->> assuming full FFT instead of rfft <<- #
        kdim=np.array(dmap.shape)

        k_list=helper.klist_fft(rsize, kdim)


	klist=mar.meshgrid( )

	Dk_list=None

        return  





"""->> NOT IN USE <<-"""
def quade_pk_single(dmap, covf, dcov, covn_vec, plist, klist, npt, npix, m_dim):
    ''' ->> construct fiducial estimator, given the fiducial pk <<- '''

    # ->> get full covariance matrix and inverse <<- #
    if not (covm.covfull(covf, dcov, covn_vec, plist, npt, npix, m_dim)):
        raise Exception('full covariance matrix error')
        
    icovf=slag.inv(covf)

    #->> quadratic estimator <<- #
    Qi=np.zeros(npt)

    for i in range(npt):
        Qi=np.einsum('ij,jk,ki', icovf, dcov[i], icovf)

    return Qi




def quade_iter(dmap, dcov, covn_vec, pfid, klist, npt, npix, m_dim, nit=0, do_cyth=True, do_mpi=False):

    covf=np.zeros((npix, npix))

    if do_cyth==True:

        # ->> first run <<- #
        if do_mpi==True:
            _Qi_=np.zeros(npt)
            cyth_cov.quad_estimator_wrapper(dmap, covf, dcov, covn_vec, pfid, \
                                            _Qi_, npt, npix, m_dim, do_mpi=True)
            Qi=mpi.gather_unify(_Qi_, root=0)


            # ->> iteration <<- #
            for it in range(nit):
                _Qi_=np.zeros(npt)
                cyth_cov.quad_estimator_wrapper(dmap, covf, dcov, covn_vec, Qi, \
                                                _Qi_, npt, npix, m_dim, do_mpi=True)
                Qi=mpi.gather_unify(_Qi_, root=0)
        else:

            Qi=np.zeros(npt)
            cyth_cov.quad_estimator_wrapper(dmap, covf, dcov, covn_vec, pfid, \
                                            Qi, npt, npix, m_dim)

            # ->> iteration <<- #
            for it in range(nit):
                cyth_cov.quad_estimator_wrapper(dmap, covf, dcov, covn_vec, Qi, \
                                                Qi, npt, npix, m_dim)

    else:
        # ->> first run <<- #
        Qi=quade_pk_single(dmap, covf, dcov, covn_vec, pfid, klist, npt, npix, m_dim)

        # ->> iteration <<- #
        for it in range(nit):
            Qi=quade_pk_single(dmap, covf, dcov, covn_vec, Qi, klist, npt, npix, m_dim)


    return Qi





''' ------------------------------------------------------------ 
              ->>  Quadratic Estimator <<-
    ------------------------------------------------------------
'''
defaultQuadestParaValueDict={
    'map_dimension':   [100, 100],
    'get_bandpower_list_type':    'from_file', 
    'bandpower_list_fname':       'x.dat', 
    'kt_list_para':               [-2, 2, 10],
    'kf_list_para':               [-2, 2, 10],
    'map_resolution':             [1, 1],
    'do_mpi':                     False,
    'calculate_dcov':             True,
    'fname_dcov':                 'y.dat',
    'map_zoom_factor':            0.5,
    }


defaultQuadestParaNameDict={
    'map_dimension':             'm_dim',
    'get_bandpower_list_type':   'get_bp_type', 
    'bandpower_list_fname':      'bp_list_fname',
    'kt_list_para':              'kt_list_para',
    'kf_list_para':              'kf_list_para',
    'map_resolution':            'dmap_res',
    'do_mpi':                    'do_mpi',
    'calculate_dcov':            'calculate_dcov',
    'fname_dcov':                'fname_dcov',
    'map_zoom_factor':           'map_zoom_factor',
    }


class QuadestPara(par.Parameters):

    def __init__(self, paramfname=None, section=None, def_var=True, 
                prog_control=None, dmap=None, skip_init=False, **pardict):

        super(QuadestPara,self).__init__(defaultQuadestParaValueDict, 
	                names_dict=defaultQuadestParaNameDict, paramfname=paramfname, 
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

	self.dt_df=np.array(self.dmap_res)  # in unit of ...

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


    def band_power_init(self):
        # ->> # of band powers <<- #
	self.npt=self.kt_list_para[-1]*self.kf_list_para[-1]

        self.klist, self.klist_low, self.klist_up, self.kt_list, \
            self.kf_list, self.Dk_list=band_power_init(self.get_bp_type, \
	    fname=self.bp_list_fname, **self.paramdict)

        return 

    def dcov_init(self, fn_dcov):

        if self.calculate_dcov==True:

            if self.do_mpi==True:
                _dcov_=covm.dcov(self.klist_low, self.klist_up, self.dt_df, \
                             self.npt, self.m_dim, speedup=True, do_mpi=self.do_mpi)
                self.dcov=mpi.gather_unify(_dcov_, root=0)
	    else:
                self.dcov=covm.dcov(self.klist_low, self.klist_up, self.dt_df, \
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



    def fid_pk_first_guess(self, guess_option='simplest'):

        if guess_option=='simplest':
            pk_fid=np.ones(self.npt)
        else:
            pk_fid=np.ones(self.npt)

        return pk_fid




    def quadest_iteration(self, pk_fid, n_iteration):
        # ->> run iterative quadratic estimator <<- #
        return quade_iter(self.dmap, self.dcov, self.covn_vec, pk_fid, \
                          self.klist, self.npt, self.npix, self.m_dim, \
                          nit=n_iteration, do_mpi=self.do_mpi)

    def get_derived(self):
        pass



''' ------------------------------------------------------------ 
              ->>                                 <<-
    ------------------------------------------------------------
'''



''' ->> other routines <<- '''
def dcov_to_cov(dcov, plist, klist):
    # ->> convert dcov and P(k) into cov <<- #

    return

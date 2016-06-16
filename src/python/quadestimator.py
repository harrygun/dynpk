import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.fftpack as sfft
import scipy.linalg as slag
import json

import genscript.parameter as par
import genscript.myarray as mar

import misc.helper as helper
import covmat as covm







def band_power_init(bp_init_type, fname=None, **pdict):

    if bp_init_type=='from_file':
        raise Exception('from_file band_power_init() NOT supported yet.')

    # ->>  <<- #
    if bp_init_type=='internal_log':
        
        kt_min, kt_max, kt_num = pdict['kt_list_para']
        kf_min, kf_max, kf_num = pdict['kf_list_para']

        kt_list=np.logspace(kt_min, kt_max, kt_num) 
        kf_list=np.logspace(kf_min, kf_max, kf_num) 

        klist=mar.meshgrid(kt_list, kf_list) 
	s=klist.shape

        return klist.reshape(2,s[1]*s[2])






def quade_pk_single(dmap, covf, dcov, covn_vec, plist, klist, npt, npix):
    ''' ->> construct fiducial estimator, given the fiducial pk <<- 
    '''
    
    qi=np.zeros(npt)

    # ->> get full covariance matrix and inverse <<- #
    if not (covm.covfull(covf, dcov, covn_vec, plist, klist, npt, npix)):
        raise Exception('full covariance matrix error')
        
    icovf=slag.inv(covf)

    #->> obtain the estimator <<- #
    for i in range(npt):
        qi=np.einsum('ij,jk,ki', (icovf, dcov[i], icovf) )

    return qi




def quade_iter(dmap, dcov, covn_vec, pfid, klist, npt, npix, nit=0):
    #raise Exception()

    covf=np.zeros((npix, npix))

    # ->> first run <<- #
    qi=quade_pk_single(dmap, covf, dcov, covn_vec, pfid, klist, npt, npix)

    # ->> iteration <<- #
    for it in range(nit):
        qi=quade_pk_single(dmap, covf, dcov, covn_vec, qi, klist, npt, npix)


    return qi





''' ------------------------------------------------------------ 
              ->>  Quadratic Estimator <<-
    ------------------------------------------------------------
'''
defaultQuadestParaValueDict={
    'map_dimension':   [100, 100],
    'get_bandpower_list_type':   'from_file', 
    'bandpower_list_fname':      'x.dat', 
    'kt_list_para':               [-2, 2, 10],
    'kf_list_para':               [-2, 2, 10],
    'map_resolution':             [1, 1],
    }


defaultQuadestParaNameDict={
    'map_dimension':             'm_dim',
    'get_bandpower_list_type':   'get_bp_type', 
    'bandpower_list_fname':      'bp_list_fname',
    'kt_list_para':              'kt_list_para',
    'kf_list_para':              'kf_list_para',
    'map_resolution':            'dmap_res',
    }


class QuadestPara(par.Parameters):

    def __init__(self, paramfname=None, section=None, def_var=True, 
                 prog_control=None, dmap=None, **pardict):

        super(QuadestPara,self).__init__(defaultQuadestParaValueDict, 
	                names_dict=defaultQuadestParaNameDict, paramfname=paramfname, 
		        section=section, def_var=def_var, **pardict)

	# ->> set the map data <<- #
	self.dmap=dmap
	#->> # of pixels <<- #
        self.npix=np.prod(dmap.shape)

        if isinstance(dmap, np.ndarray):

	    # ->> update m_dim <<- #
	    update={'m_dim': list(dmap.shape)}
            self.update_params(**update)

            #print 'dmap is np.ndarray', self.m_dim, type(self.m_dim)

	# ->> initialization band-power list <<- #
        self.band_power_init()

        # ->> initialization of quadratic estimator <<- #
        self.dcov_init()

	#->> self.covn initialization <<- #
	self.covn_vec_init()

        return


    def band_power_init(self):
        # ->> # of band powers <<- #
	self.npt=self.kt_list_para[-1]*self.kf_list_para[-1]

        self.klist=band_power_init(self.get_bp_type, fname=self.bp_list_fname, \
	                           **self.paramdict)
        return 

    def dcov_init(self):
        self.dcov=covm.dcov(self.klist, self.plist, self.npt, self.npix)
        return

    def covn_vec_init(self):
        return covm.covn_vec()


    def quadest_iteration(self, pk_fid, n_iteration):
        # ->> run iterative quadratic estimator <<- #
        return quade_iter(self.dmap, self.dcov, self.covn_vec, \
               pk_fid, self.klist, self.npt, self.npix, nit=n_iteration)

    def get_derived(self):
        pass



''' ------------------------------------------------------------ 
              ->>                                 <<-
    ------------------------------------------------------------
'''






''' ->> other routines <<- '''
def autocorr(d, auto_type='FFT'):
    ''' ->> conduct the auto-correlation of the map <<- ''' 

    if auto_type=='scipy':
        cor=sig.correlate2d(d, d) 

    if auto_type=='FFT':
       # ->> zero padding <<- #
       dx, dy=d.shape
       dd=np.zeros((dx*2, dy*2))
       dd[:dx,:dy]=d


       dk=np.fft.rfft2(dd)
       cor=np.fft.irfft2(dk*np.conjugate(dk))
       cor=np.fft.fftshift(cor)


    return cor





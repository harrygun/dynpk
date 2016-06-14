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

        plist=mar.meshgrid(kt_list, kf_list) 

        return plist





def quade_init(qe):
    # ->> initialization of quadratic estimator <<- #
    # ->> 

    return



def quade_pk(epar, d, fcm):
    ''' ->> construct fiducial estimator, given the fiducial pk <<- 
            epar:  parameters for quadratic estimator 
            fcm: fidicual covariance matrix  
    '''

    # ->> band power initialization <<- #
    #plist=


    # ->> obtain the derivative of covariance matrix <<- #
    dcov=covm.dcov(plist, pdiff)
    
    # ->> obtain the full covariance matrix <<- #
    fcov=covm.covfull(epar, d)
    fcov_inv=slag.inv(fcov)


    #->> obtain the estimator <<- #
    for i in range(npt):
        qi=np.einsum('ij,jk,ki', (fcov_inv, dcov[i], fcov_inv) )


    return








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
    }


defaultQuadestParaNameDict={
    'map_dimension':             'm_dim',
    'get_bandpower_list_type':   'get_bp_type', 
    'bandpower_list_fname':      'bp_list_fname',
    'kt_list_para':               'kt_list_para',
    'kf_list_para':               'kf_list_para',
    }

class QuadestPara(par.Parameters):

    def __init__(self, paramfname=None, section=None, def_var=True, 
                 prog_control=None, dmap=None, **pardict):

        super(QuadestPara,self).__init__(defaultQuadestParaValueDict, 
	                names_dict=defaultQuadestParaNameDict, paramfname=paramfname, 
		        section=section, def_var=def_var, **pardict)

	# ->> set the map data <<- #
	self.dmap=dmap

        if isinstance(dmap, np.ndarray):
	    self.m_dim=dmap.shape

	# ->> initialization <<- #
        self.band_power_init()

        return

    def get_derived(self):
        pass


    def band_power_init(self):
        # ->> 
        band_power_init(self.get_bp_type, fname=self.bp_list_fname, **self.paramdict)
        return

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





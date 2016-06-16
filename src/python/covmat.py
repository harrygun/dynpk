import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.special as spec

import genscript.parameter as par
import misc.helper as helper



def SinIntegral(x):
    return spec.sici(x)[0]

def Cos(x):
    return np.cos(x)


def dcov_1D(kti, Dkti, dt, dtab):
    ''' ->> kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    '''

    dcov1d=((-4*Cos(dtab*(-Dkti/2. + kti)))/(Dkti - 2*kti) + \
          (4*Cos(dt*(-Dkti/2. + kti))*Cos(dtab*(-Dkti/2. + kti)))/(Dkti-2*kti) - \
          (4*Cos(dtab*(Dkti/2. + kti)))/(Dkti + 2*kti) + \
          (4*Cos(dt*(Dkti/2. + kti))*Cos(dtab*(Dkti/2. + kti)))/(Dkti + 2*kti) + \
          (-dt + dtab)*SinIntegral((dt - dtab)*(-Dkti/2. + kti)) + \
          2*dtab*SinIntegral(dtab*(-Dkti/2. + kti)) - \
          dt*SinIntegral((dt + dtab)*(-Dkti/2. + kti)) - \
          dtab*SinIntegral((dt + dtab)*(-Dkti/2. + kti)) + \
          (dt - dtab)*SinIntegral((dt - dtab)*(Dkti/2. + kti)) - \
          2*dtab*SinIntegral(dtab*(Dkti/2. + kti)) + \
          dt*SinIntegral((dt + dtab)*(Dkti/2. + kti)) + \
          dtab*SinIntegral((dt + dtab)*(Dkti/2. + kti)))/np.pi


    return dcov1d





def dcov(klist, plist, dt, npt, npix, speedup=False):
    ''' ->> get the derivative of covariance matrix <<- 
    '''
    dcov=np.zeros(npt, npix, npix)

    if speedup==True:
        # ->> C loop <<- #
	return 

    else:
        # ->> python loop <<- #
        for i in range(npt):
            kti, kfi  =
            Dkti, Dkfi=

            for a in range(npix):
	        for b in range(npix):
                    tab=

                    dcov[i,a,b]=dcov_1D(kti, Dkti, dt, dtab)*\
                                dcov_1D(kfi, Dkfi, df, dfab)


    return dcov




def covn_vec():

    return


def covfull(covf, dcov, covn_vec, plist, klist, npt, npix):
    # ->>  from dcov and covn_vec, get the full covariance matrix <<- #

    for i in range(npt):
        covf=dcov[i]*plist[i]

    for a in range(npix):
        covf[a,a]+=covn_vec[a]

    return True




def fisher(icov, dcov, npt, npix):
    # ->> calculate the fisher matrix <<- #

    F=np.zeros((npt, npt))


    return



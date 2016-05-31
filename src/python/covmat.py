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
    ''' ->> kt_b:  boundary in k_t direction <<- 
            kf_b:  boundary in k_f direction <<- 
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


    return





def dcov(plist, pdiff, ):
    ''' ->> get the derivative of covariance matrix <<- '''

    npt=len(plist)
    dcov=np.zeros(npt, npix, npix)

    for i in range(npt):
        # ->> 
	kti, kfi  =
        Dkti, Dkfi=

        for a in range(npix):
	    for b in range(npix):
                dcov[i,a,b]=dcov_1D(kti, Dkti, dt, dtab)*
                            dcov_1D(kfi, Dkfi, df, dfab)

    return dcov




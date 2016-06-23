import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.special as spec

import genscript.parameter as par
import misc.helper as helper



def SinIntegral(x):
    return spec.sici(x)[0]

def CosIntegral(x):
    return spec.sici(x)[1]

def Cos(x):
    return np.cos(x)

def Sin(x):
    return np.sin(x)

def dcov_1D_real(kti, Dkti, dt, dtab):
    ''' ->> kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    '''

    '''
    dc1d_real=((-4*Cos(dtab*(-Dkti/2. + kti)))/(Dkti - 2*kti) + \
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
    '''

    dc1d_real=((4*(-1 + Cos(dt*(-Dkti/2. + kti)))*Cos(dtab*(-Dkti/2. + kti)))/ \
                (Dkti - 2*kti) + (4*(-1 + Cos(dt*(Dkti/2. + kti)))*  \
                      Cos(dtab*(Dkti/2. + kti)))/(Dkti + 2*kti) + \
                   (-dt + dtab)*SinIntegral((dt - dtab)*(-Dkti/2. + kti)) + \
                   2*dtab*SinIntegral(dtab*(-Dkti/2. + kti)) -  \
                   dt*SinIntegral((dt + dtab)*(-Dkti/2. + kti)) -  \
                   dtab*SinIntegral((dt + dtab)*(-Dkti/2. + kti)) +  \
                   dt*SinIntegral((dt - dtab)*(Dkti/2. + kti)) -   \
                   dtab*SinIntegral((dt - dtab)*(Dkti/2. + kti)) -   \
                   2*dtab*SinIntegral(dtab*(Dkti/2. + kti)) +   \
                   (dt + dtab)*SinIntegral((dt + dtab)*(Dkti/2. + kti)))/np.pi

    return dc1d_real



def dcov_1d_imag(kti, Dkti, dt, dtab):
    ''' ->> kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    '''
        dc1d_imag= ((dt - dtab)*(Dkti - 2*kti)*(Dkti + 2*kti)*\
             CosIntegral((dt - dtab)*(-Dkti/2. + kti)) +  \
            2*dtab*(Dkti - 2*kti)*(Dkti + 2*kti)*CosIntegral(dtab*(-Dkti/2. + kti)) - \
            dt*(Dkti - 2*kti)*(Dkti + 2*kti)*CosIntegral((dt + dtab)*(-Dkti/2. + kti)) - \
            dtab*(Dkti - 2*kti)*(Dkti + 2*kti)* \
             CosIntegral((dt + dtab)*(-Dkti/2. + kti)) - \ 
            (dt - dtab)*(Dkti - 2*kti)*(Dkti + 2*kti)* \
             CosIntegral((dt - dtab)*(Dkti/2. + kti)) -  \
            2*dtab*(Dkti - 2*kti)*(Dkti + 2*kti)*CosIntegral(dtab*(Dkti/2. + kti)) + \
            dt*(Dkti - 2*kti)*(Dkti + 2*kti)*CosIntegral((dt + dtab)*(Dkti/2. + kti)) + \
            dtab*(Dkti - 2*kti)*(Dkti + 2*kti)* \
             CosIntegral((dt + dtab)*(Dkti/2. + kti)) + \
            4*(Dkti + 2*kti)*Sin(dtab*(-Dkti/2. + kti)) - \ 
            4*(Dkti + 2*kti)*Cos(dt*(-Dkti/2. + kti))*Sin(dtab*(-Dkti/2. + kti)) + \
            4*(Dkti - 2*kti)*Sin(dtab*(Dkti/2. + kti)) -  \
            4*(Dkti - 2*kti)*Cos(dt*(Dkti/2. + kti))*Sin(dtab*(Dkti/2. + kti)))/ \
            (4.*(-(Dkti**2*Pi)/4. + kti**2*Pi))

    return dc1d_imag



def dcov(klist, Dk_list dt_df, npt, m_dim, speedup=False):
    ''' ->> get the derivative of covariance matrix <<- 
    '''
    raise Exception('restart from here. ')

    dcov=np.zeros(npt, 2*m_dim[0], 2*m_dim[1])
    dt, tf = dt_df

    if speedup==True:
        # ->> C loop <<- #
	return 

    else:
        # ->> python loop <<- #
        for i in range(npt):
            kti,  kfi  = klist[:,i]
            Dkti, Dkfi = Dk_list[:,i]

            for a in range(2*m_dim[0]):
	        for b in range(2*m_dim[1]):
                    dtab = a*dt
                    dfab = b*df
                    dcov[i,a,b]=2*(dcov_1D_real(kti, Dkti, dt, dtab)*  \
                                   dcov_1D_real(kfi, Dkfi, df, dfab) - \ 
                                   dcov_1D_imag(kti, Dkti, dt, dtab)*  \
                                   dcov_1D_imag(kfi, Dkfi, df, dfab)   )
    return dcov




def covn_vec(npix, noise_method='by_hand'):

    if noise_method=='by_hand':
        covn=np.ones(npix)

    return covn


def covfull(covf, dcov, covn_vec, plist, klist, npt, npix):
    # ->>  from dcov and covn_vec, get the full covariance matrix <<- #


    for i in range(npt):
        _covf=np.zeros((npix, npix))
        covf=dcov[i]*plist[i]

    for a in range(npix):
        covf[a,a]+=covn_vec[a]

    return True




def fisher(icov, dcov, npt, npix):
    # ->> calculate the fisher matrix <<- #

    F=np.zeros((npt, npt))


    return



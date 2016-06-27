import numpy as np
import pylab as pl
import scipy.signal as sig
import scipy.special as spec

import genscript.parameter as par
import misc.helper as helper

import cyth.covm as cyth_cov



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



def dcov_1D_imag(kti, Dkti, dt, dtab):
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
            (4.*(-(Dkti**2*np.pi)/4. + kti**2*np.pi))

    return dc1d_imag



def _dcov_NOT_use_(klist, Dk_list, dt_df, npt, m_dim, speedup=True, do_mpi=False):
    ''' ->> get the derivative of covariance matrix <<- 
    '''
    dcov=np.zeros((npt, 2*m_dim[0], 2*m_dim[1]))
    #print npt, m_dim, Dk_list.shape, klist.shape, dt_df.shape

    if speedup==True:
        # ->> C loop <<- #
        cyth_cov.get_dcov_klim(dcov, klist, Dk_list, dt_df, npt, m_dim, do_mpi=do_mpi)
	return dcov

    else:
        dt, df = dt_df

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



def dcov(klist_low, klist_up, dt_df, npt, m_dim, speedup=True, do_mpi=False):
    ''' ->> get the derivative of covariance matrix <<- 
    '''
    dcov=np.zeros((npt, 2*m_dim[0], 2*m_dim[1]))

    if speedup==True:
        # ->> cython loop <<- #
        cyth_cov.get_dcov_klim(dcov, klist_low, klist_up, dt_df, npt, \
	                       m_dim, do_mpi=do_mpi)
    else:
        raise Exception

    return dcov




def covn_vec(npix, noise_method='by_hand', noise_level='noiseless'):

    if noise_method=='by_hand':
        if noise_level=='noiseless':
            covn=np.zeros(npix)
	else:
            covn=np.ones(npix)
    else:
        raise Exception

    return covn


def covfull(covf, dcov, covn_vec, plist, npt, npix, m_dim, speedup=True):
    # ->>  from dcov and covn_vec, get the full covariance matrix <<- #

    if speedup==True:
        # ->>
        #print 'covfull ---  dcov shape:', dcov.shape
        cyth_cov.convert_cov_full(covf, dcov, covn_vec, plist, npt, npix, m_dim)

    else:
        for i in range(npt):
            covf=dcov[i]*plist[i]

        for a in range(npix):
            covf[a,a]+=covn_vec[a]

    return True




def fisher(icov, dcov, npt, npix):
    # ->> calculate the fisher matrix <<- #

    F=np.zeros((npt, npt))


    return



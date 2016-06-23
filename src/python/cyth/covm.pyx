import numpy as np
import pylab as pl
import scipy.sparse as spm
import scipy.special as spec

import genscript.mpiutil as mpi
import genscript.matrix_banded as mband
import genscript.myarray as mar

cimport numpy as cnp



cdef double SinIntegral(double x):
    return spec.sici(x)[0]

cdef double CosIntegral(double x):
    return spec.sici(x)[1]

cdef double Cos(double x):
    return np.cos(x)

cdef double Sin(double x):
    return np.sin(x)




cdef double dcov_1D_real(double kti, double Dkti, double dt, double dtab):
    ''' ->> kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    '''
    cdef double dc1d_real

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



cdef dcov_1D_imag(double kti, double Dkti, double dt, double dtab):
    ''' ->> kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    '''
    cdef double dc1d_imag

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






cpdef get_dcov(dcov, klist, Dk_list, dt_df, npt, m_dim):

    ''' ->> get the derivative of covariance matrix <<- 
    '''
    #dcov=np.zeros((npt, 2*m_dim[0], 2*m_dim[1]))

    cdef:
        int i, a, b
        double dt, df
    dt, df = dt_df

    print 'get_dcov:', npt, m_dim, Dk_list.shape, klist.shape, dt_df.shape


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





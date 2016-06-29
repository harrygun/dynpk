import numpy as np
import pylab as pl
import scipy.sparse as spm
import scipy.special as spec
import scipy.linalg as slag

import genscript.mpiutil as mpi
import genscript.matrix_banded as mband
import genscript.myarray as mar

cimport numpy as cnp
from libc.stdlib cimport malloc, free




cdef double _epsilon_ = 1e-10

cdef int mytrue=1
cdef int myfalse=0


cdef double SinIntegral(double x):
    cdef double si
    si=spec.sici(x)[0]

    if (np.isnan(si))|(np.isinf(si)):
        print 'SineIntegral nan/inf:', x
    return si

cdef double CosIntegral(double x):
    cdef double ci
    ci=spec.sici(x)[1]

    if (np.isnan(ci))|(np.isinf(ci)):
        print 'CosIntegral nan/inf:', x
    return ci


cdef double Cos(double x):
    return np.cos(x)


cdef double Sin(double x):
    return np.sin(x)


cdef void mat_multiplication(cnp.ndarray[cnp.double_t, ndim=2] m1, \
                             cnp.ndarray[cnp.double_t, ndim=2] m2, \
                             cnp.ndarray[cnp.double_t, ndim=2] m):
    cdef int i, j, k;

    return




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

    if np.isnan(dc1d_real):
        print 'real:', kti, Dkti, dt, dtab

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


    if np.isnan(dc1d_imag):
        print 'imag:', kti, Dkti, dt, dtab

    return dc1d_imag



cpdef get_dcov(dcov, klist, Dk_list, dt_df, npt, m_dim, do_mpi=False):

    ''' ->> get the derivative of covariance matrix <<- 
    '''
    #dcov=np.zeros((npt, 2*m_dim[0], 2*m_dim[1]))

    cdef:
        int i, a, b
        double dt, df

    dt, df = dt_df

    print 'get_dcov:', npt, m_dim, Dk_list.shape, klist.shape, dt_df.shape

    if do_mpi==True:
        prange=mpi.mpirange(npt)
    else:
        prange=range(npt)

    # ->> python loop <<- #
    for i in prange:
        kti,  kfi  = klist[:,i]
        Dkti, Dkfi = Dk_list[:,i]

        print i, kti, kfi, Dkti, Dkfi, dt, df, dtab, dfab
    
        for a in range(2*m_dim[0]):
            for b in range(2*m_dim[1]):
                dtab = a*dt
                dfab = b*df
                dcov[i,a,b]=2*(dcov_1D_real(kti, Dkti, dt, dtab)*  \
                               dcov_1D_real(kfi, Dkfi, df, dfab) - \
                               dcov_1D_imag(kti, Dkti, dt, dtab)*  \
                               dcov_1D_imag(kfi, Dkfi, df, dfab)   )
    return dcov




'''------------------------------------------------------------------------------'''
'''------------------------------------------------------------------------------'''
cdef double dcov1d_klim_real(double ktia, double ktib, double dt, double dtab):
    ''' ->> kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    '''
    cdef double dc1d_real

    dc1d_real= (-2*ktib*(-1 + Cos(dt*ktia))*Cos(dtab*ktia) + \
             2*ktia*(-1 + Cos(dt*ktib))*Cos(dtab*ktib) + \
             (-dt + dtab)*ktia*ktib*SinIntegral((dt - dtab)*ktia) + \
             2*dtab*ktia*ktib*SinIntegral(dtab*ktia) - \
             ktia*ktib*((dt + dtab)*SinIntegral((dt + dtab)*ktia) + \
             (-dt + dtab)*SinIntegral((dt - dtab)*ktib) + 2*dtab*SinIntegral(dtab*ktib)\
             ) + (dt + dtab)*ktia*ktib*SinIntegral((dt + dtab)*ktib))/(ktia*ktib*np.pi)
              
    if np.isnan(dc1d_real):
        print 'real:', ktia, ktib, dt, dtab

    return dc1d_real


cdef double dcov1d_klim_imag(double ktia, double ktib, double dt, double dtab):
    ''' ->> kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    '''
    cdef double dc1d_imag

    if (np.fabs(dtab-dt)<_epsilon_):
        # ->> if tab->dt
        dc1d_imag= (-2*dt*ktia*ktib*CosIntegral(dt*ktia) \
                  + 2*dt*ktia*ktib*CosIntegral(2*dt*ktia) + \
                  2*dt*ktia*ktib*CosIntegral(dt*ktib) - \
                  2*dt*ktia*ktib*CosIntegral(2*dt*ktib) + \
                  2*ktib*Sin(dt*ktia) - ktib*Sin(2*dt*ktia) - 2*ktia*Sin(dt*ktib) + \
                  ktia*Sin(2*dt*ktib))/(ktia*ktib*np.pi)

    elif (np.fabs(dtab+dt)<_epsilon_):
        #->> dtab->-dt
        dc1d_imag=(2*dt*ktia*ktib*CosIntegral(-(dt*ktia)) - \
                  2*dt*ktia*ktib*CosIntegral(2*dt*ktia) - \
                  2*dt*ktia*ktib*CosIntegral(-(dt*ktib)) + \
                  2*dt*ktia*ktib*CosIntegral(2*dt*ktib) - 2*ktib*Sin(dt*ktia) + \
                  ktib*Sin(2*dt*ktia) + 2*ktia*Sin(dt*ktib) - \
                  ktia*Sin(2*dt*ktib))/(ktia*ktib*np.pi)

    elif (np.fabs(dtab)<_epsilon_):
        # ->> dtab->0
        dc1d_imag=0.

    else:
        dc1d_imag=((-dt + dtab)*ktia*ktib*CosIntegral((dt - dtab)*ktia) - \
                  2*dtab*ktia*ktib*CosIntegral(dtab*ktia) + \
                  dt*ktia*ktib*CosIntegral((dt + dtab)*ktia) + \
                  dtab*ktia*ktib*CosIntegral((dt + dtab)*ktia) + \
                  (dt - dtab)*ktia*ktib*CosIntegral((dt - dtab)*ktib) + \
                  2*dtab*ktia*ktib*CosIntegral(dtab*ktib) - \
                  dt*ktia*ktib*CosIntegral((dt + dtab)*ktib) - \
                  dtab*ktia*ktib*CosIntegral((dt + dtab)*ktib) + 2*ktib*Sin(dtab*ktia)-\
                  2*ktib*Cos(dt*ktia)*Sin(dtab*ktia) - 2*ktia*Sin(dtab*ktib) + \
                  2*ktia*Cos(dt*ktib)*Sin(dtab*ktib))/(ktia*ktib*np.pi)



    if np.isnan(dc1d_imag):
        print 'imag:', ktia, ktib, dt, dtab

    return dc1d_imag




cpdef get_dcov_klim(dcov, klist_low, klist_up, dt_df, npt, m_dim, do_mpi=False):

    ''' ->> get the derivative of covariance matrix <<- '''

    cdef:
        int i, a, b
        double dt, df

    dt, df = dt_df

    print 'get_dcov:', npt, m_dim, klist_up.shape, klist_low.shape, dt_df.shape

    if do_mpi==True:
        prange=mpi.mpirange(npt)
    else:
        prange=range(npt)


    # ->> python loop <<- #
    for i in prange:
        #kti,  kfi  = klist[:,i]
        ktia, kfia = klist_low[:,i]
        ktib, kfib = klist_up[:,i]

        print i, ktia, ktib, kfia, kfib, dt, df
    
        for a in range(-m_dim[0], m_dim[0]):
            for b in range(-m_dim[1], m_dim[1]):
                dtab = a*dt
                dfab = b*df

                dcov[i,a,b]=2*(dcov1d_klim_real(ktia, ktib, dt, dtab)*  \
                               dcov1d_klim_real(kfia, kfib, df, dfab) - \
                               dcov1d_klim_imag(ktia, ktib, dt, dtab)*  \
                               dcov1d_klim_imag(kfia, kfib, df, dfab)   )
			    
    return dcov




'''------------------------------------------------------------------------------'''
'''------------------------------------------------------------------------------'''
cdef void mpixel_idx(int a, int mdim_t, int mdim_f, int *idx):

    idx[0]=<int>(a/<double>mdim_t) 
    idx[1]=a-mdim_t*idx[0]

    return  

cdef int mpixel_idx_inv(int idx_t, int idx_f, int mdim_t, int mdim_f):
    return idx_t*mdim_t+idx_f


cdef void full_cov_recovery(cnp.ndarray[cnp.double_t, ndim=2] covf, 
                       cnp.ndarray[cnp.double_t, ndim=3] dcov, 
                       cnp.ndarray[cnp.double_t, ndim=1] covn_vec, 
                       cnp.ndarray[cnp.double_t, ndim=1] plist, 
                       int npt, int npix, int mdim_t, int mdim_f):

    cdef: 
        int i, a, b, *idx_a, *idx_b

    idx_a=<int *>malloc(2*sizeof(int))
    idx_b=<int *>malloc(2*sizeof(int))

    # ->>  signal covariance matrix <<- #
    for a in range(npix):
        mpixel_idx(a, mdim_t, mdim_f, idx_a)

        for b in range(npix):
            mpixel_idx(b, mdim_t, mdim_f, idx_b)

            for i in range(npt):
                covf[a,b]+=dcov[i,idx_a[0]-idx_b[0], idx_a[1]-idx_b[1]]*plist[i]

        covf[a,a]+=covn_vec[a]


    free(idx_a)
    free(idx_b)

    return


cpdef convert_cov_full(covf, dcov, covn_vec, plist, npt, npix, m_dim):

    #print '--->>', covf.shape, dcov.shape, covn_vec.shape, plist.shape

    full_cov_recovery(covf, dcov, covn_vec, plist, <int>npt, <int> npix, \
                      <int> m_dim[0], <int> m_dim[1])

    return




''' ->> Quadratic Estimator <<- #'''
cdef void dcov_c_multiple(cnp.ndarray[cnp.double_t, ndim=2] icovf, 
                          cnp.ndarray[cnp.double_t, ndim=3] dcov,
                          cnp.ndarray[cnp.double_t, ndim=2] out, int i, int npix):

    # ->> give the multiplication dcov^i_ac C^-1_cb <<- #

    cdef: 
        int a, b, c
        int *idx_a

    for a in range(npix):
        for b in range(npix):

            for c in range(npix):
                out[a,b]+=  
        


    return


cdef void quad_estimator(cnp.ndarray[cnp.double_t, ndim=2] dmap, 
                        cnp.ndarray[cnp.double_t, ndim=2] covf, 
                        cnp.ndarray[cnp.double_t, ndim=3] dcov, 
                        cnp.ndarray[cnp.double_t, ndim=1] covn_vec, 
                        cnp.ndarray[cnp.double_t, ndim=1] plist, 
                        cnp.ndarray[cnp.double_t, ndim=1] Qi, 
                        int npt, int npix, int mdim_t, int mdim_f, int do_mpi):

    cdef: 
        int i, a, b, c, d, e, f, g
        int *idx_a, *idx_b, *idx_c, *idx_d


    full_cov_recovery(covf, dcov, covn_vec, plist, npt, npix, mdim_t, mdim_f)
    icovf=slag.inv(covf)
    mpi.barrier()

    print '->> preparation of icov & cov is done, now calculate Qe.', mpi.rank

    idx_a=<int *>malloc(2*sizeof(int))
    idx_b=<int *>malloc(2*sizeof(int))
    idx_c=<int *>malloc(2*sizeof(int))
    idx_d=<int *>malloc(2*sizeof(int))


    if do_mpi==True:
        prange=mpi.mpirange(npt)
    else:
        prange=range(npt)


    # ->> now start <<- #
    for i in prange:
        print 'quadratic estimator: rank-', mpi.rank, '  prange:', i

        for a in range(npix):
            mpixel_idx(a, mdim_t, mdim_f, idx_a)

            for b in range(npix):
                mpixel_idx(b, mdim_t, mdim_f, idx_b)

                for c in range(npix):
                    mpixel_idx(c, mdim_t, mdim_f, idx_c)

                    for d in range(npix):
                        mpixel_idx(d, mdim_t, mdim_f, idx_d)

                        Qi[i]+=dmap[idx_a[0],idx_a[1]]*icovf[a,b]*dcov[i,\
                               idx_b[0]-idx_c[0],idx_b[1]-idx_c[1]]*icovf[c,d]\
                               *dmap[idx_d[0],idx_d[1]]

    free(idx_a)
    free(idx_b)
    free(idx_c)
    free(idx_d)



    return




cpdef quad_estimator_wrapper(dmap, covf, dcov, covn_vec, plist, Qi, npt, npix, m_dim, do_mpi=False):

    cdef int dompi

    if do_mpi==True:
        dompi=mytrue
    else:
        dompi=myfalse

    quad_estimator(dmap, covf, dcov, covn_vec, plist, Qi, <int> npt, \
                   <int> npix, <int> m_dim[0], <int> m_dim[1], dompi)

    return 



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


cdef int matidx2(int i, int j, int ndim):
    return i*ndim+j

cdef int matidx3(int i, int j, int k, int ni, int nj, int nk):
    return i*nj*nk+j*nk+k


cdef void mat_multiplication(cnp.ndarray[cnp.double_t, ndim=2] m1, \
                             cnp.ndarray[cnp.double_t, ndim=2] m2, \
                             double *m, int ndim):
    cdef int i, j, k;

    for i in range(ndim):
        for j in range(ndim):

            m[matidx2(i,j,ndim)]=0.

            for k in range(ndim):
                m[matidx2(i,j,ndim)]+=m1[i,k]*m2[k,j] 
    return



cdef void mat_mult_vec(cnp.ndarray[cnp.double_t, ndim=2] mat, \
                       cnp.ndarray[cnp.double_t, ndim=1] vec, \
                       double *out, int ndim):
    # ->> calcuate out[i]=sum_j (m[i,j]v[j])

    cdef int i, j;

    for i in range(ndim):
        out[i]=0.
        for j in range(ndim):
            out[i]+=mat[i,j]*vec[j] 

    return


cdef void mpixel_idx(int a, int mdim_t, int mdim_f, int *idx):

    idx[0]=<int>(a/<double>mdim_t) 
    idx[1]=a-mdim_t*idx[0]

    return  

cdef int mpixel_idx_inv(int idx_t, int idx_f, int mdim_t, int mdim_f):
    return idx_t*mdim_t+idx_f





'''---------------------------------------------------------------------- '''

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

    return dc1d_real/2.



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

    return dc1d_imag/2.



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

    return dc1d_real/2.


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

    return dc1d_imag/2.




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






cpdef get_dcov_klim_1D(dcov, klist_low, klist_up, dt, npt, m_dim, do_mpi=False):

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

                dcov[i,a,b]=2*dcov1d_klim_real(ktia, ktib, dt, dtab)
			    
    return dcov






'''------------------------------------------------------------------------------'''
'''------------------------------------------------------------------------------'''

cdef void full_cov_recovery(cnp.ndarray[cnp.double_t, ndim=2] covf, \
                       cnp.ndarray[cnp.double_t, ndim=3] dcov, \
                       cnp.ndarray[cnp.double_t, ndim=1] covn_vec, \
                       cnp.ndarray[cnp.double_t, ndim=1] plist, \
                       int npt, int npix, int mdim_t, int mdim_f, int do_mpi):

    cdef: 
        int i, a, b, *idx_a, *idx_b

    idx_a=<int *>malloc(2*sizeof(int))
    idx_b=<int *>malloc(2*sizeof(int))

    do_mpi=myfalse

    # ->> now start <<- #
    if do_mpi==mytrue:
        _covf_=np.zeros((npix,npix))
        prange=mpi.mpirange(npix**2)

        for idx in prange:
            a=<int>(idx/<double>npix) 
            b=idx-npix*a

            mpixel_idx(a, mdim_t, mdim_f, idx_a)
            mpixel_idx(b, mdim_t, mdim_f, idx_b)

            for i in range(npt):
                _covf_[a,b]+=dcov[i,idx_a[0]-idx_b[0], idx_a[1]-idx_b[1]]*plist[i]/2.

            _covf_[a,a]+=covn_vec[a]

        print 'now gather full covariance matrix'
        covf=mpi.gather_unify(_covf_, root=0)
        #mpi.bcast(covf, rank=0)

        del(_covf_)

    else:

        # ->>  signal covariance matrix <<- #
        for a in range(npix):
            mpixel_idx(a, mdim_t, mdim_f, idx_a)

            for b in range(npix):
                mpixel_idx(b, mdim_t, mdim_f, idx_b)

                for i in range(npt):
                    covf[a,b]+=dcov[i,idx_a[0]-idx_b[0], idx_a[1]-idx_b[1]]*plist[i]/2.

            covf[a,a]+=covn_vec[a]


    free(idx_a)
    free(idx_b)

    return



#cpdef convert_cov_full(covf, dcov, covn_vec, plist, npt, npix, m_dim):
#
#    #print '--->>', covf.shape, dcov.shape, covn_vec.shape, plist.shape
#
#    full_cov_recovery(covf, dcov, covn_vec, plist, <int>npt, <int> npix, \
#                      <int> m_dim[0], <int> m_dim[1])
#
#    return





''' ->> Quadratic Estimator <<- #'''
cdef void icov_d_multiple(cnp.ndarray[cnp.double_t, ndim=2] icovf, \
                          cnp.ndarray[cnp.double_t, ndim=2] dmap,\
                          double *out, int npix, int mdim_t, int mdim_f):

    # ->> give the multiplication C^-1_{ab} dmap_{b} <<- #

    cdef: 
        int a, b, c, *idx_b

    idx_b=<int *>malloc(2*sizeof(int))

    for a in range(npix):
        out[a]=0.

        for b in range(npix):
            mpixel_idx(b, mdim_t, mdim_f, idx_b)
            out[a]+=icovf[a,b]*dmap[idx_b[0],idx_b[1]]

            #out[a]+=icovf[a,b]*dmap.flatten()[b]

    free(idx_b)

    return

cdef void icov_dov_multiple(cnp.ndarray[cnp.double_t, ndim=2] icovf,\
                            cnp.ndarray[cnp.double_t, ndim=3] dcov, \
                            cnp.ndarray[cnp.double_t, ndim=3] ic_dcov, \
                            int npt, int npix, int mdim_t, \
                            int mdim_f, int do_mpi):

    cdef:
        int a, b, c, i, *idx_b, *idx_c

    idx_b=<int *>malloc(2*sizeof(int))
    idx_c=<int *>malloc(2*sizeof(int))

    #print 'icov_dov start', mpi.rank, 'npt=', npt, 'npix=', npix

    _ic_dcov_=np.zeros((npt, npix, npix))

    if do_mpi==mytrue:
        prange=mpi.mpirange(npt)
    else:
        prange=range(npt)


    #for i in range(npt):
    for i in prange:
        print 'icov_dcov:', i, mpi.rank

        for a in range(npix):
            print 'icov_dcov i, a:', i, a, '(',mpi.rank,')'

            for b in range(npix):

                _ic_dcov_[i,a,b]=0.
                mpixel_idx(b, mdim_t, mdim_f, idx_b)

                for c in range(npix):
                    mpixel_idx(c, mdim_t, mdim_f, idx_c)
                    _ic_dcov_[i,a,b]+=icovf[a,c]*dcov[i,idx_c[0]-idx_b[0],idx_c[1]-idx_b[1]]

    if do_mpi==mytrue:
        ic_dcov=mpi.gather_unify(_ic_dcov_, root=0)
        #mpi.bcast(ic_dcov, rank=0)
    else:
        ic_dcov=np.copy(_ic_dcov_)


    free(idx_b)
    free(idx_c)

    del _ic_dcov_

    return


cdef void quad_estimator(cnp.ndarray[cnp.double_t, ndim=2] dmap, \
                        cnp.ndarray[cnp.double_t, ndim=2] covf, \
                        cnp.ndarray[cnp.double_t, ndim=3] dcov, \
                        cnp.ndarray[cnp.double_t, ndim=1] covn_vec, \
                        cnp.ndarray[cnp.double_t, ndim=1] plist, \
                        cnp.ndarray[cnp.double_t, ndim=1] Qi_p, \
                        cnp.ndarray[cnp.double_t, ndim=2] Fij, \
                        int npt, int npix, int mdim_t, int mdim_f, int do_mpi):

    cdef: 
        int i, a, b, *idx_a, *idx_b
        double *d_ic #, *ic_dcov

    d_ic=<double *>malloc(npix*sizeof(double))
    #ic_dcov=<double *>malloc(npt*npix*npix*sizeof(double))
    ic_dcov=np.zeros((npt, npix, npix))


    idx_a=<int *>malloc(2*sizeof(int))
    idx_b=<int *>malloc(2*sizeof(int))

    # ->> get full covariance matrix and its inverse <<- #
    print '->> now start to recover full covariance matrix.'
    full_cov_recovery(covf, dcov, covn_vec, plist, npt, npix, mdim_t, mdim_f, do_mpi)
    print '->> the inverse of covariance matrix.'
    icovf=slag.inv(covf)

    _testing_=False
    if _testing_:

        if mpi.rank0:
            fname='../../workspace/result/covariance.npz'
            np.savez(fname, covf=covf, icovf=icovf)

            print 'output the covariance matrix'

        mpi.finalize()
        quit()

    mpi.barrier()

    print '->> preparation of icov & cov is done, now calculate Qe.', mpi.rank

    # ->> get (C^{-1}.d) <<- #
    icov_d_multiple(icovf, dmap, d_ic, npix, mdim_t, mdim_f)
    print '->> icov.d is done.', mpi.rank

    # ->> get (C^{-1}.dcov[i]) <<- #
    icov_dov_multiple(icovf, dcov, ic_dcov, npt, npix, mdim_t, mdim_f, do_mpi)
    print '->> icov.dcov is done.', mpi.rank


    # ->> now start <<- #
    if do_mpi==mytrue:
        prange=mpi.mpirange(npt)
    else:
        prange=range(npt)


    for i in prange:
        print 'quadratic estimator: rank-', mpi.rank, '  prange:', i
	#->> for each i, calculate C^-1 dcov^i C^-1

        Qi_p[i]=0.
        Fij[i,j]=0.

        for a in range(npix):
            mpixel_idx(a, mdim_t, mdim_f, idx_a)

            for b in range(npix):
                mpixel_idx(b, mdim_t, mdim_f, idx_b)
                Qi_p[i]+=d_ic[a]*dcov[i,idx_a[0]-idx_b[0], idx_a[1]-idx_b[1]]*d_ic[b]/2.

        # ->> Fisher matrix <<- #
        for j in range(npt):

            for a in range(npix):
                for b in range(npix):
                    #Fij[i,j]=ic_dcov[matidx3(i,a,b,npt,npix,npix)]*ic_dcov[matidx3(j,b,a,npt,npix,npix)]/2.
                    Fij[i,j]=ic_dcov[i,a,b]*ic_dcov[j,b,a]/2.
	

        print 'i=', i, '(', mpi.rank, '),  Qi=', Qi_p[i]


    free(d_ic)
    #free(ic_dcov)
    free(idx_a)
    free(idx_b)

    del ic_dcov


    return



cdef void quad_est_fish_qi(cnp.ndarray[cnp.double_t, ndim=2] i_Fij, \
                           cnp.ndarray[cnp.double_t, ndim=1] Qi_p, \
                           cnp.ndarray[cnp.double_t, ndim=1] Qi, \
                           int npt, int npix, do_mpi=False):
    cdef:
        int i, j

    # ->> now start <<- #
    if do_mpi==True:
        prange=mpi.mpirange(npt)
    else:
        prange=range(npt)

    for i in prange:
        Qi[i]=0.
        for j in range(npt):
            Qi[i]+=i_Fij[i,j]*Qi_p[j]

    return



cpdef quad_estimator_wrapper(dmap, covf, dcov, covn_vec, plist, Qi, Fij, npt, npix, m_dim, do_mpi=False):

    cdef int dompi

    if do_mpi==True:
        dompi=mytrue
    else:
        dompi=myfalse


    _Qi_p_=np.zeros(npt)
    _Qi_=np.zeros(npt)
    _Fij_=np.zeros((npt, npt))

    quad_estimator(dmap, covf, dcov, covn_vec, plist, _Qi_p_, _Fij_, <int> npt, \
                   <int> npix, <int> m_dim[0], <int> m_dim[1], dompi)

    # ->> gather & broadcast <<- #
    if dompi:
        Qi_p=mpi.gather_unify(_Qi_p_, root=0)
        Fij=mpi.gather_unify(_Fij_, root=0)
        
        #mpi.bcast(Qi_p, rank=0)
        #mpi.bcast(Fij, rank=0)
    else:
        Qi_p=_Qi_p_
        Fij=np.copy(_Fij_)


    # ->> inverse Fij <<- #
    i_Fij=slag.inv(Fij)

    # ->> Qi=sum_j (F_ij Q_j) <<- #
    quad_est_fish_qi(i_Fij, Qi_p, _Qi_, npt, npix, dompi)

    if dompi:
        Qi=mpi.gather_unify(_Qi_, root=0)
        #mpi.bcast(Qi, rank=0)
    else:
        Qi=np.copy(_Qi_)


    del _Qi_p_, _Qi_, _Fij_
    print 'exiting quad_estimator_wrapper: rank-', mpi.rank

    return 





''' ->> some testing <<- '''

cdef void correlation_recovery(cnp.ndarray[cnp.double_t, ndim=2] covf, \
                             cnp.ndarray[cnp.double_t, ndim=3] dcov, \
                             cnp.ndarray[cnp.double_t, ndim=1] plist, \
                             int npt, int mdim_t, int mdim_f):
    cdef: 
        int i, a, b
    # ->> factor 1/2. appears becuase here we're integrating the 
    # ->> full k-space, however, dcov include both +/- k information

    # ->>  obtain correlation function matrix <<- #
    for a in range(mdim_t):
        for b in range(mdim_f):
            for i in range(npt):
                covf[a,b]+=dcov[i,a,b]*plist[i]/2.

    return 




cpdef get_correlation(covf,  dcov,  plist, npt):
    cdef int mdim_t, mdim_f

    mdim_t, mdim_f=covf.shape

    correlation_recovery(covf, dcov, plist, <int>npt, mdim_t, mdim_f)
    return



''' ->> send to Carly <<- '''
cdef return_full_dcov(cnp.ndarray[cnp.double_t, ndim=3] dcov, \
                      cnp.ndarray[cnp.double_t, ndim=2] dcov_f, \
                      int pi_idx, int npt, int mdim_t, int mdim_f):
    cdef:
        int a, b, *idx_a, *idx_b

    idx_a=<int *>malloc(2*sizeof(int))
    idx_b=<int *>malloc(2*sizeof(int))

    for a in range(npt):
        mpixel_idx(a, mdim_t, mdim_f, idx_a)

        for b in range(npt):
            #dcov_f[a,b]=0.

            mpixel_idx(b, mdim_t, mdim_f, idx_b)
            dcov_f[a,b]=dcov[pi_idx,abs(idx_a[0]-idx_b[0]), abs(idx_a[1]-idx_b[1])]

    free(idx_a)
    free(idx_b)
    
    return


cpdef get_full_dcov(dcov, dcov_f, pi_idx, npt, mdim_t, mdim_f):
    return return_full_dcov(dcov, dcov_f, <int>pi_idx, <int>npt, \
                            <int>mdim_t, <int>mdim_f)

  /*>> matrix opertion <<*/
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <iniparser.h>

  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "power.h"
  #include "cospara.h"
  #include "myinterpolate.h"

  #include "glbvarb.h"
  #include "mpinit.h"
  #include "misc.h"
  #include "io.h"
  #include "quadest.h"


#ifdef _MPI_
  #include <mpi.h>
#endif

#ifdef _OMP_
  #include <omp.h>
#endif




  double access_dcov(double *dcov, size_t n_bp, size_t npix, int i, 
                                      int a, int b, size_t map_dim) {
    // ->> i:    index of bandpower
    // ->> a/b:  indices of map pixels 

    if(map_dim==1)
      {return ArrayAccess2D_n2(dcov, n_bp, npix, i, (int)fabs((double)(a-b)));}

    else if(map_dim==2) {abort();}

    else {abort();}

    return -99.;
    }





  //void Fisher(MPIpar *mpi, double *dcov, double *icov, double *F,
  //                              size_t npix, size_t n_bp, size_t map_dim)  {

  double *Fisher(MPIpar *mpi, double *dcov, double *icov, size_t npix, 
                 size_t n_bp, size_t map_dim)  {

    int i, j, a, b, c, d, idx, id, irk, nrun;
    double *Fs, *Frev, *F;

    #ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    mpi->max=n_bp*n_bp;    mpi->start=0;
    mpi_loop_init(mpi, "Fisher");
    //mpi_loop_init(mpi, NULL);

    Fs=(double *)malloc(mpi->ind_run*sizeof(double));

    for(idx=0; idx<mpi->ind_run; idx++) {

      id=mpi_id(mpi, idx);
      i=(int)(id/(double)n_bp);
      j=id-i*n_bp;

      //#ifdef _OMP_
      //#pragma omp parallel for private(a,b,c,d)
      //#endif

      Fs[idx]=0.;
      for(a=0; a<npix; a++)
          for(b=0; b<npix; b++)
            for(c=0; c<npix; c++)
              for(d=0; d<npix; d++) {

                Fs[idx]+=access_dcov(dcov, n_bp, npix, i, a, b, map_dim)*
                         ArrayAccess2D(icov, npix, b, c)*access_dcov(dcov, n_bp, npix, 
			 j, c, d, map_dim)*ArrayAccess2D(icov, npix, d, a)/2.;
                }

        //printf("Fij[%d, %d]=%lg\n", i, j, Fs[idx]);
        //fflush(stdout);
      }
    printf("Fisher is done. (rank %d)\n", mpi->rank); fflush(stdout);


    /*
    if(mpi->rank==0){ Frev=(double *)malloc(sizeof(double)*n_bp*n_bp); }

    // ->> gather all data by root <<- //
    #ifdef _MPI_
    MPI_Gather(Fs, mpi->ind_run, MPI_DOUBLE, Frev, mpi->ind_run, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #else

    #endif

    if(mpi->rank==0){

      for(irk=0; irk<mpi->ntask; irk++) {
        nrun=mpi_nrun(mpi->max, irk, mpi->ntask);

        for(i=0; i<nrun; i++)
          F[mpi_get_id(irk, mpi->ntask, i)]=Frev[irk*nrun+i];
	  }

      printf("got here. (rank-%d)\n", mpi->rank); fflush(stdout);
      free(Frev);
      }

    #ifdef _MPI_
    MPI_Bcast(F, n_bp*n_bp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif
    */

    //mpi_gather_dist_double(mpi, Fs, F, mpi->ind_run, mpi->max);
    F=mpi_gather_dist_double(mpi, Fs, mpi->ind_run, mpi->max);

    printf("Existing Fisher. (rank-%d)\n", mpi->rank); fflush(stdout);

    char *fn="result/r1d/Fij_ins.dat";
    write_data(mpi, fn, Fs, sizeof(double), n_bp*n_bp);


    #ifndef _MPI_
    free(Fs);
    #endif

    return F;
    }



  void cov_noise(MPIpar *mpi, double *covn_v, size_t npix, char *type){
    int i;

    if(type==NULL) {
      for(i=0; i<npix; i++)
        covn_v[i]=0;
      }
    else abort();


    return;
    }


  //void full_covmat_recov(MPIpar *mpi, double *dcov, double *cov, double *covn_v, 
  //                       double *plist, size_t n_bp, size_t npix, size_t map_dim) {

  double *full_covmat_recov(MPIpar *mpi, double *dcov, double *covn_v, double *plist, 
                            size_t n_bp, size_t npix, size_t map_dim) {

    // ->> obtain full covariance matrix from dvoc and pk_list <<- //

    int ip, i, j, a, b, c, d, idx, id, irk, nrun;
    double *cov_s, *crev, *cov;
    char *fn;

    #ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    mpi->max=npix*npix;    mpi->start=0;
    mpi_loop_init(mpi, "cov");

    cov_s=(double *)malloc(mpi->ind_run*sizeof(double));

    for(idx=0; idx<mpi->ind_run; idx++) {

      // ->> pixel location <<- //
      id=mpi_id(mpi, idx);

      // ->> 1D map <<- //
      if(map_dim==1)  {

        // ->> convert to pixel index <<- //
        a=(int)(id/(double)npix);
        b=(int)(id-a*npix);

        cov_s[idx]=0.;
        for(ip=0; ip<n_bp; ip++){
          // ->> summing over all bandpowr <<- //
          cov_s[idx]+=access_dcov(dcov, n_bp, npix, ip, a, b, map_dim)*plist[ip];
          }

        //if(a==b){ cov_s[idx]+=covn_v[idx]; }
        }

      // ->> 2D map <<- //
      if(map_dim==2)  {
        abort();   // ->> DEFINITELY some error here. <<- //

        cov_s[idx]=0.;
        for(ip=0; ip<n_bp; ip++){
          // ->> summing over all bandpowr <<- //
          cov_s[idx]+=access_dcov(dcov, n_bp, npix, ip, a, b, map_dim)*plist[ip];
          }

        //if(a==b){ cov_s[idx]+=covn_v[idx]; }
        }

      }

    /*
    // ->> gather all data by root <<- //
    if(mpi->rank==0){ crev=(double *)malloc(sizeof(double)*npix*npix); }

    // ->> gather <<- //
    MPI_Gather(cov_s, mpi->ind_run, MPI_DOUBLE, crev, mpi->ind_run, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // ->> re-organize <<- //
    if(mpi->rank==0){

      for(irk=0; irk<mpi->ntask; irk++) {
        nrun=mpi_nrun(mpi->max, irk, mpi->ntask);

        for(i=0; i<nrun; i++)
          cov[mpi_get_id(irk, mpi->ntask, i)]=crev[irk*nrun+i];
	  }

      free(crev);
      }

    // ->> broadcast <<- //
    MPI_Bcast(cov, npix*npix, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    */


    cov=mpi_gather_dist_double(mpi, cov_s, mpi->ind_run, mpi->max);
    printf("pointer: %p  %p  (%d  %d)\n", cov, cov_s, mpi->ind_run, mpi->max);


    fn="result/r1d/cov_out_ss.dat";
    write_data(mpi, fn, cov_s, sizeof(double), npix*npix);

    free(cov_s);
    return cov;
    }





  void quad_est(MPIpar *mpi, QEpar *qe) {
    // ->> calculate quadratic estimator <<- //
    int a, b, i, j, idx, id;
    char *fn;

    // ->> first recover the full covariance matrix <<- //
    //full_covmat_recov(mpi, qe->dcov, qe->cov, qe->covn_v, qe->plist, 
    //                  qe->n_bp, qe->npix, qe->map_dim);

    qe->cov=full_covmat_recov(mpi, qe->dcov, qe->covn_v, qe->plist, 
                              qe->n_bp, qe->npix, qe->map_dim);

    printf("qe->cov pointer: %p\n", qe->cov);
    printf("Full covariance matrix done.(%d)\n", mpi->rank); fflush(stdout);

    fn="result/r1d/cov_out.dat";
    write_data(mpi, fn, qe->cov, sizeof(double), qe->npix*qe->npix);

    //MPI_Barrier(MPI_COMM_WORLD);
    //abort();

    // ->> inverse matrix <<- //
    mat_inv(mpi, qe->cov, qe->icov, qe->npix);
    fn="result/r1d/icov_out.dat";
    write_data(mpi, fn, qe->icov, sizeof(double), qe->npix*qe->npix);

    //fn="result/r1d/icov.dat";
    //import_data_double(mpi, fn, qe->icov, sizeof(double), qe->npix*qe->npix);
    //printf("inverse of covariance matrix done.\n"); fflush(stdout);

    /*

    // ->> calculate Fisher matrix and its inverse <<- //
    qe->Fij=Fisher(mpi, qe->dcov, qe->icov, qe->npix, qe->n_bp, qe->map_dim);
    fn="result/r1d/Fij.dat";
    write_data(mpi, fn, qe->Fij, sizeof(double), qe->n_bp*qe->n_bp);


    mat_inv(mpi, qe->Fij, qe->iFij, qe->npix);
    fn="result/r1d/inv_Fij.dat";
    write_data(mpi, fn, qe->iFij, sizeof(double), qe->n_bp*qe->n_bp);

    */

    printf("import inv_fisher\n"); fflush(stdout);
    fn="result/r1d/inv_Fij.dat";
    //fn="result/r1d/Carly/Fij_inv.dat";
    qe->iFij=(double *)malloc(sizeof(double)*qe->n_bp*qe->n_bp);
    import_data_double(mpi, fn, qe->iFij, sizeof(double), qe->n_bp*qe->n_bp);


    printf("(inv)-Fisher Matrix done.\n"); fflush(stdout);
  
    // ->> some pre-calculation of Qi <<- //
    // ->> (C^{-1}.d) <<- //
    double *d_ic=(double *)malloc(sizeof(double)*qe->npix);
    mat_vec_mult(qe->icov, qe->map, d_ic, qe->npix, qe->npix);


    #ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    // ->> calculate Qi' <<- //
    double *Qip_s, *Qi_s;
    Qip_s=(double *)malloc(mpi->ind_run*sizeof(double));
    Qi_s=(double *)malloc(mpi->ind_run*sizeof(double));

    mpi->max=qe->n_bp;    mpi->start=0;
    mpi_loop_init(mpi, "Qi");

    for(idx=0; idx<mpi->ind_run; idx++) {
      id=mpi_id(mpi, idx);

      Qip_s[idx]=0.;
      for(a=0; a<qe->npix; a++)
        for(b=0; b<qe->npix; b++) {
          Qip_s[idx]+=d_ic[a]*access_dcov(qe->dcov, qe->n_bp, qe->npix, id, a, b, qe->map_dim)*d_ic[b]/2.;
          }

      Qi_s[idx]=0;
      for(j=0; j<qe->n_bp; j++) {
        Qi_s[idx]+=Qip_s[idx]*ArrayAccess2D(qe->iFij, qe->n_bp, id, j);
        }
      }


    // ->> gather all data by root <<- //
    
    qe->Qip=(double *)malloc(qe->n_bp*sizeof(double));
    qe->Qi =(double *)malloc(qe->n_bp*sizeof(double));

    //mpi_gather_dist_double(mpi, Qip_s, qe->Qip, mpi->ind_run, mpi->max);
    //mpi_gather_dist_double(mpi, Qi_s, qe->Qi, mpi->ind_run, mpi->max);

    qe->Qip=mpi_gather_dist_double(mpi, Qip_s, mpi->ind_run, mpi->max);
    qe->Qi=mpi_gather_dist_double(mpi, Qi_s, mpi->ind_run, mpi->max);

    free(Qip_s); free(Qi_s); free(d_ic);
    return;
    }





  /*>> matrix opertion <<*/
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
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


#ifdef _MPI_
  #include <mpi.h>
#endif

#ifdef _OMP_
  #include <omp.h>
#endif




  double access_dcov(double *dcov, size_t n_bp, size_t npix, size_t i, 
                                      size_t a, size_t b, int map_dim) {
    // ->> i:    index of bandpower
    // ->> a/b:  indices of map pixels 

    if(map_dim==1)
      {return ArrayAccess2D_n2(dcov, n_bp, npix, i, (size_t)fabs(a-b));}

    else if(map_dim==2) {abort();}

    else {abort();}
    }





  void Fisher(MPIpar *mpi, double *dcov, double *icov, double *F,
                                size_t npix, size_t n_bp, int map_dim)  {
    size_t i, j, a, b, c, d, idx, id, irk, nrun;
    double *Fs, *Frev;


    mpi->max=n_bp*n_bp;    mpi->start=0;
    mpi_loop_init(mpi, "Fisher");

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

                Fs[idx]+=access_dcov(dcov, n_bp, npix, i, a, b, map_dim)*icov[b,c]
		        *access_dcov(dcov, n_bp, npix, j, c, d, map_dim)*icov[d,a];
                }

        //printf("Fij[%d, %d]=%lg\n", i, j, Fs[idx]);
        //fflush(stdout);
      }
    printf("Fisher is done. (rank %d)\n", mpi->rank); fflush(stdout);


    if(mpi->rank==0){ Frev=(double *)malloc(sizeof(double)*n_bp*n_bp); }

    // ->> gather all data by root <<- //
    MPI_Gather(Fs, mpi->ind_run, MPI_DOUBLE, Frev, mpi->ind_run, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(mpi->rank==0){

      for(irk=0; irk<mpi->ntask; irk++) {
        nrun=mpi_nrun(mpi->max, irk, mpi->ntask);

        for(i=0; i<nrun; i++)
          F[mpi_get_id(irk, mpi->ntask, i)]=Frev[irk*nrun+i];
	  }

      free(Frev);
      }

    free(Fs);
    return;
    }



  void full_covmat_recov(MPIpar *mpi, double *dcov, double *cov, double *covn_v, 
                         double *plist, size_t n_bp, size_t npix, int map_dim) {
    // ->> obtain full covariance matrix from dvoc and pk_list <<- //

    size_t ip, i, j, a, b, c, d, idx, id, irk, nrun;
    double *cov_s, *crev;


    mpi->max=npix*npix;    mpi->start=0;
    mpi_loop_init(mpi, "covm");

    cov_s=(double *)malloc(mpi->ind_run*sizeof(double));

    for(idx=0; idx<mpi->ind_run; idx++) {

      id=mpi_id(mpi, idx);
      a=(int)(id/(double)npix);
      b=id-a*npix;

      cov_s[idx]=0.;

        for(ip=0; ip<n_bp; ip++){
	  // ->> summing over all bandpowr <<- //
          cov_s[idx]+=access_dcov(dcov, n_bp, npix, ip, a, b, map_dim)*plist[i];
	  }

      if(a==b){ cov_s[idx]+=covn_vec[idx]; }

      }


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


    free(cov_s);
    return;
    }





  void quad_est(MPIpar *mpi, QEpar *qe) {
    // ->> calculate quadratic estimator <<- //


    // ->> first recover the full covariance matrix <<- //
    full_covmat_recov(mpi, qe->dcov, qe->cov, qe->covn_v, qe->plist, 
                      qe->n_bp, qe->npix, qe->map_dim);

    return;
    }





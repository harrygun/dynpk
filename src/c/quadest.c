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
    mpi_loop_init(mpi, "output");

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

                //Fs[idx]+=ArrayAccess2D_n2(dcov, n_bp, npix, i, (size_t)fabs(a-b))*icov[b,c]
		//       *ArrayAccess2D_n2(dcov, n_bp, npix, j, (size_t)fabs(c-d))*icov[d,a];

                Fs[idx]+=access_dcov(dcov, n_bp, npix, i, a, b, map_dim)*icov[b,c]
		        *access_dcov(dcov, n_bp, npix, j, c, d, map_dim)*icov[d,a];
                }

        //printf("Fij[%d, %d]=%lg\n", i, j, Fs[idx]);
        //fflush(stdout);
      }
    printf("Fisher is done. %d\n"); fflush(stdout);


    if(mpi->rank==0){ Frev=(double *)malloc(sizeof(double)*n_bp*n_bp); }

    // ->> gather all data by root <<- //
    MPI_Gather(Fs, mpi->ind_run, MPI_DOUBLE, Frev, mpi->ind_run, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(mpi->rank==0){
      printf("re-store data.\n"); fflush(stdout);

      for(irk=0; irk<mpi->ntask; irk++) {
        nrun=mpi_nrun(mpi->max, irk, mpi->ntask);

        for(i=0; i<nrun; i++)
          F[mpi_get_id(irk, mpi->ntask, i)]=Frev[irk*nrun+i];
	  }

      printf("re-store data finished.\n"); fflush(stdout);
      free(Frev);
      }

    if(mpi->rank==0){ printf("before existing.\n"); fflush(stdout); }

    free(Fs);
    return;
    }

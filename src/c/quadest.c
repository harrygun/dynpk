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









  void Fisher(MPIpar *mpi, double *dcov, double *icov, double *F,
                                         size_t npix, size_t n_bp)  {
    size_t i, j, a, b, c, d, idx;


    mpi->max=n_bp*n_bp;    mpi->start=0;
    mpi_loop_init(mpi, "output");


    for(idx=0; idx<mpi->ind_run; idx++) {

      i=(int)(idx/(double)n_bp);
      j=idx-i*n_bp;

      //#ifdef _OMP_
      //#pragma omp parallel for private(a,b,c,d)
      //#endif

      //abort("need to optimize the parallization here.");

      F[idx]=0.;
      for(a=0; a<npix; a++)
          for(b=0; b<npix; b++)
            for(c=0; c<npix; c++)
              for(d=0; d<npix; d++) {

                F[idx]+=ArrayAccess2D_n2(dcov, n_bp, npix, i, (size_t)fabs(a-b))*icov[b,c]
		       *ArrayAccess2D_n2(dcov, n_bp, npix, j, (size_t)fabs(c-d))*icov[d,a];
                }

      printf("Fij[%d, %d]=%lg\n", i, j, F[idx]);
      fflush(stdout);
      }


    return;
    }

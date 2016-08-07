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


#ifdef _MPI_
  #include "mpi.h"
#endif


  void mpi_loop_init(MPIpar *mpi, char *prefix) {

  #ifdef _MPI_
    asprintf(&mpi->fname, "%s_%d.dat", prefix, mpi->rank);
    //printf("proc(%d): %s\n", mpi->rank, mpi->fname); 
    //fflush(stdout);

    /* total run needed */
    mpi->totrun=mpi->max-mpi->start;
    /* run 'mpi->ind_run' times for each individual process. */
    mpi->ind_run= mpi->totrun/mpi->ntask;
    if( mpi->rank - (mpi->totrun-mpi->ntask*(int)(mpi->totrun/mpi->ntask) )<0 ) 
      mpi->ind_run +=1;

    printf("proc(%d): nrun=%d, totrun=%d\n", mpi->rank, mpi->ind_run, mpi->totrun);
    fflush(stdout);

    /*if(mpi_rank==0) {
      // ???? what's this for??
      } */

   #else 
     mpi->ind_run = mpi->max-mpi->start; 
     asprintf(&mpi->fname, "%s.dat", prefix);
   #endif

    return;
  }


  int mpi_id(MPIpar *mpi, int i) {
    return (i*mpi->ntask + mpi->rank + mpi->start);
  }




  int mpi_nrun(int totrun, int rank, int ntask){
    int ind_run;

    ind_run= totrun/ntask;
    if(rank - (totrun-ntask*(int)(totrun/ntask) )<0 ) 
      ind_run +=1;

    return ind_run;
    }


  int mpi_get_id(int rank, int ntask, int i) {

    return i*ntask+rank;
    }

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
  #include "mpi.h"
#endif


  void mpi_gather_dist_double(MPIpar *mpi, double *in, double *out, 
                                       size_t count_pp, size_t count_tot ) {
    // ->> gather & redistribution:  count_pp: count per process 
    
    size_t irk, i, nrun;
    double *rev;
    
    if(mpi->rank==0)
      rev=(double *)malloc(sizeof(double)*count_tot); 

    // ->> gather <<- //
    MPI_Gather(in, count_pp, MPI_DOUBLE, rev, count_pp, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    printf("MPI_gather (rank=%d)\n", mpi->rank); fflush(stdout);

    // ->> re-organize <<- //
    if(mpi->rank==0){

      printf("Before calculating... \n"); fflush(stdout);

      for(irk=0; irk<mpi->ntask; irk++) {

        printf("irk=%d, ", irk); fflush(stdout);

        nrun=mpi_nrun(count_tot, irk, mpi->ntask);

	printf("%d\n", nrun); fflush(stdout);

        for(i=0; i<nrun; i++)
          out[mpi_get_id(irk, mpi->ntask, i)]=rev[irk*nrun+i];
	  }

      printf("MPI_gather before free.\n"); fflush(stdout);

    free(rev);
    }

    // ->> broadcast <<- //
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(out, count_tot, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return;
    }


  void mpi_loop_init(MPIpar *mpi, char *prefix) {

  #ifdef _MPI_
    //asprintf(&mpi->fname, "%s_%d.dat", prefix, mpi->rank);
    //printf("proc(%d): %s\n", mpi->rank, mpi->fname); 
    //fflush(stdout);

    /* total run needed */
    mpi->totrun=mpi->max-mpi->start;
    /* run 'mpi->ind_run' times for each individual process. */
    mpi->ind_run= (size_t)(mpi->totrun/mpi->ntask);
    if( mpi->rank - (mpi->totrun-mpi->ntask*(size_t)(mpi->totrun/mpi->ntask) )<0 ) 
      mpi->ind_run +=1;

    printf("proc(%d): nrun=%d, totrun=%d\n", mpi->rank, mpi->ind_run, mpi->totrun);
    fflush(stdout);

    /*if(mpi_rank==0) {
      // ???? what's this for??
      } */

   #else 
     mpi->ind_run = mpi->max-mpi->start; 
     //asprintf(&mpi->fname, "%s.dat", prefix);
   #endif

    return;
  }


  size_t mpi_id(MPIpar *mpi, size_t i) {
    return (i*mpi->ntask + mpi->rank + mpi->start);
  }




  size_t mpi_nrun(size_t totrun, size_t rank, size_t ntask){
    size_t ind_run;

    ind_run= (size_t)(totrun/ntask);
    if(rank - (totrun-ntask*(size_t)(totrun/ntask) )<0 ) 
      ind_run +=1;

    return ind_run;
    }


  size_t mpi_get_id(size_t rank, size_t ntask, size_t i) {

    return i*ntask+rank;
    }

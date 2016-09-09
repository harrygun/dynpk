  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>

  #include "glbvarb.hpp"
  #include "io.hpp"
  #include "mpinit.hpp"

  using namespace std;

  #ifdef _MPI_
  #include "mpi.h"
  #endif




  void mpi_init(MPIpar &mpi, int argc, char *argv[]){

    #ifdef _MPI_
      //int mpi_ntask, mpi_rank, mpi_rc;
      mpi.rc = MPI_Init(&argc, &argv);

      if (mpi.rc != MPI_SUCCESS) {
           cout << "Error starting MPI program. Terminating." << endl; 
           MPI_Abort(MPI_COMM_WORLD, mpi.rc);
           }

      MPI_Comm_size(MPI_COMM_WORLD,&mpi.ntask);
      MPI_Comm_rank(MPI_COMM_WORLD,&mpi.rank);

      printf ("Number of tasks= %d My rank= %d\n", mpi.ntask, mpi.rank);

  /* MPI initialization ends. */
      if(mpi.rank==0) 
        printf("%d Sending parameter filename %s to other processes.\n", mpi.rank, ini_name);

      MPI_Bcast( ini_name, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

      if(mpi.rank!=0) 
        printf("%d Received parameter filename %s.\n", mpi.rank, ini_name);
    #else
      mpi.ntask = 1;
      mpi.rank = 0;
    #endif

    return;
    }


  double *mpi_gather_dist_double(MPIpar *mpi, double *in, size_t count_pp, 
                                 size_t count_tot ) {
    // ->> gather & redistribution:  count_pp: count per process 
    
    int irk, i, nrun;
    double *rev, *out;

    #ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);

    if(mpi->ntask==1){return in;}
    
    out=(double *)malloc(sizeof(double)*count_tot); 

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

    #else
    // ->> no MPI <<- //
    out=in;
    #endif

    return out;
    }






  void mpi_loop_init(MPIpar *mpi, char *prefix) {

  #ifdef _MPI_
    //asprintf(&mpi->fname, "%s_%d.dat", prefix, mpi->rank);
    //printf("proc(%d): %s\n", mpi->rank, mpi->fname); 
    //fflush(stdout);

    /* total run needed */
    mpi->totrun=mpi->max-mpi->start;
    /* run 'mpi->ind_run' times for each individual process. */
    mpi->ind_run= (int)(mpi->totrun/mpi->ntask);
    if( mpi->rank - (mpi->totrun-mpi->ntask*(int)(mpi->totrun/mpi->ntask) )<0 ) 
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


  int mpi_id(MPIpar *mpi, int i) {
    return (i*mpi->ntask + mpi->rank + mpi->start);
  }




  int mpi_nrun(int totrun, int rank, int ntask){
    int ind_run;

    ind_run= (int)(totrun/ntask);
    if(rank - (totrun-ntask*(int)(totrun/ntask) )<0 ) 
      ind_run +=1;

    return ind_run;
    }


  int mpi_get_id(int rank, int ntask, int i) {

    return i*ntask+rank;
    }

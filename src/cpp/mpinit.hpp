  #ifndef _H_MPI_INIT_
  #define _H_MPI_INIT_

  #ifdef _MPI_
    #include <mpi.h>
  #endif
  




  //void mpi_init(MPIpar &glmpi);




  double *mpi_gather_dist_double(MPIpar *mpi, double *in, size_t count_pp, 
                                 size_t count_tot );

  void mpi_loop_init(MPIpar *mpi, char *prefix);


  int mpi_id(MPIpar *mpi, int i);


  int mpi_nrun(int totrun, int rank, int ntask);
  int mpi_get_id(int rank, int ntask, int i);


  #endif

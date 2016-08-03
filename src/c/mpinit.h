  #ifndef _H_MPI_INIT_
  #define _H_MPI_INIT_

  void mpi_loop_init(MPIpar *mpi, char *prefix);


  int mpi_id(MPIpar *mpi, int i);


  int mpi_nrun(int totrun, int rank, int ntask);
  int mpi_get_id(int rank, int ntask, int i);

  #endif

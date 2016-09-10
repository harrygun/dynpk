  #ifndef _H_MPI_INIT_
  #define _H_MPI_INIT_


    #include "El.hpp"

    class MPIpar {

      public:
        int ntask, rank, rc, root;
        int max, start, totrun, ind_run;
	bool rank0;

        std::string fname, extfname[5];
        El::mpi::Comm world;

        MPIpar () {

          ntask=El::mpi::Size(world);
          rank=El::mpi::Rank(world);

          std::cout << "Number of tasks= " << ntask << "; My rank= " << rank << std::endl;

          if(rank==0)  {rank0=true;}
	  else {rank0=false;}

          }


      };





    double *mpi_gather_dist_double(MPIpar *mpi, double *in, size_t count_pp, 
                                   size_t count_tot );

    void mpi_loop_init(MPIpar *mpi, char *prefix);


    int mpi_id(MPIpar *mpi, int i);


    int mpi_nrun(int totrun, int rank, int ntask);
    int mpi_get_id(int rank, int ntask, int i);


  #endif

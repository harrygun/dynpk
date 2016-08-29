  #include <iostream>
  #include <cmath>
  #include <cstring>
  using namespace std;


  #include "glbvarb.hpp"
  //#include "mpinit.h"


#ifdef _MPI_
  #include <mpi.h>
#endif

#ifdef _OMP_
  #include <omp.h>
#endif




  void import_data_double(MPIpar *mpi, const char *fn, void *d, size_t size, size_t count) {
    FILE *fp;

    #ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    if(mpi->rank==0){

      if(!(fp=fopen(fn, "r"))) {
        printf("can't open file `%s`\n", fn); fflush(stdout);
        exit(0);
        }

      if(!(fread(d, size, count, fp)) ) {
        printf("File '%s' import Error.\n", fn); fflush(stdout);
        exit(0);
        }

      fclose(fp);
      }
    
    #ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(d, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif

    return;
    }


  void write_data(MPIpar *mpi, const char *fn, const void *d, size_t size, size_t count){
    FILE *fp;

    //printf("write_data: %p\n", d);

    #ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    if(mpi->rank==0){
      if(!(fp=fopen(fn, "wb"))) {
        printf("can't open file `%s`\n", fn); fflush(stdout);
        exit(0);
        }

      fwrite(d, size, count, fp);

      fclose(fp);
      }

    return;
    }

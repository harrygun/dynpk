  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  using namespace std;


  #include "glbvarb.hpp"
  #include "mpinit.hpp"




  void import_data_double(MPIpar *mpi, const char *fn, void *d, size_t size, size_t count) {
    FILE *fp;

    #ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Barrier(mpi->world);
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

    //MPI_Barrier(mpi->world);
    //MPI_Bcast(d, count, MPI_DOUBLE, 0, mpi->world);
    #endif

    return;
    }


 /*
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
    */

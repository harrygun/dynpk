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




  void import_data_double(MPIpar *mpi, char *fn, void *d, size_t size, size_t count) {
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


  void write_data(MPIpar *mpi, char *fn, const void *d, size_t size, size_t count){
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

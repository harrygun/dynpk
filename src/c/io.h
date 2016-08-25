  #ifndef _H_IO_
  #define _H_IO_

  #include "glbvarb.h"

  //void import_data(char *fn, void *d, size_t size, size_t count);
  void import_data_double(MPIpar *mpi, char *fn, void *d, size_t size, size_t count);

  //void write_data(char *fn, void *d, size_t size, size_t count);
  //void write_data(MPIpar *mpi, char *fn, void *d, size_t size, size_t count);
  void write_data(MPIpar *mpi, char *fn, const void *d, size_t size, size_t count);

  #endif

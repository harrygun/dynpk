  #ifndef _H_IO_
  #define _H_IO_

  #include "glbvarb.hpp"

  void import_data_double(MPIpar *mpi, const char *fn, void *d, size_t size, size_t count);

  void write_data(MPIpar *mpi, const char *fn, const void *d, size_t size, size_t count);

  #endif

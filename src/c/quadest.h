#ifndef _H_QUADEST_
#define _H_QUADEST_

  void Fisher(MPIpar *mpi, double *dcov, double *icov, double *F,
                                size_t npix, size_t n_bp, int map_dim);
#endif

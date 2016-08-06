#ifndef _H_QUADEST_
#define _H_QUADEST_

  double access_dcov(double *dcov, size_t n_bp, size_t npix, size_t i, 
                                      size_t a, size_t b, int map_dim);

  void Fisher(MPIpar *mpi, double *dcov, double *icov, double *F,
                                size_t npix, size_t n_bp, int map_dim);
#endif

#ifndef _H_QUADEST_
#define _H_QUADEST_

  #include "glbvarb.h"

  double access_dcov(double *dcov, size_t n_bp, size_t npix, int i, 
                                      int a, int b, size_t map_dim);

  void cov_noise(MPIpar *mpi, double *covn_v, size_t npix, char *type);

  //void Fisher(MPIpar *mpi, double *dcov, double *icov, double *F,
  //                              size_t npix, size_t n_bp, int map_dim);
  double *Fisher(MPIpar *mpi, double *dcov, double *icov, size_t npix, 
                 size_t n_bp, size_t map_dim);

  double *full_covmat_recov(MPIpar *mpi, double *dcov, double *covn_v, double *plist, 
                            size_t n_bp, size_t npix, size_t map_dim);


  void quad_est(MPIpar *mpi, QEpar *qe);

#endif

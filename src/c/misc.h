
  #ifndef _H_MISC_
  #define _H_MISC_


  void mat_inv(MPIpar *mpi, double *matx, double *matx_inv, int n);


  void mat_vec_mult(double *m, double *v, double *out, size_t nm, size_t nv);

  #endif

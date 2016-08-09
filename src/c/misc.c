  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <iniparser.h>
  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_linalg.h>
  #include <gsl/gsl_blas.h>


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



  
  void mat_inv(MPIpar *mpi, double *matx, double *matx_inv, size_t n) {

    size_t i, j, k; 
    int signum;
    double ele;
    
    gsl_permutation *perm=gsl_permutation_alloc(n);
    gsl_matrix_view mat = gsl_matrix_view_array (matx, n, n);

    
    gsl_vector *vec_e=gsl_vector_alloc(n);;
    gsl_vector_view vec_inv;


    // ->> parallized matrix inversion <<- //
    //mpi->max=n_bp*n_bp;    mpi->start=0;
    //mpi_loop_init(mpi, "Fisher");
    // ->> !!! LET's parallize the code later !!! <<- //
    
    // ->> LU decomposition <<- //
    gsl_linalg_LU_decomp(&mat.matrix, perm, &signum);


    // ->> get inverse <<- //
    for(i=0; i<n; i++) {

      vec_inv=gsl_vector_view_array( &(matx_inv[n*i]), n);

      for(j=0; j<n; j++) {
        if(j==i)
          gsl_vector_set(vec_e, j, 1);
        else
          gsl_vector_set(vec_e, j, 0);
        }
    
      gsl_linalg_LU_solve(&mat.matrix, perm, vec_e, &vec_inv.vector);
      }
    

    /*
    for(j=0;j<n;j++)
      for(k=0;k<n;k++) {
        ele=gsl_vector_get(vec_inv[k],j);
        matx_inv[j][k] = ele;
        }
    */
    
    

    gsl_permutation_free(perm);
    //gsl_matrix_free(mat_dest);
    //gsl_matrix_free(mat);
    
    gsl_vector_free(vec_e);

    //for(j=0;j<n;j++) {
    //  gsl_vector_free(vec_inv[j]);
    //  }
    

    return;
    }
  



  void mat_vec_mult(double *m, double *v, double *out, size_t nm, size_t nv) {
    // ->> m:  (nm, nv);   v:  (nv);    out: (nm)
    // ->> could be easily parallized by OpenMP <<- #

    size_t i, j;

    for(i=0; i<nm; i++) {
      out[i]=0.;

      for(j=0; j<nv; j++) 
        out[i]+=ArrayAccess2D_n2(m, nm, nv, i, j)*v[j];

      }

    return;
    }

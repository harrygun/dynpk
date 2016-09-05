
  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  #include "El.hpp"

  #include "glbvarb.hpp"
  #include "io.hpp"
  #include "mpinit.hpp"
  #include "quadest.hpp"

  using namespace std;
  using namespace El;
  typedef double Real;


  #ifdef _MPI_
  #include <mpi.h>
  #endif



  void DistMatrix_init(double *d, int pidx) {
    
    DistMatrix<double> dcov(nrows,ncols);

    return;
    }

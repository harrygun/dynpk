
  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  #include <string>

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



  void QE_init(QEpar &qe, const string &mat_type="DistMatrix") {
     
    if (mat_type.compare("DistMatrix") != 0)
      throw runtime_error("Error: Only DistMatrix type is Supported.");

    DistMatrix<double> qe.dcov_vec(qe.n_bp, qe.npix), qe.dcov_i(qe.npix, qe.npix), 
                       qe.cov(qe.npix, qe.npix),   qe.icov(qe.npix, qe.npix);

    return;
    }




  void dcov_recovery( ) {


    return;
    }



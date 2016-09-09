
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
  //#include "quadest.hpp"

  using namespace std;
  using namespace El;

  typedef double Real;

  #ifdef _MPI_
  #include <mpi.h>
  #endif



  void QE_init(QEpar &qe, const string &mat_type="DistMatrix") {
     
    if (mat_type.compare("DistMatrix") != 0)
      throw runtime_error("Error: Only DistMatrix type is Supported.");

    //DistMatrix<double> qe.dcov_vec(qe.nbp, qe.npix), qe.dcov_i(qe.npix, qe.npix), 
    //                   qe.cov(qe.npix, qe.npix),   qe.icov(qe.npix, qe.npix);

    //Read(qe.dcov_vec, fname); 

    return;
    }



  void dcov_recovery(QEpar &qe, DistMatrix<double> dcov[], DistMatrix<double> dcov_vec) {
    //->>
    vector <double> colvectorsymm;

    for(int i=0; i<qe.nbp; i++) {
      dcov[i]=DistMatrix<double>(qe.npix,qe.npix);
      Toeplitz(dcov[i], qe.npix, qe.npix, colvectorsymm);
      }

    return;
    }




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


  double SinIntegral(double x){
    cdef double si;
    si=spec.sici(x)[0];

    if (np.isnan(si))|(np.isinf(si)):
        print 'SineIntegral nan/inf:', x
    return si;
    }

  double CosIntegral(double x){
    cdef double ci
    ci=spec.sici(x)[1]

    if (np.isnan(ci))|(np.isinf(ci)):
        print 'CosIntegral nan/inf:', x
    return ci
  }


  double Cos(double x){
    return np.cos(x);
    }


  double Sin(double x){
    return np.sin(x)
    }









  void dcov_vector(){

    }




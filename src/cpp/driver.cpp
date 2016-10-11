  /*--------------------------------------------------------------------------
          >>>         Driver for Quadratic Estimator       <<<
  --------------------------------------------------------------------------*/
  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  #include <El.hpp>
  #include <boost/property_tree/ptree.hpp>
  #include <boost/property_tree/ini_parser.hpp>

  #include "glbvarb.hpp"
  #include "io.hpp"
  #include "mpinit.hpp"
  #include "quadest.hpp"

  using namespace std;
  using namespace El;



  int main(int argc, char *argv[])  {

    char *ini_name;
    ini_name=argv[1];

    Environment env( argc, argv );
    ProcessInput();

    /*------------------------------------------------
               Global MPI initialization
      ------------------------------------------------*/
    MPIpar glmpi;

    if(glmpi.rank==0)
      cout << "Opening File:  " << ini_name << endl;

    /*-----------------------------------------------
         ->>   Parameters Initialization   <<- 
      -----------------------------------------------*/
    string sec="Quadratic_Estimator_Boost";
    QEpar qe(ini_name, sec, &glmpi);

    if(glmpi.rank0)
      cout << "rank-" << glmpi.rank << "; nbp=" << qe.nbp 
           << "; npix=" << qe.npix << endl;



    string fn_Qi;
    //fn_dcov="result/r1d/dcov.dat";
    //fn_map ="result/r1d/dmap.dat";
    //fn_plist="result/r1d/plist.dat";

    fn_out ="result/r1d/Qi.dat";




    // ->> call Quadratic Estimator <<- //
    int n_it=1;
    qe.Quad_Estimator(qe.pfid, n_it);


    // output //
    Write(qe.Qi, fn_Qi);


    return 0;
    }





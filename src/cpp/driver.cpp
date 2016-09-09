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
  //#include "quadest.hpp"
  #include "quadest_init.hpp"

  using namespace std;
  using namespace El;

  #ifdef _MPI_
  #include <mpi.h>
  #endif








  int main( int argc, char *argv[])  {

      char *ini_name;
      Environment env( argc, argv );
      ProcessInput();

      ini_name=argv[1];

      boost::property_tree::ptree pt;
      boost::property_tree::ini_parser::read_ini(ini_name, pt);
      // examples //
      //cout << pt.get<string>("Section1.Value1") << endl;
      //cout << pt.get<string>("Section1.Value2") << endl;
  
    /*------------------------------------------------
                MPI initialization.
    ------------------------------------------------*/
      MPIpar mpi;

  /*--------------------------------------
         End of MPI initialization.
  --------------------------------------*/
      cout << "Opening File:  " << ini_name << endl;


  /*-----------------------------------------------
        initialize the QE parameters
  -----------------------------------------------*/
      QEpar qe;
      /*
      char *output_prefix, *klam_fname, *plin_name;
      qe.= iniparser_getstring(dict, "General:", NULL);
      output_prefix=iniparser_getstring(dict,"General:output_prefix", NULL);
      qe.= iniparser_getint(dict, "General:", 0);
      qe.= iniparser_getdouble(dict, "General:", 1);
      */


    /*-----------------------------------------------
               Here begin the calculation.
    -----------------------------------------------*/

      /* Calculating Eulerian Biasing Model */
      mpi.start = 0;   // mpi.max = 1;
      cout << "==================================" << endl; fflush(stdout);

      const char *fn_dcov, *fn_cov, *fn_icov, *fn_out, *fn_map, *fn_plist;

      fn_dcov="result/r1d/dcov.dat";
      fn_map ="result/r1d/dmap.dat";
      fn_plist="result/r1d/plist.dat";
      fn_out ="result/r1d/Qi.dat";


      // ->>   <<- //
      //size_t mdim, npix, nbp;
      qe.mdim=50; 
      qe.npix=qe.mdim; // 1D

      qe.nbp=24;   //mdim*mdim;
      qe.map_dim=1;


      // ->> initialization <<- //
      DistMatrix<double> dcov[qe.nbp], cov(qe.npix,qe.npix), icov(qe.npix,qe.npix);



      QE_init(qe);
      abort();

      // ->> call Quadratic Estimator <<- //
      //quad_est_dismat(&mpi, &qe);

   
      // ->> output <<- //
      //printf("Output data.\n"); fflush(stdout);
      //write_data(&mpi, fn_out, qe.Qi, sizeof(double), qe.n_bp);


  /*-----------------------------------------------------
                     free all
  -----------------------------------------------------*/
    stop:

      //free(qe.dcov);   free(qe.icov); 
      //free(qe.cov);    free(qe.map);
      //free(qe.plist);  free(qe.covn_v);


    return 0;
    }





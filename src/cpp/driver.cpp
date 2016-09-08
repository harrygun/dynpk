  /*--------------------------------------------------------------------------
          >>>         Driver for Quadratic Estimator       <<<
  --------------------------------------------------------------------------*/
  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  //#include <El.hpp>

  //using namespace El;

  #include "glbvarb.hpp"
  #include "io.hpp"
  #include "mpinit.hpp"
  //#include "quadest.hpp"
  #include "quadest_init.hpp"

  using namespace std;

  #ifdef _MPI_
  #include <mpi.h>
  #endif



  int main( int argc, char *argv[])  {

      char *ini_name;
      Environment env( argc, argv );
      ProcessInput();

      //ini_name=argv[1];
  
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

      qe.n_bp=24;   //mdim*mdim;
      qe.map_dim=1;


      //qe.dcov=(double *)malloc(sizeof(double)*qe.n_bp*qe.npix);
      //qe.cov= (double *)malloc(sizeof(double)*qe.npix*qe.npix);
      //qe.icov=(double *)malloc(sizeof(double)*qe.npix*qe.npix);
      //qe.plist=(double *)malloc(sizeof(double)*qe.n_bp);
      //qe.map=(double *)malloc(sizeof(double)*qe.npix);



      //import_data_double(&mpi, fn_map, qe.map, sizeof(double), qe.npix);
      //import_data_double(&mpi, fn_dcov, qe.dcov, sizeof(double), qe.n_bp*qe.npix);
      //import_data_double(&mpi, fn_plist, qe.plist, sizeof(double), qe.n_bp);

      // ->> noise covariance matrix <<- //
      //qe.covn_v=(double *)malloc(sizeof(double)*qe.npix);
      //cov_noise(&mpi, qe.covn_v, qe.npix, NULL);




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





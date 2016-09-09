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

      // ->> Initialization <<- //
      ini_name=argv[1];


    /*------------------------------------------------
                MPI initialization.
    ------------------------------------------------*/
      MPIpar mpiwd;
      mpi_init(mpiwd, argc, argv);


    /*-----------------------------------------------
             Parameters Initialization 
    -----------------------------------------------*/
      cout << "Opening File:  " << ini_name << endl;
      boost::property_tree::ptree pt;
      boost::property_tree::ini_parser::read_ini(ini_name, pt);


      // examples //
      std::string sec="Quadratic_Estimator";
  

      QEpar qe;
      if(mpiwd.rank==0){
        //char *output_prefix, *klam_fname, *plin_name;
        //output_prefix=iniparser_getstring(dict,"General:output_prefix", NULL);

        qe.mdim=pt.get<double>(sec+".map_resolution_val");
	qe.nbp=pt.get<double>(sec+".num_band_power");
	qe.npix=qe.mdim;
        qe.map_dim=1;
        }

      cout << "mdim=" << qe.mdim << "; nbp=" qe.nbp << "; npix=" << qe.npix +1 << endl;

    /*-----------------------------------------------
               Here begin the calculation.
    -----------------------------------------------*/

      /* Calculating Eulerian Biasing Model */
      cout << "==================================" << endl; fflush(stdout);

      const char *fn_dcov, *fn_cov, *fn_icov, *fn_out, *fn_map, *fn_plist;

      fn_dcov="result/r1d/dcov.dat";
      fn_map ="result/r1d/dmap.dat";
      fn_plist="result/r1d/plist.dat";
      fn_out ="result/r1d/Qi.dat";


      // ->>   <<- //
      //size_t mdim, npix, nbp;
      //qe.mdim=50; 
      //qe.npix=qe.mdim; // 1D

      //qe.nbp=24;   //mdim*mdim;


      // ->> initialization <<- //
      DistMatrix<double> dcov[qe.nbp], cov(qe.npix,qe.npix), icov(qe.npix,qe.npix);



      QE_init(qe);

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





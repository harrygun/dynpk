  /*--------------------------------------------------------------------------
          >>>         Driver for Quadratic Estimator       <<<
  --------------------------------------------------------------------------*/

  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <iniparser.h>

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
  #include "io.h"
  #include "quadest.h"
  #include "misc.h"


#ifdef _MPI_
  #include <mpi.h>
#endif

#ifdef _OMP_
  #include <omp.h>
#endif



  int main( int argc, char *argv[])  {

      int debug= 20, i, j;
      char *ini_name;

      //if(argc!=2)  //
      //  myerr("Input parameters file is needed.", FALSE);

      ini_name=argv[1];
/*------------------------------------------------
            MPI initialization.
------------------------------------------------*/
      MPIpar mpi;
    #ifdef _MPI_
      //int mpi_ntask, mpi_rank, mpi_rc;
      mpi.rc = MPI_Init(&argc, &argv);

      if (mpi.rc != MPI_SUCCESS) {
           printf ("Error starting MPI program. Terminating.\n");
           MPI_Abort(MPI_COMM_WORLD, mpi.rc);
           }

      MPI_Comm_size(MPI_COMM_WORLD,&mpi.ntask);
      MPI_Comm_rank(MPI_COMM_WORLD,&mpi.rank);

      printf ("Number of tasks= %d My rank= %d\n", mpi.ntask, mpi.rank);

  /* MPI initialization ends. */
      if(mpi.rank==0) 
        printf("%d Sending parameter filename %s to other processes.\n", mpi.rank, ini_name);

      MPI_Bcast( ini_name, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

      if(mpi.rank!=0) 
        printf("%d Received parameter filename %s.\n", mpi.rank, ini_name);
    #else
      mpi.ntask = 1;
      mpi.rank = 0;
    #endif

  /*--------------------------------------
         End of MPI initialization.
  --------------------------------------*/

      printf("Opening File '%s'.\n", ini_name);
      dictionary * dict;
      dict = iniparser_load(ini_name);

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
      printf("==================================\n"); fflush(stdout);

      char *fn_dcov, *fn_cov, *fn_icov, *fn_out, *fn_map, *fn_plist;

      fn_dcov="result/r1d/dcov.dat";
      fn_map ="result/r1d/dmap.dat";
      fn_plist="result/r1d/plist.dat";
      fn_out ="result/r1d/Qi.dat";

      // ->>   <<- //
      //size_t mdim, npix, nbp;

      qe.mdim=50; 
      //npix=mdim*mdim;

      qe.npix=qe.mdim;
      qe.n_bp=24;   //mdim*mdim;
      qe.map_dim=1;


      qe.dcov=(double *)malloc(sizeof(double)*qe.n_bp*qe.npix);
      qe.cov= (double *)malloc(sizeof(double)*qe.npix*qe.npix);
      qe.icov=(double *)malloc(sizeof(double)*qe.npix*qe.npix);
      qe.plist=(double *)malloc(sizeof(double)*qe.n_bp);
      qe.map=(double *)malloc(sizeof(double)*qe.npix);


     
      import_data_double(&mpi, fn_map, qe.map, sizeof(double), qe.npix);
      import_data_double(&mpi, fn_dcov, qe.dcov, sizeof(double), qe.n_bp*qe.npix);
      import_data_double(&mpi, fn_plist, qe.plist, sizeof(double), qe.n_bp);


      // ->> noise covariance matrix <<- //
      qe.covn_v=(double *)malloc(sizeof(double)*qe.npix);
      cov_noise(&mpi, qe.covn_v, qe.npix, NULL);

      // ->> call Quadratic Estimator <<- //
      quad_est(&mpi, &qe);

   
      // ->> output <<- //
      printf("Output data.\n"); fflush(stdout);
      write_data(&mpi, fn_out, qe.Qi, sizeof(double), qe.n_bp);


  /*-----------------------------------------------------
                     free all
  -----------------------------------------------------*/
    stop:
      iniparser_freedict(dict);

      //free(qe.dcov);   free(qe.icov); 
      //free(qe.cov);    


      #ifdef _MPI_
      MPI_Finalize();
      #endif
      
    }





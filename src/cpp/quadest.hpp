
  #ifndef _H_QUADEST_
  #define _H_QUADEST_

    #include "El.hpp"
    #include "mpinit.hpp"

    class QEpar {

      public:
        int  debug;
        size_t mdim, map_dim, nbp, npix;

        double m_dim[2], kt_list_para[3], kf_list_para[3], dmap_res[2];
        double map_zoom_factor;

        char *ini_name;
        std::string get_bp_type, bp_list_fname;

	// ->> IO 
        std::string output_prefix;
      

        //double *dcov, *cov, *icov, *covn_v, *plist, *map;
        El::vector<double> Qip, Qi, covn_vec, plist, map;

        // ->> DistMatrix <<- //
	El::DistMatrix<double>  dcov_vec, Fij, iFij;
        //El::DistMatrix<double> cov, icov, dcov[nbp];


        //QEpar();  // constructor //
        //QEpar(char *ini_name, std::string sec="Quadratic_Estimator");
        QEpar(char *ini_name, std::string sec);

        //void QE_parameter(char *ini_name, std::string sec="Quadratic_Estimator");
        void QE_parameter(char *ini_name, std::string sec);

      };





  #endif

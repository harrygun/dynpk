
  #ifndef _H_QUADEST_
  #define _H_QUADEST_

    #include "El.hpp"
    #include "mpinit.hpp"

    class QEpar {

      public:
        int  debug;
        size_t ndim, nbp, npix, map_size[2];
	size_t map1d_f, map_dim;  // 1d map freq index
	bool  import_dcov;

        double kt_list_para[3], kf_list_para[3], dmap_res[2];
	double map_zoom_factor;

        //char *ini_name, *bp_list_fname, *data_fname;
        std::string bp_init_type, output_prefix, data_fname, 
                    bp_list_fname, dcov_fname;

	// ->> global MPI 
        MPIpar *glmpi;
      
        El::vector<double> covn_vec, pfid, klow, kup, klist;

        // ->> DistMatrix <<- //
	El::DistMatrix<double>  dcov_vec, Fij, iFij;
	El::DistMatrix<double> Qip, Qi, *dcov;

        // ->> The following will be defined internally <<- //
        //El::DistMatrix<double> cov, icov, dcov[nbp];

	// Matrix //
	El::Matrix<double> bpk, dmap;

    
        // ->> Methods <<- //
        // constructor //
        QEpar(char *ini_name, std::string sec, MPIpar *glmpi);

        // destructor //
        ~QEpar();

        // ->> initialization <<- //
        void QE_parameter(char *ini_name, std::string sec);

        void rawdata_init(std::string data_fname);

        void band_power_init(std::string bp_init_type, std::string bp_fname);

        void dcov_init(bool from_file, std::string fname);

        void noise_init();

        void fdcov_recovery();

	// ->> calculation <<- //
        //void Quad_Estimator(El::vector<double> pk_fid, int n_it);



      };









  #endif

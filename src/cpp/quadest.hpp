
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
      

        El::vector<double> Qip, Qi, covn_vec, map, pfid, klow, kup, klist;

        // ->> DistMatrix <<- //
	El::DistMatrix<double>  dcov_vec, Fij, iFij;
        //El::DistMatrix<double> cov, icov, dcov[nbp];
	//
	El::Matrix<double> bpk, dmap;


        // constructor //
        QEpar(char *ini_name, std::string sec, MPIpar *glmpi);

        void QE_parameter(char *ini_name, std::string sec);

        void rawdata_init(std::string data_fname);

        void band_power_init(std::string bp_init_type, std::string bp_fname);

        void dcov_init(bool from_file, std::string fname);


        // destructor //
        ~QEpar();


      };





  #endif


  #ifndef _H_GLB_VARB_
    #define _H_GLB_VARB_

    #include "El.hpp"



    class QEpar {

      public:
        int  debug;
        size_t mdim;
        size_t map_dim, nbp, npix;

        double m_dim[2], kt_list_para[3], kf_list_para[3], dmap_res[2];
        double map_zoom_factor;

        std::string get_bp_type, bp_list_fname;

        //double *dcov, *cov, *icov, *covn_v, *plist, *map;
        El::vector<double> Qip, Qi, covn_vec, plist, map;

        // ->> DistMatrix <<- //
	El::DistMatrix<double>  dcov_vec, Fij, iFij;
	//El::DistMatrix<double> cov, icov, dcov[nbp];

      };



    class MPIpar {

      public:
        int ntask, rank, rc, root;
        int max, start, totrun, ind_run;

        std::string fname, extfname[5];

        El::mpi::Comm world;

      };







  // ->>
  #define MemIdx3D(n, i1, i2, i3)  (i3+n*(i2+n*i1))
  #define MemIdx3D_n3(n1, n2, n3, i1, i2, i3)  (i3+n3*(i2+n2*i1))

  // ->> for equal-length cubic array <<- //
  #define ArrayAccess2D(a, n, i, j) (a)[ j+n*i ]
  #define ArrayAccess3D(a, n, i, j, k) ((a)[(i)*(n)*(n)+(j)*(n)+(k)])

  // ->> non-equal-length cubic array <<- //
  #define ArrayAccess3D_n3(a, n1, n2, n3, i, j, k) ((a)[ ((n2)*(i)+(j))*(n3)+(k) ])

  #define ArrayAccess5D_n5(a, n1, n2, n3, n4, n5, i1, i2, i3, i4, i5) ((a)[ i5+n5*(i4+n4*(i3+n3*(i2+n2*i1))) ])

  #define ArrayAccess4D_n4(a, n1, n2, n3, n4, i1, i2, i3, i4) ((a)[ i4+n4*(i3+n3*(i2+n2*i1)) ])
   
  #define ArrayAccess2D_n2(a, n1, n2, i1, i2) (a)[ i2+n2*i1 ]

  #define ArrayAccess2D_n2_list(a, n1, n2, nlist, i1, i2) (a)[ (i2+n2*i1)*nlist ]
   




  #endif



  /*
    typedef struct {
      int  init, debug;
      size_t mdim;
      size_t map_dim, nbp, npix;

      double m_dim[2], kt_list_para[3], kf_list_para[3], dmap_res[2];
      double map_zoom_factor;

      char *get_bp_type, *bp_list_fname;


      //double *dcov, *cov, *icov, *covn_v, *plist, *map;
      //double *Qip, *Qi, *Fij, *iFij;

      // ->> DistMatrix <<- //
      //El::DistMatrix<double> &dcov_vec, &dcov, &cov, &icov;
      //El::DistMatrix<double> dcov[nbp], cov, icov;

      }QEpar;


     typedef struct {
      int ntask, rank, rc, root;
      int max, start, totrun, ind_run;
      char *fname, *extfname[5];
      } MPIpar;
  */


  #ifndef _H_GLB_VARB_
    #define _H_GLB_VARB_

    #include "varb.h"
    #include "myinterpolate.h"

    typedef struct {
      int  init, debug;
      size_t ndim;

      double m_dim[2], kt_list_para[3], kf_list_para[3], dmap_res[2];
      double map_zoom_factor;

      char *get_bp_type, *bp_list_fname;


      double *dcov, *cov, *icov;

      
      //Interpar *initp;
      }QEpar;



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

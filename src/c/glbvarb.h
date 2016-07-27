
  #ifndef _H_GLB_VARB_
    #define _H_GLB_VARB_

    #include "varb.h"
    #include "myinterpolate.h"

    typedef struct {
      int  init, debug;
      double m_dim[2], kt_list_para[3], kf_list_para[3], dmap_res[2];
      double map_zoom_factor;

      char *get_bp_type, *bp_list_fname;


      double *dcov, *;

      
      //Interpar *initp;
      }QEpar;





  #endif


  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  #include <string>
  #include <boost/property_tree/ptree.hpp>
  #include <boost/property_tree/ini_parser.hpp>

  #include "El.hpp"

  #include "glbvarb.hpp"
  #include "io.hpp"
  #include "mpinit.hpp"
  #include "quadest.hpp"
  #include "quadest_init.hpp"
  #include "cov.hpp"

  using namespace El;
  using namespace std;




  void QEpar::QE_parameter(char *ini_name, string sec) {

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(ini_name, pt);

    output_prefix=pt.get<string>(sec+".output_prefix");
    nbp=pt.get<size_t>(sec+".num_band_power");
    map_dim=pt.get<size_t>(sec+".map_dimension");

    mdim=pt.get<size_t>(sec+".map_resolution_val");

    if(map_dim==1)
      npix=mdim;
    else if (map_dim==2)
      npix=mdim;

    return;
    }


  void QEpar::band_power_init(string bp_init_type, char *bp_fname){
    // ->> import band_power data <<- // 

    if(bp_init_type=="import") {
      // ->> import klist & plist from file <<- //

      bpk=Matrix<double>(4, nbp);  // import k, k_low, k_up, and band power P(k)
      Read(bpk, bp_fname);

      plist=vector<double>(nbp);
      klist=vector<double>(nbp);
      klow=vector<double>(nbp);
      kup=vector<double>(nbp);

      for(int i=0; i<nbp; i++) {

        klist[i] = bpk.Get(0, i);
        klow[i]  = bpk.Get(1, i);
        kup[i]   = bpk.Get(2, i);
        plist[i] = bpk.Get(3, i);
        }

      }
    else 
        throw runtime_error("Error: band_power_init()");


    return;
    }




  void QEpar::dcov_init() {
    //     ->> get the derivative of covariance matrix <<-   //
    // ->> I'd like to have a copy of dcov_vec for every process //

    int iloc, jloc, iglo, jglo; 
    double dc;

    DistMatrix<double> *dcov_vdist = new DistMatrix<double>(nbp, npix); 

    const int localHeight = dcov_vdist->LocalHeight();
    const int localWidth = dcov_vdist->LocalWidth();


    for(jloc=0; jloc<localWidth; ++jloc) {
      for(iloc=0; iloc<localHeight; ++iloc) {

        jglo=dcov_vdist->GlobalCol(jloc);
        iglo=dcov_vdist->GlobalRow(iloc);

        // ->> band power k-list <<- //
        //kt_low = ;
        //kt_up = ;

        dc=get_dcov_klim_r1d(klow[iglo], kup[iglo], dt, dtab);

        dcov_vdist->SetLocal(iloc, jloc, dc);
        }
      }

    // redistribution //
    dcov_vec=DistMatrix<double, STAR, STAR>(nbp, npix); 
    dcov_vec=*(dcov_vdist);
    delete dcov_vdist;


    Display(dcov_vec, "dcov_vec");

    return;
    }



  QEpar::QEpar(char *ini_name, string sec) {
    //->>  QEpar constructor <<- //

    QE_parameter(ini_name, sec);

    // ->> data import <<- //

    //  band_power initialization  //
    

    // dcov initialzation //
    dcov_init();


    }


  QEpar::~QEpar() {
    }




  /*
  void QE_init(QEpar &qe, const string &mat_type="DistMatrix") {
     
    if (mat_type.compare("DistMatrix")!=0)
      throw runtime_error("Error: Only DistMatrix type is Supported.");

    //DistMatrix<double> qe.dcov_vec(qe.nbp, qe.npix), qe.dcov_i(qe.npix, qe.npix), 
    //                   qe.cov(qe.npix, qe.npix),   qe.icov(qe.npix, qe.npix);

    //Read(qe.dcov_vec, fname); 

    return;
    }

  */


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

  using namespace El;
  using namespace std;




  //void QEpar::QE_parameter(char *ini_name, string sec="Quadratic_Estimator") {
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







  //QEpar::QEpar(char *ini_name, string sec="Quadratic_Estimator") {
  QEpar::QEpar(char *ini_name, string sec) {
    //->>  QEpar constructor <<- //

    QEpar::QE_parameter(ini_name, sec);

    //  band_power initialization  //


    // dcov initialzation //



    }









  void QE_init(QEpar &qe, const string &mat_type="DistMatrix") {
     
    if (mat_type.compare("DistMatrix")!=0)
      throw runtime_error("Error: Only DistMatrix type is Supported.");

    //DistMatrix<double> qe.dcov_vec(qe.nbp, qe.npix), qe.dcov_i(qe.npix, qe.npix), 
    //                   qe.cov(qe.npix, qe.npix),   qe.icov(qe.npix, qe.npix);

    //Read(qe.dcov_vec, fname); 

    return;
    }



  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  #include <boost/property_tree/ptree.hpp>
  #include <boost/property_tree/ini_parser.hpp>

  #include "glbvarb.hpp"
  #include "io.hpp"
  #include "mpinit.hpp"
  #include "quadest.hpp"

  using namespace std;




  void para_importer(char *ini_name){

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(ini_name, pt);
    string sec="Quadratic_Estimator";






    return;
    }

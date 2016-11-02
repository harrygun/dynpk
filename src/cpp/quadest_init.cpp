
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

    // debug //
    /*
    try{
      debug=pt.get<int>(sec+".debug"); 
      throw 1;}
    catch(int e){
      debug=e; }
    */
    debug=pt.get<int>(sec+".debug"); 
    cout << "Param init, debug=" << debug << endl;

    // map parameters //
    //cout << "Map init" << endl;
    data_fname=pt.get<string>(sec+".input_data_fname");
    ndim=pt.get<size_t>(sec+".number_of_dimension_map");
    map_zoom_factor=pt.get<double>(sec+".map_zoom_factor");

    //TODO >> delete map_dim <<- //
    map_dim=ndim;

    map_size[0]=pt.get<size_t>(sec+".map_size_0");
    dmap_res[0]=pt.get<double>(sec+".map_resolution_0");

    if(ndim==1){
      npix=map_size[0];
      map_size[1]=1; dmap_res[1]=1; 
      }
    else if (ndim==2){
      map_size[1]=pt.get<size_t>(sec+".map_size_1");
      dmap_res[1]=pt.get<double>(sec+".map_resolution_1");
      npix=map_size[0]*map_size[1];
      }

    //cout << "Bandpower init" << endl;

    // band power parameters //
    nbp=pt.get<size_t>(sec+".num_band_power");
    bp_list_fname=pt.get<string>(sec+".band_power_list_fname");
    bp_init_type =pt.get<string>(sec+".band_power_init_type");


    // band power klist //
    kt_list_para[0]=pt.get<double>(sec+".kt_list_min");
    kt_list_para[1]=pt.get<double>(sec+".kt_list_max");
    kt_list_para[2]=pt.get<double>(sec+".kt_list_num");

    if (ndim==2){
      kf_list_para[0]=pt.get<double>(sec+".kf_list_min");
      kf_list_para[1]=pt.get<double>(sec+".kf_list_max");
      kf_list_para[2]=pt.get<double>(sec+".kf_list_num");
      }
    else{
      kf_list_para[0]=0; kf_list_para[1]=0; kf_list_para[2]=0;
      }

    // dcov //
    import_dcov=pt.get<bool>(sec+".import_dcov_from_file");
    dcov_fname=pt.get<string>(sec+".dcov_fname");


    //cout << "others" << endl;

    // output //
    output_prefix=pt.get<string>(sec+".output_prefix");
    

    if(debug>=50){
      cout << "QuaDest Parameters: (debug=" << debug << ")" << endl;
      cout << "  map size: " << map_size[0] << ", " << map_size[1] << endl;
      cout << "  n of dim: " << ndim << endl;
      cout << "  map_resolution: " << dmap_res[0] << ", " << dmap_res[1] << endl;
      cout << "data_fname: " << data_fname << endl;
      }

    return;
    }


  void QEpar::rawdata_init(string data_fname){

    if(ndim>=3) 
      throw runtime_error("Error: Doesn't support higher-dim map yet.");

    //Matrix<double> *dmap_=new Matrix<double>(map_size[0], map_size[1]);

    dmap=DistMatrix<double, STAR, STAR>(npix, 1);
    Read(dmap, data_fname);

    cout << "Imported raw map." << endl;
    
    //Display(dmap);
    return;
    }


  void QEpar::band_power_init(string bp_init_type, string bp_fname){
    // ->> import band_power data <<- // 

    if(bp_init_type=="import") {
      // ->> import klist & pfid from file <<- //

      bpk=Matrix<double>(4, nbp);  // import k, k_low, k_up, and band power P(k)
      Read(bpk, bp_fname);

      cout << "Imported bandpower." << endl;

      pfid=vector<double>(nbp);
      klist=vector<double>(nbp);
      klow=vector<double>(nbp);
      kup=vector<double>(nbp);

      for(int i=0; i<nbp; i++) {

        klist[i]= bpk.Get(0, i);
        klow[i] = bpk.Get(1, i);
        kup[i]  = bpk.Get(2, i);
        pfid[i] = bpk.Get(3, i);
        }

      }
    else
        throw runtime_error("Error: band_power_init()");


    return;
    }




  void QEpar::dcov_init(bool from_file, string fname) {
    //     ->> get the derivative of covariance matrix <<-   //
    // ->> I'd like to have a copy of dcov_vec for every process //

    int iloc, jloc, iglo, jglo; 
    double dc, dt, df, dtab;

    dcov_vec=DistMatrix<double, STAR, STAR>(nbp, npix); 

    // import from file //
    if (from_file) {
      cout << "importing dcov from file '" << fname << "'"<< endl;

      Read(dcov_vec, fname);

      //Display(dcov_vec, "dcov_vec");
      return;
      }


    // otherwise, calculate internally //

    DistMatrix<double> *dcov_vdist = new DistMatrix<double>(nbp, npix); 

    const int localHeight = dcov_vdist->LocalHeight();
    const int localWidth = dcov_vdist->LocalWidth();

    for(jloc=0; jloc<localWidth; jloc++) {
      for(iloc=0; iloc<localHeight; iloc++) {

        jglo=dcov_vdist->GlobalCol(jloc);
        iglo=dcov_vdist->GlobalRow(iloc);

        if(ndim==1){
          // ->> band power k-list <<- //
          dt = dmap_res[0];
          dtab = jglo*dt;

          dc=get_dcov_klim_r1d(klow[iglo], kup[iglo], dt, dtab);
	  }

	else if(ndim==2){
          dc=0;//get_dcov_klim_r1d(klow[iglo], kup[iglo], dt, dtab);
	  }

        dcov_vdist->SetLocal(iloc, jloc, dc);
        }
      }

    // redistribution //
    dcov_vec=*(dcov_vdist);
    delete dcov_vdist;


    Display(dcov_vec, "dcov_vec");

    return;
    }


 



  void QEpar::fdcov_recovery(){

    int a, iglo, jglo, iloc, jloc, localHeight, localWidth; 
    double dcc;

    // allocate an array of DistMatrix //
    dcov = new DistMatrix<double>[nbp]; 


    for(a=0; a<nbp; a++) {
      dcov[a]=DistMatrix<double>(npix, npix);

      localHeight=dcov[a].LocalHeight();
      localWidth =dcov[a].LocalWidth();

      // full recovery //
      for(jloc=0; jloc<localWidth; jloc++){
        // ->> get global index <<- //
        jglo=dcov[a].GlobalCol(jloc);

        for(iloc=0; iloc<localHeight; iloc++) {
          // ->> get global index <<- //
          iglo=dcov[a].GlobalRow(iloc);

          if (ndim==1){
            dcc=dcov_vec.Get(a, (int)fabs(iglo-jglo)); }
          else {
            throw runtime_error("Error: fdcov n>1 NOT supported yet."); }

          dcov[a].SetLocal(iloc, jloc, dcc);
	  }
        }
     
      }

    if(debug>=50) {
      string fn_dcov;
      //ostringstream convert;

      for(a=0; a<nbp; a++) {
        fn_dcov="result/r1d/dcov/dcovf_" + to_string(a);
        Write(*dcov, fn_dcov);
	}
      }




    return;
    }




  void QEpar::noise_init(){

    return;
    }


  QEpar::QEpar(char *ini_name, string sec, MPIpar *mpi_add) {

    //->>  QEpar constructor <<- //
    glmpi=mpi_add;

    //if (glmpi->rank0) // ->> add this later
    QE_parameter(ini_name, sec);

    // ->> data import <<- //
    rawdata_init(data_fname);
    
    //  band_power initialization  //
    band_power_init(bp_init_type, bp_list_fname);
    
    // dcov initialzation //
    dcov_init(import_dcov, dcov_fname);

    // ->> recovery of the full dcov <<- //
    fdcov_recovery();

    cout << "QEpar initialization is complete." << endl;
    return;
    }



    // noise initialization //
   void noise_init(){
      }


  QEpar::~QEpar() {
    delete[] dcov;
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

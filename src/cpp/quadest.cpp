  /*>> matrix opertion <<*/

  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  #include "El.hpp"

  #include "glbvarb.hpp"
  #include "io.hpp"
  #include "mpinit.hpp"
  #include "quadest.hpp"
  #include "quadest_init.hpp"

  using namespace std;
  using namespace El;





  void fcov_recovery(vector<double> &pk, DistMatrix<double> *covf, 
                     DistMatrix<double> *dcov, size_t nbp, int debug) {
    // ->> obtain full covariance matrix from dvoc and pk_list <<- //

    int a, iglo, jglo, iloc, jloc; //localHeight, localWidth; 
    double val;

    const int localHeight=covf->LocalHeight();
    const int localWidth =covf->LocalWidth();


    for(jloc=0; jloc<localWidth; jloc++){
      for(iloc=0; iloc<localHeight; iloc++) {

        val=0;
        for(a=0; a<nbp; a++) {
          val+=dcov[a].GetLocal(iloc,jloc)*pk[a]; 
          }

        covf->SetLocal(iloc, jloc, val);
	}
      }

    return;
    }



  void iFisher(DistMatrix<double> *dcov, DistMatrix<double> *icov, 
       DistMatrix<double> *iFij, size_t npix, size_t nbp, size_t ndim, int debug)  {

    int a, b, iglo, jglo, iloc, jloc, localHeight, localWidth; 
    double Fval;

    DistMatrix<double> *ic_dcov = new DistMatrix<double>[nbp];
    DistMatrix<double> *Fij = new DistMatrix<double>(nbp, nbp);

    for(a=0; a<nbp; a++) {
      ic_dcov[a]=DistMatrix<double>(npix, npix);
      Gemm(NORMAL, NORMAL, double(1.), *icov, dcov[a], double(0.), ic_dcov[a]);
      }


    for(a=0; a<nbp; a++) {
      for(b=0; b<nbp; b++) {
        Fval=HilbertSchmidt(ic_dcov[a], ic_dcov[b])*0.5;

	// ?? local set ?? //
        Fij->Set(a, b, Fval);
        }
      }

    // inverse //
    *iFij=*Fij;
    SymmetricInverse(LOWER, *iFij);

    if(debug>=50) {
      string fn_Fij, fn_iFij;
      fn_Fij="result/r1d/Fij.out";
      fn_iFij="result/r1d/iFij.out";

      Write(*Fij, fn_Fij);
      Write(*iFij, fn_iFij);
      }



    delete Fij;
    delete[] ic_dcov;

    cout << "Fisher Matrix is done." << endl; 
    return;
    }







  void QEpar::Quad_Estimator(vector<double> pk_fid, int n_it) {
    // ->> method for calculating the quadratic estimator <<- //

    int a, b; 
    vector<double> pk;
    //DistMatrix<double> covf(npix, npix), covf_inv(npix, npix);
    DistMatrix<double> *covf, *covf_inv;


    if (n_it>1)
      throw runtime_error("Error: Iterations NOT supported yet.");
    else if (n_it==1)
      pk=pk_fid;

    Display(pk);


    covf     = new DistMatrix<double>(npix, npix);
    covf_inv = new DistMatrix<double>(npix, npix);

    // ->> first recover the full covariance matrix <<- //
    fcov_recovery(pk, covf, dcov, nbp, debug);
    cout << "Full covariance matrix done." << endl;

    if(debug>=50) {
      string fn_out="result/r1d/covf.dat";
      Write(*covf, fn_out);
      }

    // inversion //
    *covf_inv=*covf;
    SymmetricInverse(LOWER, *covf_inv);

    cout << "inverse of covariance matrix done." << endl; 


    // ->> calculate and its inverse <<- //
    DistMatrix<double> *iFij, *Uni, *Ml, *Qpar, *d_ic;
    double Qval;
    
    iFij = new DistMatrix<double>(nbp, nbp);
    Uni = new DistMatrix<double>(nbp, 1);
    Ml= new DistMatrix<double>(nbp, 1);

    iFisher(dcov, covf_inv, iFij, npix, nbp, ndim, debug);

    Ones(*Uni, nbp, 1);
    Gemv(NORMAL, double(1.), *iFij, *Uni, double(0.), *Ml);

    delete Uni;


    // ->> some pre-calculation of Qi <<- //
    // ->> (C^{-1}.d) <<- //
    cout << "now calculate quadratic estimator" << endl;
    d_ic=new DistMatrix<double>(npix, 1);
    Gemv(NORMAL, double(1.), *covf_inv, dmap, double(0.), *d_ic);

    Qpar= new DistMatrix<double>(npix, 1);
    Qi=DistMatrix<double>(nbp, 1);

    for(a=0; a<nbp; a++) {
      //Zeros(*Qpar);
      Gemv(NORMAL, double(1.), dcov[a], *d_ic, double(0.), *Qpar);
      Qi.Set(a, 0, Dot(*d_ic, *Qpar)*Ml->Get(a, 0)*0.5);
      }

    delete iFij, d_ic, Qpar, Ml;
    return;
    }





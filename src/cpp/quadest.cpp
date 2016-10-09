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




  double access_dcov(double *dcov, size_t nbp, size_t npix, int i, 
                                      int a, int b, size_t ndim) {
    // ->> i:    index of bandpower
    // ->> a/b:  indices of map pixels 

    if(ndim==1)
      {return ArrayAccess2D_n2(dcov, nbp, npix, i, (int)fabs((double)(a-b)));}

    else if(ndim==2) {abort();}

    else {abort();}

    return -99.;
    }






  void cov_noise(MPIpar *mpi, double *covn_v, size_t npix, char *type){
    int i;

    if(type==NULL) {
      for(i=0; i<npix; i++)
        covn_v[i]=0;
      }
    else abort();


    return;
    }





  //----------------------------------------------//





  void fcov_recovery(vector<double> &pk, DistMatrix<double> &covf) {
    // ->> obtain full covariance matrix from dvoc and pk_list <<- //

    int a, iglo, jglo, iloc, jloc, localHeight, localWidth; 
    double val;

    const int localHeight=covf.LocalHeight();
    const int localWidth =covf.LocalWidth();


    for(jloc=0; jloc<localWidth; jloc++){
      for(iloc=0; iloc<localHeight; iloc++) {

        val=0;
        for(a=0; a<nbp; a++) {
          val+=dcov[a].GetLocal(iloc,jloc)*pk[a]; 
          }

        covf.SetLocal(iloc, jloc, val);
	}
      }

    return;
    }



  void iFisher(DistMatrix<double> *dcov, DistMatrix<double> &icov, 
       DistMatrix<double> *iFij, size_t npix, size_t nbp, size_t ndim)  {

    int a, b, iglo, jglo, iloc, jloc, localHeight, localWidth; 
    double Fval;
    DistMatrix<double> *ic_dcov;

    ic_dcov=new DistMatrix<double>[nbp];
    DistMatrix<double> *Fij=new DistMatrix<double>(nbp, nbp);

    for(a=0; a<nbp; a++) {
      ic_dcov[a]=DistMatrix<double>(npix, npix);
      Gemm(NORMAL, NORMAL, double(1.), icov, dcov[a], double(0.), ic_dcov[a]);
      }


    for(a=0; a<nbp; a++) {
      for(b=0; b<nbp; b++) {
        Fval=HilbertSchmidt(ic_dcov[a], ic_dcov[b])*0.5;

	// ?? local set ?? //
        Fij.Set(a, b, Fval);
        }
      }

    // inverse //
    iFij=Fij;
    HPDInverse(LOWER, iFij);

    delete[] Fij;
    delete[] ic_dcov;

    return;
    }







  void QEpar::Quad_Estimator(vector<double> pk_fid, int n_it) {
    // ->> method for calculating the quadratic estimator <<- //

    //int a, b, i, j, idx, id;
    string fn;
    vector<double> pk;
    DistMatrix<double> covf(npix, npix), covf_inv(npix, npix);

    if (n_it>1)
      throw runtime_error("Error: Iterations NOT supported yet.");
    else if (n_it==1)
      pk=pk_fid;

    // ->> first recover the full covariance matrix <<- //
    fcov_recovery(pk, covf);
    cout << "Full covariance matrix done." << endl;

    // inversion //
    // ??
    covf_inv=covf;
    HPDInverse(LOWER, covf_inv);
    //MakeHermitian(LOWER, covf_inv);

    cout << "inverse of covariance matrix done." << endl; 


    // ->> calculate and its inverse <<- //
    iFij = new DistMatrix<double>(nbp, nbp);
    iFisher(dcov, covf_inv, iFij, npix, nbp, ndim);


    // ->> some pre-calculation of Qi <<- //
    // ->> (C^{-1}.d) <<- //
    Matrix<double>* d_ic=new Matrix<double>();
    Gemv(NORMAL, double(1.), covf_inv, dmap, double(0.), d_ic);


    // 
    Qi= ;


    // ->> calculate Qi' <<- //
    double *Qip_s, *Qi_s;
    Qip_s=(double *)malloc(mpi->ind_run*sizeof(double));
    Qi_s=(double *)malloc(mpi->ind_run*sizeof(double));

    mpi->max=qe->nbp;    mpi->start=0;
    mpi_loop_init(mpi, "Qi");

    for(idx=0; idx<mpi->ind_run; idx++) {
      id=mpi_id(mpi, idx);

      Qip_s[idx]=0.;
      for(a=0; a<qe->npix; a++)
        for(b=0; b<qe->npix; b++) {
          Qip_s[idx]+=d_ic[a]*access_dcov(qe->dcov, qe->nbp, qe->npix, id, a, b, qe->ndim)*d_ic[b]/2.;
          }

      Qi_s[idx]=0;
      for(j=0; j<qe->nbp; j++) {
        Qi_s[idx]+=Qip_s[idx]*ArrayAccess2D(qe->iFij, qe->nbp, id, j);
        }
      }


    // ->> gather all data by root <<- //
    
    qe->Qip=(double *)malloc(qe->nbp*sizeof(double));
    qe->Qi =(double *)malloc(qe->nbp*sizeof(double));


    qe->Qip=mpi_gather_dist_double(mpi, Qip_s, mpi->ind_run, mpi->max);
    qe->Qi=mpi_gather_dist_double(mpi, Qi_s, mpi->ind_run, mpi->max);

    free(Qip_s); free(Qi_s); free(d_ic);
    return;
    }





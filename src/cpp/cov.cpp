  #include <iostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include <cstdio>
  #include <string>

  #include <gsl/gsl_sf_expint.h>

  #include "El.hpp"

  #include "glbvarb.hpp"
  #include "mpinit.hpp"
  #include "io.hpp"


  using namespace std;
  using namespace El;




  inline double SinIntegral(const double x){
    return gsl_sf_Si(x);
    }

  inline double CosIntegral(const double x){
    return gsl_sf_Ci(x);
  }


  #define Cos cos
  #define Sin sin
  #define _epsilon_ 1e-10






  double dcov1d_klim_real(double ktia, double ktib, double dt, double dtab) {
    /* ->>  kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    */
    double dc1d_real;

    dc1d_real= (-2*ktib*(-1 + Cos(dt*ktia))*Cos(dtab*ktia) + 
             2*ktia*(-1 + Cos(dt*ktib))*Cos(dtab*ktib) + 
             (-dt + dtab)*ktia*ktib*SinIntegral((dt - dtab)*ktia) + 
             2*dtab*ktia*ktib*SinIntegral(dtab*ktia) - 
             ktia*ktib*((dt + dtab)*SinIntegral((dt + dtab)*ktia) + 
             (-dt + dtab)*SinIntegral((dt - dtab)*ktib) + 2*dtab*SinIntegral(dtab*ktib)
             ) + (dt + dtab)*ktia*ktib*SinIntegral((dt + dtab)*ktib))/(ktia*ktib*M_PI);
              
    //if (isnan(dc1d_real)) {
    //  cout << "real:", ktia, ktib, dt, dtab << endl; 
    //  fflush(stdout);
    //  }

    return dc1d_real/2.;
    }



  double dcov1d_klim_imag(double ktia, double ktib, double dt, double dtab){
    /* ->>  kti:   center value of kt_i <<- 
            Dkti:  Delta kt_i <<- 
	    dt:    Delta t, pixel size
	    dtab:  Delta t_ab
    */
    double dc1d_imag;

    if (fabs(dtab-dt)<_epsilon_){
      // ->> if tab->dt
      dc1d_imag= (-2*dt*ktia*ktib*CosIntegral(dt*ktia) 
                + 2*dt*ktia*ktib*CosIntegral(2*dt*ktia) + 
                2*dt*ktia*ktib*CosIntegral(dt*ktib) - 
                2*dt*ktia*ktib*CosIntegral(2*dt*ktib) + 
                2*ktib*Sin(dt*ktia) - ktib*Sin(2*dt*ktia) - 2*ktia*Sin(dt*ktib) + 
                ktia*Sin(2*dt*ktib))/(ktia*ktib*M_PI);
      }

    else if (fabs(dtab+dt)<_epsilon_){
      //->> dtab->-dt
      dc1d_imag=(2*dt*ktia*ktib*CosIntegral(-(dt*ktia)) - 
                2*dt*ktia*ktib*CosIntegral(2*dt*ktia) - 
                2*dt*ktia*ktib*CosIntegral(-(dt*ktib)) + 
                2*dt*ktia*ktib*CosIntegral(2*dt*ktib) - 2*ktib*Sin(dt*ktia) + 
                ktib*Sin(2*dt*ktia) + 2*ktia*Sin(dt*ktib) - 
                ktia*Sin(2*dt*ktib))/(ktia*ktib*M_PI);
      }

    else if (fabs(dtab)<_epsilon_){
      // ->> dtab->0
      dc1d_imag=0.;
      }

    else{
      dc1d_imag=((-dt + dtab)*ktia*ktib*CosIntegral((dt - dtab)*ktia) - 
                2*dtab*ktia*ktib*CosIntegral(dtab*ktia) + 
                dt*ktia*ktib*CosIntegral((dt + dtab)*ktia) + 
                dtab*ktia*ktib*CosIntegral((dt + dtab)*ktia) + 
                (dt - dtab)*ktia*ktib*CosIntegral((dt - dtab)*ktib) + 
                2*dtab*ktia*ktib*CosIntegral(dtab*ktib) - 
                dt*ktia*ktib*CosIntegral((dt + dtab)*ktib) - 
                dtab*ktia*ktib*CosIntegral((dt + dtab)*ktib) + 2*ktib*Sin(dtab*ktia)-
                2*ktib*Cos(dt*ktia)*Sin(dtab*ktia) - 2*ktia*Sin(dtab*ktib) + 
                2*ktia*Cos(dt*ktib)*Sin(dtab*ktib))/(ktia*ktib*M_PI);
      }

    //if (isnan(dc1d_imag)){
    //  cout << "imag:", ktia, ktib, dt, dtab << endl; fflush(stdout);
    //  }

    return dc1d_imag/2.;
    }



  double get_dcov_klim_r1d(double kt_low, double kt_up, double dt, double dtab)  {
    return 2.*dcov1d_klim_real(kt_low, kt_up, dt, dtab);
    }





  /* <<<< TO BE modified << -
  void get_dcov_klim_2D(vector<double> dcov, klist_low, klist_up, dt_df, npt, map_size, bool do_mpi=true) {
    // ->> get the derivative of covariance matrix <<- //
    // ->> to parallize 

    int i, a, b;
    double dt, df;

    dt, df = dt_df

    print 'get_dcov:', npt, map_size, klist_up.shape, klist_low.shape, dt_df.shape

    if do_mpi==True:
        prange=mpi.mpirange(npt)
    else:
        prange=range(npt)


    # ->> python loop <<- #
    for i in prange:
        #kti,  kfi  = klist[:,i]
        ktia, kfia = klist_low[:,i]
        ktib, kfib = klist_up[:,i]

        print i, ktia, ktib, kfia, kfib, dt, df
    
        for a in range(-map_size[0], map_size[0]):
            for b in range(-map_size[1], map_size[1]):
                dtab = a*dt
                dfab = b*df

                dcov[i][a][b]=2*(dcov1d_klim_real(ktia, ktib, dt, dtab)*  
                               dcov1d_klim_real(kfia, kfib, df, dfab) - 
                               dcov1d_klim_imag(ktia, ktib, dt, dtab)* 
                               dcov1d_klim_imag(kfia, kfib, df, dfab) );
			    
    return dcov;
    }

  */



  /*
  void dcov_recovery(QEpar &qe, DistMatrix<double> dcov[], DistMatrix<double> dcov_vec) {
    //->>
    vector <double> colvectorsymm;

    for(int i=0; i<qe.nbp; i++) {
      dcov[i]=DistMatrix<double>(qe.npix,qe.npix);
      Toeplitz(dcov[i], qe.npix, qe.npix, colvectorsymm);
      }



    // ->> initialization <<- //
    //DistMatrix<double> dcov[qe.nbp], cov(qe.npix,qe.npix), icov(qe.npix,qe.npix);
    //DistMatrix<double> dcov[qe.nbp], cov, icov;

    return;
    }


  */


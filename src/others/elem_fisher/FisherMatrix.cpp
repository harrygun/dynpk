#include <iostream>
#include <El.hpp>
#include "./include/FisherHelperFunctions.hpp"
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
  Environment env( argc, argv );

  const int nrows = Input("--nrows","nrows",50);        // #rows in covariance matrix
  const int ncols = Input("--ncols","ncols",50);        // #cols in covariance matrix
  const int ndcov_i = Input("--ndcov_i","ndcov_i",24);  // # derivative matrices
  const std::string Cfile =                             // location of covariance matrix
    Input("--Cfile","Cfile",std::string(""));
  const std::string dcovi_dir =                         // location of derivative matrix
    Input("--dcovi_dir","dcovi_dir",std::string(""));
  const std::string FisherOut =
    Input("--FisherOut","FisherOut",std::string(""));   // filepath where Fisher Matrix will be written
  ProcessInput();
  
  // Read in the covariance matrix
  DistMatrix<double> Cov(nrows,ncols);   // initialize covariance matrix "Cov"
  //Read(Cov,Cfile, bool sequential=true);     // fill covariance matrix with contents of "Cfile"
  Read(Cov,Cfile);     // fill covariance matrix with contents of "Cfile"



  // ->> testing <<- //
  int flag=false;

  if(flag==true) {
    DistMatrix<double> dcov_arr[ndcov_i];
    for(int i=0; i<ndcov_i; i++) {
      dcov_arr[i]=DistMatrix<double>(nrows,ncols);
      Zeros(dcov_arr[i], nrows, ncols);
      }

    Write(dcov_arr[0], "./ZERO_test_0", MATRIX_MARKET);
    Write(dcov_arr[10], "./ZERO_test_10", MATRIX_MARKET);
    Write(dcov_arr[23], "./ZERO_test_23", MATRIX_MARKET);
  }


  
  // Get the Fisher Matrix
  DistMatrix<double> FisherMatrix(ndcov_i,ndcov_i),dcovfull(nrows,ncols*ndcov_i);  //initialize matrices
  //DistMatrix<double> FisherMatrix(ndcov_i,ndcov_i),dcovfull(nrows*ncols, ndcov_i); 

  Combine_dcovi(dcovi_dir, ndcov_i, nrows, ncols, dcovfull);   // makes matrices out of derivative matrix (turns rows into Toeplitz matrices)
  GetFisherMatrix(dcovfull, ndcov_i, Cov, FisherMatrix);       // performs inversion of Cov and matrix multiplication w/ derivatives to get Fisher   

  std::cout << "Fisher is Doen." << endl;
  Write(FisherMatrix,FisherOut,MATRIX_MARKET);           // write FisherMatrix to file as MATRIX_MARKET type (read in python w/ scipy.io's mmread)

  return 0;
}

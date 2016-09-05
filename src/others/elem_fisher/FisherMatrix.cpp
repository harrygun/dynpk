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
  Read(Cov,Cfile, sequential=true);     // fill covariance matrix with contents of "Cfile"
  
  // Get the Fisher Matrix
  //DistMatrix<double> FisherMatrix(ndcov_i,ndcov_i),dcovfull(nrows,ncols*ndcov_i);   // initialize Fisher matrix, derivative matrix


  DistMatrix<double> FisherMatrix(ndcov_i,ndcov_i),dcovfull(ndcov_i,nrows*ncols);   // initialize Fisher matrix, derivative matrix

  Combine_dcovi(dcovi_dir, ndcov_i, nrows, ncols, dcovfull);   // makes matrices out of derivative matrix (turns rows into Toeplitz matrices)
  GetFisherMatrix(dcovfull, ndcov_i, Cov, FisherMatrix);       // performs inversion of Cov and matrix multiplication w/ derivatives to get Fisher   
  Write(FisherMatrix,FisherOut,MATRIX_MARKET);           // write FisherMatrix to file as MATRIX_MARKET type (read in python w/ scipy.io's mmread)

  return 0;
}

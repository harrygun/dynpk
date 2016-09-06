#include <iostream>
#include <fstream>
#include <El.hpp>
#include "../include/MatrixSetUpFunctions.hpp"
#include "../include/VectorFunctions.hpp"

using namespace El;
typedef double Real;

void mpixel_idx(int pixel,int timedimension, int freqdimension, vector<int>& idx)
{
  // fill 2-d vector idx with distance in time, frequency indices
  idx[0] = int(pixel/double(timedimension));
  idx[1] = pixel-freqdimension*idx[0];
}

double GetCovarianceMatrix(int pixeldimension, int submatrixlength, int dx0, int dx1, DistMatrix<double>& bandpower, DistMatrix<double>& dcov)
{
  // Get the covariance matrix element "C_element_out"
  double C_element_out=0;
  int shift=0;
 
  for (int i=0;i<pixeldimension;i++)
    {
      C_element_out += dcov.Get(dx0,dx1+i*submatrixlength)*bandpower.Get(i,0)/2.0;  //this works but is slow
      // below is a the same line rewritten with the faster function IR(start_index,end_index) vs. "Get" -- it is untested
      //C_element_out += dcov(IR(dx0,dx0+submatrixlength),IR(dx1+i*submatrixlength,dx1+(i+1)*submatrixlength)*bandpower.Get(i,0)/2.0;
    }
  return C_element_out;
}

void GetSymmetricCirculantVector(vector<double>& v, vector<double>&v_out)
{
  // Accepts correlation vector of length n                                                                                                                    
  // Returns full correlation vector (length 2n-1)                                                                                                             
  // length 2n-1 in order to generate Toeplitz matrix                                                                                                          
  // leads to symmetric matrix                                                                                                                                 
  int n = v.size();
  for (int i=0;i<2*n-1;i++)
    {
      v_out[i]= v[n-i-1];
      if (i==(n-1))
	{
	  v_out[i]=v[0];
	}
      else if (i>=n)
	{
	  v_out[i]=v[i+1-n];
	}
      Output("v_out(",i,")=",v_out[i]);
    }
}



void Combine_dcovi(std::string basename, int nmatrices, int nrows, int ncols, DistMatrix<double>& full_dcovi)
{
  //DistMatrix<double> dCk(nrows,nrows);
  DistMatrix<double> dCk(2*nrows-1,2*nrows-1);

  DistMatrix<double> dcovcols(nmatrices, nrows), colmatrix(nrows,1);
  vector<double> colvector(nrows),colvectorsymm(2*nrows-1);

  Zeros(colmatrix,nrows,1);

  //Read(dcovcols,basename, bool sequential=true);
  Read(dcovcols,basename);

  // ->> column major <<- //
  Transpose(dcovcols,dcovcols);

  ofstream fout;

  for (int k=0;k<nmatrices;k++)
    {
      colmatrix=dcovcols(IR(0,nrows),IR(k,k+1));
      DistMatrixToVector(colmatrix,colvector);
      GetSymmetricCirculantVector(colvector, colvectorsymm);
      Toeplitz(dCk,nrows,nrows,colvectorsymm);
      full_dcovi(IR(0,nrows),IR(k*nmatrices,k*nmatrices+ncols))=dCk(IR(0,nrows),IR(0,ncols));

      if(k==0){
        //Write(colvectorsymm,"./dC_ext.dat"); 

	fout.open("dC_ext.dat", std::ios::in | std::ios::binary);
	std::count<< colvectorsymm.size() << endl;

        fout.write(colvectorsymm[0], colvectorsymm.size());
        fout.close();

	
	std::ofstream output_file("./example.txt");
        Write(dCk(IR(0,nrows),IR(0,ncols)),"./dC0.dat", MATRIX_MARKET); 
        }

    }

}


void GetFisherMatrix(DistMatrix<double>& dcovfull, int nmatrices, DistMatrix<double>& Cov, DistMatrix<double>& FisherMatrix)
{
  double FisherElement=0;
  int nrows=Cov.Height();
  int ncols=Cov.Width();
  DistMatrix<double> dCi(nrows,ncols),dCj(nrows,ncols);
  DistMatrix<double> Cinv(nrows,ncols),CinvdCi(nrows,ncols),CinvdCj(nrows,ncols),Product(nrows,ncols);

  Cinv=Cov;
  HPDInverse(LOWER,Cinv);
  MakeHermitian(LOWER,Cinv);
  Zeros(dCi,nrows,ncols);
  Zeros(dCj,nrows,ncols);
  Zeros(FisherMatrix,nmatrices,nmatrices);
  Zeros(CinvdCi,nrows,ncols);
  Zeros(CinvdCj,nrows,ncols);
  Zeros(Product,nrows,ncols);

  for (int i=0;i<nmatrices;i++)
    {
      dCi=dcovfull(IR(0,nrows),IR(i*nmatrices,i*nmatrices+ncols));
      Gemm(NORMAL,NORMAL,Real(1),Cinv,dCi,Real(0),CinvdCi);                                                                                           
      for (int j=0;j<nmatrices;j++)
        {
	  dCj=dcovfull(IR(0,nrows),IR(j*nmatrices,j*nmatrices+ncols));
          Gemm(NORMAL,NORMAL,Real(1),Cinv,dCj,Real(0),CinvdCj);
          Gemm(NORMAL,NORMAL,Real(1),CinvdCi,CinvdCj,Real(0),Product);                                                                                    
          FisherElement=Trace(Product)*0.5;
          FisherMatrix.Set(i,j,FisherElement);
	  Output("FisherMatrix(",i,",",j,")=",FisherElement);
        }
    }
}

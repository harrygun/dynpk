#include "El.hpp"
#include "../include/VectorFunctions.hpp"
using namespace El;
typedef double Real;


void GetAutocovarianceMatrix(DistMatrix<double>& v_data, DistMatrix<double>&  CovarMatrix)
{
	// take read matrix (column matrix) to serve as vector, 
	// output covariance matrix CovarMatrix    
	    									
	int n = v_data.Height();
	vector<double> v_data_vector(n),v_covariance(n),v_covar_full(2*n-1);

	DistMatrixToVector(v_data, v_data_vector);
	Autocovariance1DReal(v_data_vector,v_covariance);
	GetAutocovarianceVector(v_covariance,v_covar_full);
	Toeplitz(CovarMatrix,n,n,v_covar_full);
}

void GetAveragedAutocovarianceMatrix(DistMatrix<double>& M_data, DistMatrix<double>&  CovarMatrix)
{                                                                                                                                           
  int n = M_data.Height();
  DistMatrix<double> sum(n,1),matcovar(n,1);
  Zeros(sum,n,1);
  Zeros(matcovar,n,1);
  vector<double> meancovar(n),v_data_vector(n),v_covariance(n),v_covar_full(2*n-1);
  
  // method 1                                    
  for (int i=0;i<n;i++)
    {
      GetRow(M_data,i,v_data_vector);
      Autocovariance1DReal(v_data_vector,v_covariance);
      VectorToDistMatrix(v_covariance,matcovar);
      sum+=matcovar;
    }
  sum*=(1/double(n));
  DistMatrixToVector(sum,meancovar);
  GetAutocovarianceVector(meancovar,v_covar_full);
  Toeplitz(CovarMatrix,n,n,v_covar_full);
}



void FillBlockToeplitz(DistMatrix<double>& A,DistMatrix<double>& X)
{                                                   
	// function overwrites A, which becomes block Toeplitz                              
	// each generating vector for each submatrix C, Cdiag, is a row of X                
	int nsubmatrices = X.Height(); 
	int blocksize = X.Width();
	int nskip = X.Width();   // submatrices same size as row vectors of X
	int sizeA = A.Height();
	int i,j,k,iloc,jloc;
	DistMatrix<double> C,Cdiag;
	vector<double> vrow(blocksize);
	vector<double> vrow2(blocksize*2-1);
	
	Zeros(C,blocksize,blocksize);
	Zeros(Cdiag,blocksize,blocksize);
	
	iloc=0;
	jloc=0;

	// handle diagonal matrices separately         
	GetRow(X,0,vrow);
	GetCirculantGeneratingVector(vrow,vrow2);
	Toeplitz(Cdiag,blocksize,blocksize,vrow2);
		
	for (k=0;k<sizeA;k+=nskip)
	{
		A(IR(k,k+blocksize),IR(k,k+blocksize))=Cdiag;
	}
	
	// fill lower off-diagonal blocks              
	int rownum=0;
	for (i=nskip;i<pow(nskip,2);i+=nskip)
	{
		rownum+=1;
		for (j=0;j<(pow(nskip,2)-i);j+=nskip)
		{
			GetRow(X,rownum,vrow);
			GetCirculantGeneratingVector(vrow,vrow2);
			Toeplitz(C,blocksize,blocksize,vrow2);
			A(IR(i+j,i+j+blocksize),IR(j,j+blocksize))=C;
		}
	} 

	// fill upper off-diagonal blocks                                         
	rownum=nskip;
	for (i=nskip;i<pow(nskip,2);i+=nskip)
	{
		rownum-=1;
		for (j=0;j<(pow(nskip,2)-i);j+=nskip)
		  {
			GetRow(X,rownum,vrow);
			GetCirculantGeneratingVector(vrow,vrow2);
			Toeplitz(C,blocksize,blocksize,vrow2);
			A(IR(j,j+blocksize),IR(j+i,j+i+blocksize))=C;
		}
	}
}

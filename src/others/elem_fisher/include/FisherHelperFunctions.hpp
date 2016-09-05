#ifndef FISHERHELPERFUNCTIONS_H
#define FISHERHELPERFUNCTIONS_H

#include <El.hpp>
using namespace El;

void mpixel_idx(int pixel,int timedimension, int freqdimension, vector<int>& idx);
double GetCovarianceMatrix(int pixeldimension, int submatrixlength, int dx0, int dx1, DistMatrix<double>& bandpower, DistMatrix<double>& dcov);
void GetSymmetricCirculantVector(vector<double>& v, vector<double>&v_out);
void Combine_dcovi(std::string basename, int nmatrices, int nrows, int ncols, DistMatrix<double>& full_dcovi);
void GetFisherMatrix(DistMatrix<double>& dcovfull, int nmatrices, DistMatrix<double>& Cov, DistMatrix<double>& FisherMatrix);
#endif

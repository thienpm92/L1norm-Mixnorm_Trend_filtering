#ifndef MIXLP_FILTER_H
#define MIXLP_FILTER_H
#include <armadillo>
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/CholmodSupport>
#include <Eigen/SuperLUSupport>
#include <Eigen/SPQRSupport>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace arma;
class MixLpfilter{
    private:
        const double MAXITER = 20;   /* IPM parameter: max iter. of IPM */
        double lambda1;
        double lambda2;
        int size;
        const double ABSTOL = 1e-4;
        const double RELTOL = 1e-2;
        const double alpha = 1.0;
        double rho = 10.0;
        double rho2 = 10.0;
        arma::vec y_vec;
        arma::mat D1;
        arma::mat D2;

        arma::mat soft_threshold(arma::mat& x,double k);
    public:
        MixLpfilter(const double* y, double lamda1, double lamda2, const int n, double ro);
        ~MixLpfilter();
        void Mxfilter_L21(double* output);
        void Mxfilter_L12(double* output);
        void Mxfilter_L2(double *x);
        void Mxfilter_L1(double *x);

};

#endif // MIXLP_FILTER_H

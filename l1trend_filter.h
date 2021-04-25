#ifndef L1TREND_FILTER_H
#define L1TREND_FILTER_H
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
class L1filter{
    private:
        const double ALPHA = 0.01; /* linesearch parameter (0,0.5] */
        const double BETA = 0.5;  /* linesearch parameter (0,1) */
        const double MU = 2;    /* IPM parameter: t update */
        const double MAXITER = 20;   /* IPM parameter: max iter. of IPM */
        const double MAXLSITER = 20;   /* IPM parameter: max iter. of linesearch */
        const double TOL = 1e-4; /* IPM parameter: tolerance */

        arma::sp_mat D1;
        arma::sp_mat D2;
        arma::vec y_vec;
        double lambda1;
        double lambda2;
        int size;
    public:
        L1filter(const double *y,const int n);
        ~L1filter();
        void l1tf(const double lambda, double *x, int mode);
        void lambdamax1( double* lmax);
        void lambdamax2(double* lmax);
        void l1tf_mixed(const double lambda1,const double lambda2, double *x);
        void HP_filter_1(const int n, const double *y, const double lambda, double *x);
        void HP_filter_2(const int n, const double *y, const double lambda, double *x);
};



#endif // L1TREND_FILTER_H

#ifndef L1TREND_ADMM_H
#define L1TREND_ADMM_H
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

class l1_admm
{
    private:
        const double MAXITER = 20;   /* IPM parameter: max iter. of IPM */
        double lambda;
        int size;
        const double ABSTOL = 1e-4;
        const double RELTOL = 1e-2;
        const double alpha = 1.0;
        double rho=10.0;
        arma::vec y_vec;
        arma::mat D1;
        arma::mat D2;

        arma::mat soft_threshold(arma::mat& x,double k);
    public:
        l1_admm(const double* y,double lambda, const int n );
        ~l1_admm();
        void TV_lasso(double* output, int m);

};

#endif // L1TREND_ADMM_H

#include <QCoreApplication>
#include "l1trend_filter.h"
#include "l1trend_admm.h"
#include "mixlp_filter.h"
using namespace std;
using namespace arma;

void HP_filter_1(const int n, const double *y, const double lambda, double *x)
{
    const int    m = n - 1;  /* length of Dx */
    vec y_vec = vec(y, n);

    mat I2 = eye(m, m);
    mat O2 = zeros(m, 1);
    sp_mat D = sp_mat(join_horiz(-1.0*I2, O2)+ join_horiz(O2,I2));
    sp_mat DTD = D.t()*D*2*lambda;
    sp_mat I = speye<sp_mat>(n,n);
    mat invD(I+DTD);
    vec x_vec=inv(invD)*y_vec;
    ::memcpy(x, x_vec.memptr(), y_vec.n_elem * sizeof(double));
    return;
};

void HP_filter(const int n, const double *y, const double lambda1,const double lambda2, double *x)
{
    vec y_vec = vec(y, n);
    const int    m = n - 1;  /* length of Dx */
    const int    p = n - 2;  /* length of Dx */
    sp_mat Unit = speye<sp_mat>(n,n);
    arma::mat I1 = eye(m, m);
    arma::mat O1 = zeros(m, 1);
    sp_mat D1 = arma::sp_mat(join_horiz(-1.0*I1, O1)+ join_horiz(O1,I1));

    arma::mat I2 = eye(p, p);
    arma::mat O2 = zeros(p, 1);
    sp_mat D2 = arma::sp_mat(join_horiz(I2, join_horiz(O2, O2)) + join_horiz(O2, join_horiz(-2.0 * I2, O2)) + join_horiz(O2, join_horiz(O2, I2)));

    sp_mat D2xD2 = D2.t()*D2*lambda2;
    sp_mat D1xD1 = D1.t()*D1*lambda1;
    mat invD(D1xD1+D2xD2+Unit);
    vec x_vec=inv(invD)*y_vec;
    ::memcpy(x, x_vec.memptr(), y_vec.n_elem * sizeof(double));
    return;
}

int main(int argc, char *argv[])
{
    arma::vec y;
    y.load("/home/eric/qt_project/Mixed_filtering/videodata.txt");
    arma::vec x(y.n_rows);
    double lambda1 = 50;
    double lambda2 = 100;
    double rho=10;
    int mode=1;
    wall_clock w;
    w.tic();

//    L1filter filter( y.memptr(),y.n_rows);
////    filter.l1tf(lambda1, x.memptr(), mode);
//    filter.l1tf_mixed(lambda1,lambda2, x.memptr());
//    cout << "Time to L1 Trend filter " << w.toc() << endl;
//    arma::mat r = join_horiz(y, x);
//    r.eval().save("result.csv", csv_ascii);


//   l1_admm test(y.memptr(),lambda1, y.n_rows);
//   test.TV_lasso(x.memptr(),mode);
//   arma::mat r = join_horiz(y, x);
//   r.eval().save("result_admm.csv", csv_ascii);

   MixLpfilter test(y.memptr(),lambda1,lambda2, y.n_rows,rho);
   test.Mxfilter_L21(x.memptr());
   arma::mat r = join_horiz(y, x);
   r.eval().save("result_mxlp.csv", csv_ascii);

//   HP_filter_1(y.n_rows, y.memptr(),lambda1, x.memptr());
//   arma::mat r = join_horiz(y, x);
//   r.eval().save("result_L2.csv", csv_ascii);



   return 0;
}

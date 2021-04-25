#include"mixlp_filter.h"
#include <armadillo>
#include <limits>

MixLpfilter::MixLpfilter(const double *y,double lamda1, double lamda2, const int n,double ro){
    y_vec=vec(y,n);
    size = n;
    lambda1=lamda1;
    lambda2=lamda2;
    rho=ro;
    const int    m = n - 1;  /* length of Dx */
    const int    p = n - 2;  /* length of Dx */
    arma::mat I1 = eye(m, m);
    arma::mat O1 = zeros(m, 1);
    arma::sp_mat temp1 = arma::sp_mat(join_horiz(-1.0*I1, O1)+ join_horiz(O1,I1));
    D1 = conv_to<mat>::from(temp1);

    arma::mat I2 = eye(p, p);
    arma::mat O2 = zeros(p, 1);
    arma::sp_mat temp2 = arma::sp_mat(join_horiz(I2, join_horiz(O2, O2)) + join_horiz(O2, join_horiz(-2.0 * I2, O2)) + join_horiz(O2, join_horiz(O2, I2)));
    D2 = conv_to<mat>::from(temp2);
};

MixLpfilter::~MixLpfilter()
{

};

arma::mat MixLpfilter::soft_threshold(arma::mat& x, double k)
{
    double vi;
    arma::mat z=zeros(x.n_rows, 1);

    for(int i=0;i<x.n_rows;i++){
        vi=x.at(i,0);
        if(vi>k)
            z.at(i,0)=vi-k;
        else if(vi<-k)
            z.at(i,0)=vi+k;
        else
            z.at(i,0)=0;
    }
    return z;
};

void MixLpfilter::Mxfilter_L21(double *output)
{
    int m=size-1;
    int p=size-2;

    arma::mat z = zeros(p, 1);
    arma::mat zold = zeros(p, 1);
    arma::mat u = zeros(p, 1);
    arma::mat x = zeros(size, 1);
    arma::mat I = eye<arma::mat>(size,size);
    arma::mat DTD1 = D1.t()*D1;
    arma::mat DTD2 = D2.t()*D2;
    int iter = 0;
    while (iter < MAXITER) {
        //x-update
        arma::mat temp1 = arma::inv((I + lambda1*DTD1 + rho*DTD2));
        x=temp1*(y_vec + rho*D2.t()*(z-u));

        //z-update with relaxation
        zold=z;
        arma::mat A_hat = alpha*D2*x + (1-alpha)*zold ;
        arma::mat temp2 = A_hat + u;
        z = soft_threshold(temp2,lambda2/rho);

        //u-update
        u = u + A_hat - z;

        /* Termination checks */
        double r_norm = norm(D2*x - z);
        double s_norm = norm(-rho*D2.t()*(z-zold));
        double eps_pri = sqrt(size)*ABSTOL + RELTOL*std::max(norm(D2*x),norm(-z));
        double eps_dual = sqrt(size)*ABSTOL + RELTOL*norm(rho*D2.t()*u);

        if(r_norm<=eps_pri && s_norm<=eps_dual)
            break;
        cout<<iter<<endl;
        iter++;
    }
    if (iter >= MAXITER) {
        arma::colvec x_vec = conv_to< colvec >::from(x);
        ::memcpy(output, x_vec.memptr(), y_vec.n_elem * sizeof(double));
        return;
    }
};

void MixLpfilter::Mxfilter_L12(double *output)
{
    int m=size-2;
    int p=size-1;

    arma::mat z = zeros(p, 1);
    arma::mat zold = zeros(p, 1);
    arma::mat u = zeros(p, 1);
    arma::mat x = zeros(size, 1);
    arma::mat I = eye<arma::mat>(size,size);
    arma::mat DTD1 = D2.t()*D2;
    arma::mat DTD2 = D1.t()*D1;
    int iter = 0;
    while (iter < MAXITER) {
        //x-update
        arma::mat temp1 = arma::inv((I + lambda1*DTD1 + rho*DTD2));
        x=temp1*(y_vec + rho*D1.t()*(z-u));

        //z-update with relaxation
        zold=z;
        arma::mat A_hat = alpha*D1*x + (1-alpha)*zold ;
        arma::mat temp2 = A_hat + u;
        z = soft_threshold(temp2,lambda2/rho);

        //u-update
        u = u + A_hat - z;

        /* Termination checks */
        double r_norm = norm(D1*x - z);
        double s_norm = norm(-rho*D1.t()*(z-zold));
        double eps_pri = sqrt(size)*ABSTOL + RELTOL*std::max(norm(D1*x),norm(-z));
        double eps_dual = sqrt(size)*ABSTOL + RELTOL*norm(rho*D1.t()*u);

        if(r_norm<=eps_pri && s_norm<=eps_dual)
            break;
        cout<<iter<<endl;
        iter++;
    }
    if (iter >= MAXITER) {
        arma::colvec x_vec = conv_to< colvec >::from(x);
        ::memcpy(output, x_vec.memptr(), y_vec.n_elem * sizeof(double));
        return;
    }
};

void MixLpfilter::Mxfilter_L2(double *x)
{
    const int    m = size - 1;  /* length of Dx */
    const int    p = size - 2;  /* length of Dx */
    mat Unit = eye<mat>(size,size);
    mat D2xD2 = D2.t()*D2*lambda2;
    mat D1xD1 = D1.t()*D1*lambda1;
    mat invD(D1xD1+D2xD2+Unit);
    vec x_vec=inv(invD)*y_vec;
    ::memcpy(x, x_vec.memptr(), y_vec.n_elem * sizeof(double));
    return;
}

void MixLpfilter::Mxfilter_L1(double *output)
{
    int m=size-1;
    int p=size-2;

    arma::mat z1 = zeros(m, 1);
    arma::mat z1old = zeros(m, 1);
    arma::mat z2 = zeros(p, 1);
    arma::mat z2old = zeros(p, 1);
    arma::mat u1 = zeros(m, 1);
    arma::mat u2 = zeros(p, 1);
    arma::mat x = zeros(size, 1);
    arma::mat I = eye<arma::mat>(size,size);
    arma::mat DTD1 = D1.t()*D1;
    arma::mat DTD2 = D2.t()*D2;
    int iter = 0;
    while (iter < MAXITER) {
        //x-update
        arma::mat temp = arma::inv((I + rho*DTD1 + rho2*DTD2));
        x=temp*(y_vec + rho*D1.t()*(z1-u1) + rho2*D2.t()*(z2-u2));

        //z-update with relaxation
        z1old=z1;
        z2old=z2;
        arma::mat A_hat1 = alpha*D1*x + (1-alpha)*z1old ;
        arma::mat temp1 = A_hat1 + u1;
        arma::mat A_hat2 = alpha*D2*x + (1-alpha)*z2old ;
        arma::mat temp2 = A_hat2 + u2;
        z1 = soft_threshold(temp1,lambda1/rho);
        z2 = soft_threshold(temp2,lambda2/rho2);
        //u-update
        u1 = u1 + A_hat1 - z1;
        u2 = u2 + A_hat2 - z2;
        /* Termination checks */
        double r1_norm = norm(D1*x - z1);
        double s1_norm = norm(-rho*D1.t()*(z1-z1old));
        double eps1_pri = sqrt(size)*ABSTOL + RELTOL*std::max(norm(D1*x),norm(-z1));
        double eps1_dual = sqrt(size)*ABSTOL + RELTOL*norm(rho*D1.t()*u1);

        double r2_norm = norm(D2*x - z2);
        double s2_norm = norm(-rho2*D2.t()*(z2-z2old));
        double eps2_pri = sqrt(size)*ABSTOL + RELTOL*std::max(norm(D2*x),norm(-z2));
        double eps2_dual = sqrt(size)*ABSTOL + RELTOL*norm(rho2*D2.t()*u2);

        if(r1_norm<=eps1_pri && s1_norm<=eps1_dual && r2_norm<=eps2_pri && s2_norm<=eps2_dual)
            break;
        cout<<iter<<endl;
        iter++;
    }
    if (iter >= MAXITER) {
        arma::colvec x_vec = conv_to< colvec >::from(x);
        ::memcpy(output, x_vec.memptr(), y_vec.n_elem * sizeof(double));
        return;
    }
};

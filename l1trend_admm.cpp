#include "l1trend_admm.h"

l1_admm:: l1_admm(const double *y,const double lamda, const int n)
{
   y_vec=vec(y,n);
   size = n;
   lambda=lamda;
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

l1_admm::~l1_admm()
{

}

arma::mat l1_admm::soft_threshold(arma::mat& x, double k)
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

void l1_admm::TV_lasso(double *output,int mode)
{
    int m;
    arma::mat D;
    if(mode==0){
        m=size-1;
        D=D1;
    }
    else{
        m=size-2;
        D=D2;
    }
    arma::mat z = zeros(m, 1);
    arma::mat zold = zeros(m, 1);
    arma::mat u = zeros(m, 1);
    arma::mat x = zeros(size, 1);
    arma::mat I = eye<arma::mat>(size,size);
    arma::mat DTD = D.t()*D;

    int iter = 0;
    while (iter < MAXITER) {
        //x-update
        arma::mat temp1 = arma::inv((I + rho*DTD));
        x=temp1*(y_vec + rho*D.t()*(z-u));

        //z-update with relaxation
        zold=z;
        arma::mat A_hat = alpha*D*x + (1-alpha)*zold ;
        arma::mat temp2 = A_hat + u;
        z = soft_threshold(temp2,lambda/rho);

        //u-update
        u = u + A_hat - z;

        /* Termination checks */
        double r_norm = norm(D*x - z);
        double s_norm = norm(-rho*D.t()*(z-zold));
        double eps_pri = sqrt(size)*ABSTOL + RELTOL*std::max(norm(D*x),norm(-z));
        double eps_dual = sqrt(size)*ABSTOL + RELTOL*norm(rho*D.t()*u);

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

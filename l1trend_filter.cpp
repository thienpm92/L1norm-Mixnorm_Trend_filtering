#include"l1trend_filter.h"
#include <armadillo>
#include <limits>

L1filter::L1filter(const double *y,const int n){
  y_vec = vec(y, n);
  size = n;
  const int    m = n - 1;  /* length of Dx */
  const int    p = n - 2;  /* length of Dx */
  arma::mat I1 = eye(m, m);
  arma::mat O1 = zeros(m, 1);
  D1 = arma::sp_mat(join_horiz(-1.0*I1, O1)+ join_horiz(O1,I1));

  arma::mat I2 = eye(p, p);
  arma::mat O2 = zeros(p, 1);
  D2 = arma::sp_mat(join_horiz(I2, join_horiz(O2, O2)) + join_horiz(O2, join_horiz(-2.0 * I2, O2)) + join_horiz(O2, join_horiz(O2, I2)));

};

L1filter::~L1filter()
{

}

void L1filter::l1tf(const double lambda, double *x,int mode)
{

    double t = 1e-10;
    double step = std::numeric_limits<double>::infinity();
    double dobj = 0.0;
    unsigned int iter = 0;
    int m;
    int n=size;
    arma::sp_mat D;
    arma::sp_mat DDT;
    arma::mat Dy;
    if(mode==0){
        m = size-1;
        D=D1;
        DDT = D * D.t();
        Dy = D * y_vec;
    }
    else{
        m = size-2;
        D=D2;
        DDT = D * D.t();
        Dy = D * y_vec;
    }


    arma::mat z = zeros(m, 1);
    arma::mat mu1 = ones(m, 1);
    arma::mat mu2 = ones(m, 1);

    arma::mat f1 = z - lambda;
    arma::mat f2 = -z - lambda;
    arma::mat DTz(n, 1);
    arma::mat DDTz(m, 1);
    arma::mat w(m, 1);
    arma::mat rz(m, 1);
    arma::sp_mat S(m, m);
    arma::mat r(m, 1);
    arma::mat dz(m, 1);
    arma::mat dmu1(m, 1);
    arma::mat dmu2(m, 1);
    arma::mat resDual(m, 1);
    arma::mat newResDual(m, 1);
    arma::mat resCent(2 * m, 1);
    arma::mat newresCent(2 * m, 1);
    arma::mat residual(3 * m, 1);
    arma::mat newResidual(3 * m, 1);
    arma::mat newz(m, 1);
    arma::mat newmu1(m, 1);
    arma::mat newmu2(m, 1);
    arma::mat newf1(m, 1);
    arma::mat newf2(m, 1);
    for (; iter < MAXITER; ++iter)
    {
        DTz = (z.t() * D).t();
        DDTz = D * DTz;
        w = Dy - (mu1 - mu2);

        // two ways to evaluate primal objective :
        // 1) using dual variable of dual problem
        // 2) using optimality condition
        arma::vec xw = spsolve(DDT, w);

        arma::mat pobj1 = (0.5 * w.t() * (xw)) + lambda * arma::sum(mu1 + mu2);

        arma::mat pobj2 = ((0.5 * DTz.t() * DTz)) + lambda * arma::sum(abs(Dy - DDTz));
        arma::mat pobjm = arma::min(pobj1, pobj2);
        double pobj = pobjm.at(0, 0);
        dobj = std::max((-0.5 * DTz.t() * DTz + Dy.t() * z)[0,0], dobj);
        double gap = pobj - dobj;

        //Stopping criteria
        if (gap <= TOL)
        {
            arma::vec x_vec = y_vec - D.t() * z;
            ::memcpy(x, x_vec.memptr(), sizeof(double)* y_vec.n_elem);
            return;
        }

        if (step >= 0.2)
        {
            t = std::max(2.0 * m * MU/gap, 1.2 * t);
        }

        // Calculate Newton Step
        rz = DDTz - w;
        S = DDT - diagmat(mu1/f1 + mu2/f2);
        r = -DDTz + Dy + ((1 / t) / f1) - ((1 / t) / f2);
        dz = arma::mat(spsolve(S, r));
        dmu1 = -(mu1 + ((dz % mu1) + (1 / t)) / f1);
        dmu2 = -(mu2 + ((dz % mu2) + (1 / t)) / f2);

        resDual = rz;
        resCent =  arma::join_vert((-mu1 % f1) - 1 / t, (-mu2 % f2) - 1 / t);
        residual =  arma::join_vert(resDual, resCent);

        // Backtracking linesearch.
        arma::umat  negIdx1 = all(dmu1 < 0.0);
        arma::umat negIdx2 = all(dmu2 < 0.0);
        step = 1.0;

        if (any(vectorise(negIdx1))){
            step = std::min(step, 0.99*arma::min(-mu1(negIdx1)/dmu1(negIdx1)));
        }

        if (any(vectorise(negIdx2)))
        {
            step = std::min(step, 0.99*arma::min(-mu2(negIdx2)/ dmu2(negIdx2)));
        }

        for (unsigned int liter = 0; liter < MAXLSITER; ++liter)
        {
            newz = z + step * dz;
            newmu1 = mu1 + step * dmu1;
            newmu2 = mu2 + step * dmu2;
            newf1 = newz - lambda;
            newf2 = -newz - lambda;

            // Update residual

            //% UPDATE RESIDUAL
            newResDual = DDT * newz - Dy + newmu1 - newmu2;
            newresCent = arma::join_vert((-newmu1 % newf1) - 1 / t, (-newmu2 % newf2) - 1 / t);
            newResidual = arma::join_vert(newResDual, newresCent);

            if ((std::max(arma::max(vectorise(newf1)), arma::max(vectorise(newf2))) < 0.0) && norm(newResidual) <= (1 - ALPHA*step)*norm(residual)) {
                break;
            }

            step = BETA * step;
        }
        z = newz; mu1 = newmu1; mu2 = newmu2; f1 = newf1; f2 = newf2;
    }

    // The solution may be close at this point, but does not meet the stopping
    // criterion(in terms of duality gap).

    if (iter >= MAXITER) {
        arma::vec x_vec = y_vec - D.t() *z;
        ::memcpy(x, x_vec.memptr(), y_vec.n_elem * sizeof(double));
        return;
    }
};

void L1filter::l1tf_mixed(const double lambda1,const double lambda2, double *x)
{
    /* dimension */
    const int n=size;
    const int m = size - 1;  /* length of Dx */
    const int p = size - 2;  /* length of Dx */
    double t = 1e-10;
    double step = std::numeric_limits<double>::infinity();
    double dobj = 0.0;
    unsigned int iter = 0;

    arma::sp_mat D= join_vert( D1, D2);
    arma::sp_mat DDT = D * D.t();
    arma::mat Dy = D * y_vec;

    arma::mat z = zeros(m+p, 1);
    arma::mat z_d1 = zeros(m, 1);
    arma::mat z_d2 = zeros(p, 1);
    arma::mat f1_d1 = z_d1 - lambda1;
    arma::mat f2_d1 = -z_d1 - lambda1;
    arma::mat f1_d2 = z_d2 - lambda2;
    arma::mat f2_d2 = -z_d2 - lambda2;
    arma::mat f1 = join_vert( f1_d1, f1_d2 );
    arma::mat f2 = join_vert( f2_d1, f2_d2 );

    arma::mat mu1_d1 = lambda1*ones(m, 1);
    arma::mat mu2_d1 = lambda1*ones(m, 1);
    arma::mat mu1_d2 = lambda2*ones(p, 1);
    arma::mat mu2_d2 = lambda2*ones(p, 1);
    arma::mat mu1 = join_vert(mu1_d1, mu1_d2);
    arma::mat mu2 = join_vert(mu2_d1, mu2_d2);
    arma::mat lambda = join_vert(mu1_d1, mu1_d2);

    arma::mat DTz(n, 1);
    arma::mat DDTz(m+p, 1);
    arma::mat w(m+p,1);
    arma::mat rz(m+p,1);
    arma::sp_mat S(m+p, m+p);
    arma::mat r(m+p, 1);
    arma::mat dz(m+p, 1);
    arma::mat dmu1(m+p, 1);
    arma::mat dmu2(m+p, 1);
    arma::mat resDual(m+p, 1);
    arma::mat newResDual(m+p, 1);
    arma::mat resCent(2 * (m+p), 1);
    arma::mat newresCent(2 * (m+p), 1);
    arma::mat residual(3 * (m+p), 1);
    arma::mat newResidual(3 * (m+p), 1);
    arma::mat newz(m+p, 1);
    arma::mat newmu1(m+p, 1);
    arma::mat newmu2(m+p, 1);
    arma::mat newf1(m+p, 1);
    arma::mat newf2(m+p, 1);
    arma::mat xw(m+p,1);

    for (; iter < MAXITER; ++iter)
    {
        cout<<iter<<endl;
        DTz = (z.t() * D).t();
        DDTz = D * DTz;
        w = Dy - (mu1 - mu2);

        // two ways to evaluate primal objective :
        // 1) using dual variable of dual problem
        // 2) using optimality condition
        mat tempDDT = conv_to<mat>::from(DDT);
        arma::mat invDDT = pinv(tempDDT);
        vec xw = invDDT* w;
//        vec xw = spsolve(DDT, w);

        mat pobj1 = (0.5 * w.t() * (xw)) + arma::sum(mu1 + mu2);

        arma::mat pobj2 = ((0.5 * DTz.t() * DTz)) + arma::sum(abs(Dy - DDTz));
        arma::mat pobjm =  pobj2;
        double pobj = pobjm.at(0, 0);
        dobj = std::max((-0.5 * DTz.t() * DTz + Dy.t() * z)[0,0], dobj);
        double gap = pobj - dobj;

    //        Stopping criteria
        if (gap <= TOL)
        {
            arma::vec x_vec = y_vec - D.t() * z;
            ::memcpy(x, x_vec.memptr(), sizeof(double)* y_vec.n_elem);
            return;
        }
        if (step >= 0.2)
            t = std::max(2.0 * (m+p) * MU/gap, 1.2 * t);

        // Calculate Newton Step
        rz = DDTz - w;
        S = DDT - diagmat(mu1/f1 + mu2/f2);
        r = -DDTz + Dy + ((1 / t) / f1) - ((1 / t) / f2);

        dz = arma::mat(spsolve(S, r));


        dmu1 = -(mu1 + ((dz % mu1) + (1 / t)) / f1);
        dmu2 = -(mu2 + ((dz % mu2) + (1 / t)) / f2);

        resDual = rz;
        resCent =  join_vert((-mu1 % f1) - 1 / t, (-mu2 % f2) - 1 / t);
        residual =  join_vert(resDual, resCent);

        // Backtracking linesearch.
        arma::umat  negIdx1 = all(dmu1 < 0.0);
        arma::umat negIdx2 = all(dmu2 < 0.0);
        step = 1.0;

        if (any(vectorise(negIdx1))){
            step = std::min(step, 0.99*arma::min(-mu1(negIdx1)/dmu1(negIdx1)));
        }

        if (any(vectorise(negIdx2)))
        {
            step = std::min(step, 0.99*arma::min(-mu2(negIdx2)/ dmu2(negIdx2)));
        }

        for (unsigned int liter = 0; liter < MAXLSITER; ++liter)
        {
            newz = z + step * dz;
            newmu1 = mu1 + step * dmu1;
            newmu2 = mu2 + step * dmu2;
            newf1 = newz - lambda;
            newf2 = -newz - lambda;

            // Update residual

            //% UPDATE RESIDUAL
            newResDual = DDT * newz - Dy + newmu1 - newmu2;
            newresCent = join_vert((-newmu1 % newf1) - 1 / t, (-newmu2 % newf2) - 1 / t);
            newResidual = join_vert(newResDual, newresCent);

            if ((std::max(arma::max(vectorise(newf1)), arma::max(vectorise(newf2))) < 0.0) && norm(newResidual) <= (1 - ALPHA*step)*norm(residual)) {
                break;
            }

            step = BETA * step;
        }
        z = newz; mu1 = newmu1; mu2 = newmu2; f1 = newf1; f2 = newf2;
    }

    // The solution may be close at this point, but does not meet the stopping
    // criterion(in terms of duality gap).

    if (iter >= MAXITER) {
        arma::vec x_vec = y_vec - D.t() *z;
        ::memcpy(x, x_vec.memptr(), y_vec.n_elem * sizeof(double));
        return;
    }
};

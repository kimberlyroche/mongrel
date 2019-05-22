#include <Rcpp.h>
#include <RcppEigen.h>
#include "MatDist.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::ArrayXXd;
using Eigen::Map;

// these should be passed as references

// explicit references to standard namespace here are motivated by the observation that
// these calls don't properly resolve on some Windows installs
// [[Rcpp::export]]
Eigen::MatrixXd power_G(Eigen::MatrixXd G, 
                        int it_begin,
                        int it_end){
  if(it_begin == it_end) { return(G); }
  int p = G.rows();
  int power_it = std::abs(it_end - it_begin) + 1;
  // powering must explicitly be carried out on complex components (if they exist) using
  // complex exponents; doesn't look like these a vector wrapper for this; perform
  // element-wise powering by hand
  std::complex<double> power_complex(power_it, 0);
  Eigen::ComplexEigenSolver<MatrixXcd> Gdec(G);
  VectorXcd powered_eigenvalues = Gdec.eigenvalues();
  for(int j=0; j<p; j++) {
    powered_eigenvalues(j) = std::pow(powered_eigenvalues(j), power_complex);
  }
  MatrixXcd Lambda = MatrixXcd::Zero(p,p);
  // use of .diagonal() in assignment doesn't work here for type reasons
  // troubleshoot this
  for(int j=0; j<p; j++) {
    Lambda(j,j) = powered_eigenvalues(j,1);
  }
  MatrixXcd reconstituted = Gdec.eigenvectors()*Lambda*(Gdec.eigenvectors().inverse());
  // in R the imaginary part was occasionally vanishingly non-zero
  // doesn't seem to happen here but chop it off just in case
  return reconstituted.real();
}

// build A matrix (covariance over samples) assuming time-invariant parameters F, G, W, gamma
// [[Rcpp::export]]
Eigen::MatrixXd dlm_cov(int T, double gamma, Eigen::VectorXd F, Eigen::MatrixXd G, Eigen::MatrixXd W, Eigen::MatrixXd C0) {
  // check T >= 1
  MatrixXd res = MatrixXd::Zero(T,T);
  MatrixXd Ft = F.transpose();
  MatrixXd Gt = G.transpose();
  for(int i=0; i<T; i++) {
    for(int j=0; j<T; j++) {
      if(i==j) {
        int t = j+1;
        // diagonal
        res(j,j) += gamma;
        res(j,j) += (Ft*W*F)(0,0); // 1x1 vector product
        for(int ell=t; ell>=2; ell--) {
          res(j,j) += (Ft*power_G(G, t, ell)*W*power_G(Gt, ell, t)*F)(0,0); // 1x1 vector
        }
        res(j,j) += (Ft*power_G(G, t, 1)*C0*power_G(Gt, 1, t)*F)(0,0); // 1x1 vector
      } else {
        // off-diagonal
        // symmetric: just copy lower triangular from upper triangular
        if(i > j) {
          res(i,j) = res(j,i);
        } else {
          int t = j+1;
          int tk = i+1;
          res(i,j) += (Ft*power_G(G, t, tk+1)*W*F)(0,0);
          for(int ell=tk; ell>=2; ell--) {
            res(i,j) += (Ft*power_G(G, t, ell)*W*power_G(Gt, ell, tk)*F)(0,0);
          }
          res(i,j) += (Ft*power_G(G, t, 1)*C0*power_G(Gt, 1, tk)*F)(0,0);
        }
      }
    }
  }
  return res;
}

// returns Sigma but we need to have it return an object containing Theta, etc.
// data is N rows x D columns; observations gives the days-after-baseline for the observation
// in the ith row
// again, we're assuming time-invariant F, G, W, and gamma
// [[Rcpp::export]]
Eigen::MatrixXd filter(Eigen::MatrixXd eta, Eigen::VectorXd F, Eigen::MatrixXd G,
                       Eigen::MatrixXd W, double gamma, int upsilon, Eigen::MatrixXd Xi,
                       Eigen::MatrixXd M0, Eigen::MatrixXd C0, Eigen::VectorXi observations) {
  //int D = eta.cols();
  //int system_dim = G.cols();
  int T = observations.maxCoeff(); // we're presuming this exists
  // init the objects we'll iteratively overwrite where they've got an initial value
  int upsilon_t = upsilon;
  Eigen::MatrixXd Gt = G.transpose();
  Eigen::VectorXd Ft = F.transpose();
  Eigen::MatrixXd Xi_t = Xi;
  Eigen::MatrixXd M_t = M0;
  Eigen::MatrixXd C_t = C0;
  // instantiate the others
  Eigen::MatrixXd A_t;
  Eigen::MatrixXd R_t;
  Eigen::VectorXd ft_t;
  double q_t;
  Eigen::VectorXd et_t;
  Eigen::MatrixXd S_t;
  for(int t=1; t<T; t++) {
    // system prior at t
    A_t = G*M_t;
    R_t = G*C_t*Gt + W;
    // one-step ahead observation forecast
    ft_t = Ft*A_t;
    q_t = gamma + (Ft*R_t*F)(0,0);
    // system posterior at t
    et_t = eta.row(t-1) - ft_t;
    S_t = R_t*F/q_t;
    M_t = A_t + S_t*et_t;
    C_t = R_t - q_t*S_t*(S_t.transpose());
    upsilon_t += 1;
    Xi += ((et_t.transpose())*et_t)/q_t;
    // sample as:
    // Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
    // Thetas.t[,,t] <- rmatrixnormal(1, M.t, C.t, Sigma.t)[,,1]
  }
  MatrixXd U = rInvWishRevCholesky(upsilon_t, Xi_t);
  return(U*U.transpose());
}




















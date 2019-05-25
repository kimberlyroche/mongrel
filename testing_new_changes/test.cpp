#include <Rcpp.h>
#include <RcppEigen.h>
#include "MatDist.h"

// [[Rcpp::depends(RcppEigen)]]

// note: tensors are available via #include <unsupported/Eigen/CXX11/Tensor>
// but strictly require compilation with C++11
// usage: Eigen::Tensor<double, 3> Thetas(system_dim, D, T);

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

// [[Rcpp::export]]
Eigen::MatrixXd dlm_B(Eigen::MatrixXd F, Eigen::MatrixXd G, Eigen::MatrixXd M0, int T) {
  int D = M0.cols();
  Eigen::MatrixXd B(D, T);
  Eigen::MatrixXd alpha(1, D);
  for(int t=1; t<=T; t++) {
    // column-wise
    alpha = (F.transpose())*power_G(G, t, 1)*M0;
    B.block(0,t-1,D,1) = alpha.transpose();
  }
  return(B);
}

// build A matrix (covariance over samples) assuming time-invariant parameters F, G, W, gamma
// [[Rcpp::export]]
Eigen::MatrixXd dlm_A(int T, double gamma, Eigen::VectorXd F, Eigen::MatrixXd G, Eigen::MatrixXd W, Eigen::MatrixXd C0) {
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
List filter(Eigen::MatrixXd eta, Eigen::MatrixXd F, Eigen::MatrixXd G,
                       Eigen::MatrixXd W, double gamma, int upsilon, Eigen::MatrixXd Xi,
                       Eigen::MatrixXd M0, Eigen::MatrixXd C0, Eigen::VectorXi observations) {
  List out(6);
  out.names() = CharacterVector::create("Thetas", "upsilon.T", "Xi.T", "Rs", "Ms", "Cs");
  int D = eta.cols();
  int system_dim = G.cols();
  int T = observations.maxCoeff();
  // return objects
  Eigen::MatrixXd Thetas(system_dim*D, T);
  Eigen::MatrixXd Rs(system_dim*system_dim, T);
  Eigen::MatrixXd Ms(system_dim*D, T);
  Eigen::MatrixXd Cs(system_dim*system_dim, T);
  // init the objects we'll iteratively overwrite where they've got an initial value
  int upsilon_t = upsilon;
  Eigen::MatrixXd Gt = G.transpose();
  Eigen::MatrixXd Ft = F.transpose();
  Eigen::MatrixXd Xi_t = Xi;
  Eigen::MatrixXd M_t = M0;
  Eigen::MatrixXd C_t = C0;
  // instantiate the others
  Eigen::MatrixXd A_t(system_dim, D);
  Eigen::MatrixXd R_t(system_dim, system_dim);
  Eigen::MatrixXd ft_t(1, D);
  double q_t;
  Eigen::MatrixXd et_t(1, D);
  Eigen::MatrixXd S_t(system_dim, 1);
  Eigen::MatrixXd LV(D, D);
  Eigen::MatrixXd Theta_t(system_dim, system_dim);
  Eigen::VectorXd res(T);
  int curr_obs_idx = -1;
  for(int t=1; t<=T; t++) {
    // find the index of this time point in the observation vector (if it exists)
    bool found = false;
    for(int k=curr_obs_idx+1; k<observations.size(); k++) {
      // this will skip all but the first observation on a given date
      if(observations(k) == t) {
        curr_obs_idx = k;
        found = true;
      }
    }
    // not found at all
    if(!found) curr_obs_idx = -1;
    // system prior at t
    A_t = G*M_t;
    R_t = G*C_t*Gt + W;
    if(curr_obs_idx < 0) {
      M_t = A_t;
      C_t = R_t;
    } else {
      // one-step ahead observation forecast
      ft_t = Ft*A_t;
      q_t = gamma + (Ft*R_t*F)(0,0);
      // system posterior at t
      et_t = eta.row(curr_obs_idx) - ft_t;
      S_t = R_t*F/q_t;
      M_t = A_t + S_t*et_t;
      C_t = R_t - q_t*S_t*(S_t.transpose());
      upsilon_t += 1;
      Xi_t += ((et_t.transpose())*et_t)/q_t;
    }
    Eigen::LLT<Eigen::MatrixXd> lltOfCt(C_t);
    Eigen::MatrixXd LU = lltOfCt.matrixL();
    // sample Sigma(t)
    LV = rInvWishRevCholesky(upsilon_t, Xi_t).matrix().transpose(); // note this is Chol. upper triangular!
    // sample Theta(t)
    Theta_t = rMatNormalCholesky(M_t, LU, LV);
    // pack up the samples to return -- R_t and C_t
    for(int i=0; i<system_dim; i++) {
      for(int j=0; j<D; j++) {
        Thetas(i+(j*system_dim),t-1) = Theta_t(i,j);
        Ms(i+(j*system_dim),t-1) = M_t(i,j);
      }
    }
    // pack up the samples to return -- Theta_t and M_t
    for(int i=0; i<system_dim; i++) {
      for(int j=0; j<system_dim; j++) {
        Rs(i+(j*system_dim),t-1) = R_t(i,j);
        Cs(i+(j*system_dim),t-1) = C_t(i,j);
      }
    }
  }
  out[0] = Thetas;
  out[1] = upsilon_t;
  out[2] = Xi_t;
  out[3] = Rs;
  out[4] = Ms;
  out[5] = Cs;
  return(out);
}

Eigen::MatrixXd unpack_sample(Eigen::MatrixXd samples, int no_rows, int no_cols, int sample_idx) {
  // add check for bad dimensions here
  Eigen::MatrixXd sample(no_rows, no_cols);
  for(int i=0; i<no_rows; i++) {
    for(int j=0; j<no_cols; j++) {
      sample(i,j) = samples(i+(j*no_rows),sample_idx);
    }
  }
  return(sample);
}

// currently only included for diagnostic purposes
// [[Rcpp::export]]
List simulation_smooth(Eigen::MatrixXd eta, Eigen::MatrixXd G, int upsilon_T, Eigen::MatrixXd Xi_T,
                       Eigen::MatrixXd Rs, Eigen::MatrixXd Ms, Eigen::MatrixXd Cs) {
  List out(1);
  out.names() = CharacterVector::create("Thetas");
  int D = eta.cols();
  int system_dim = G.cols();
  int T = Rs.cols(); // Theta is column-collapsed here
  Eigen::MatrixXd smoothed_Thetas(system_dim*D, T);
  Eigen::MatrixXd LV = rInvWishRevCholesky(upsilon_T, Xi_T).matrix().transpose(); // note this is Chol. upper triangular!
  Eigen::MatrixXd M_t = unpack_sample(Ms, system_dim, D, T-1); // M_T
  Eigen::MatrixXd C_t = unpack_sample(Cs, system_dim, system_dim, T-1); // C_T
  Eigen::LLT<Eigen::MatrixXd> lltOfCt(C_t);
  Eigen::MatrixXd LU = lltOfCt.matrixL();
  Eigen::MatrixXd smoothed_Theta_t = rMatNormalCholesky(M_t, LU, LV); // Theta_T
  for(int i=0; i<system_dim; i++) {
    for(int j=0; j<D; j++) {
      smoothed_Thetas(i+(j*system_dim),T-1) = smoothed_Theta_t(i,j);
    }
  }
  Eigen::MatrixXd R_t(system_dim, system_dim);
  Eigen::MatrixXd R_t_inv(system_dim, system_dim);
  Eigen::MatrixXd Z_t(system_dim, system_dim);
  Eigen::MatrixXd A_t(system_dim, D);
  Eigen::MatrixXd M_t_star(system_dim, D);
  Eigen::MatrixXd C_t_star(system_dim, system_dim);
  for(int t=T-1; t>0; t--) {
    // note to self: 1-indexed loop to 0-indexed data structure
    R_t = unpack_sample(Rs, system_dim, system_dim, t); // R_t+1
    R_t_inv = R_t.inverse();
    M_t = unpack_sample(Ms, system_dim, D, t-1);
    C_t = unpack_sample(Cs, system_dim, system_dim, t-1);
    Z_t = C_t*(G.transpose())*R_t_inv;
    A_t = G*M_t;
    M_t_star = M_t + Z_t*(smoothed_Theta_t - A_t);
    C_t_star = C_t - Z_t*R_t*(Z_t.transpose());
    lltOfCt = Eigen::LLT<Eigen::MatrixXd>(C_t_star); // reuse
    LU = lltOfCt.matrixL();
    smoothed_Theta_t = rMatNormalCholesky(M_t_star, LU, LV); // reuse
    for(int i=0; i<system_dim; i++) {
      for(int j=0; j<D; j++) {
        smoothed_Thetas(i+(j*system_dim),t-1) = smoothed_Theta_t(i,j);
      }
    }
  }
  out[0] = smoothed_Thetas;
  return(out);
}
















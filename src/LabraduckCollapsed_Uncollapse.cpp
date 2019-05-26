#include <stray.h>
#include <Rcpp/Benchmark/Timer.h>
#include <boost/random/mersenne_twister.hpp>

#ifdef STRAY_USE_PARALLEL
#include <omp.h>
#endif 

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXXd;
using Eigen::Map;
using Eigen::Lower;

// returns Sigma but we need to have it return an object containing Theta, etc.
// data is N rows x D columns; observations gives the days-after-baseline for the observation
// in the ith row
// again, we're assuming time-invariant F, G, W, and gamma
// no need to export
List apply_Kalman_filter(Eigen::MatrixXd eta, Eigen::MatrixXd F, Eigen::MatrixXd G,
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

// currently only included for diagnostic purposes; no need to export
List apply_simulation_smoother(Eigen::MatrixXd eta, Eigen::MatrixXd G, int upsilon_T, Eigen::MatrixXd Xi_T,
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

// TODO -- eta is passed from the optimization as a collapsed vector (I think)!
// TODO -- insert const where appropriate
// [[Rcpp::export]]
List uncollapseLabraduck(Eigen::Map<Eigen::MatrixXd> eta, // note this is essentially eta
                    Eigen::Map<Eigen::MatrixXd> F, 
                    Eigen::Map<Eigen::MatrixXd> G, 
                    Eigen::Map<Eigen::MatrixXd> W, 
                    double gamma,
                    int upsilon,
                    Eigen::Map<Eigen::MatrixXd> Xi, 
                    Eigen::Map<Eigen::MatrixXd> M0, 
                    Eigen::Map<Eigen::MatrixXd> C0, 
                    Eigen::Map<Eigen::VectorXd> observations, 
                    long seed, 
                    int ncores=-1){
  //List fit_filter = apply_Kalman_filter(eta, F, G, W, gamma, upsilon, Xi, M0, C0, observations);
  List out(2);
  out.names() = CharacterVector::create("Sigma", "Timer");
/*
  #ifdef STRAY_USE_PARALLEL
    Eigen::initParallel();
    if (ncores > 0) Eigen::setNbThreads(ncores);
    if (ncores > 0) {
      omp_set_num_threads(ncores);
    } else {
      omp_set_num_threads(omp_get_max_threads());
    }
  #endif 
  Timer timer;
  timer.step("Overall_start");
  List out(3);
  out.names() = CharacterVector::create("Lambda", "Sigma", "Timer");
  int Q = Gamma.rows();
  int D = Xi.rows()+1;
  int N = X.cols();
  int iter = eta.size()/(N*(D-1)); // assumes result is an integer !!!
  double upsilonN = upsilon + N;
  const MatrixXd GammaInv(Gamma.lu().inverse());
  const MatrixXd GammaInvN(GammaInv + X*X.transpose());
  const MatrixXd GammaN(GammaInvN.lu().inverse());
  const MatrixXd LGammaN(GammaN.llt().matrixL());
  //const Map<const MatrixXd> Eta(NULL);
  const MatrixXd ThetaGammaInvGammaN(Theta*GammaInv*GammaN);
  const MatrixXd XTGammaN(X.transpose()*GammaN);

  // Storage for output
  MatrixXd LambdaDraw0((D-1)*Q, iter);
  MatrixXd SigmaDraw0((D-1)*(D-1), iter);
  
  //iterate over all draws of eta - embarrassingly parallel with parallel rng
  #ifdef STRAY_USE_PARALLEL
    Eigen::setNbThreads(1);
    //Rcout << "thread: "<< omp_get_max_threads() << std::endl;
  #endif 
  #pragma omp parallel shared(D, N, Q, LambdaDraw0, SigmaDraw0)
  {
  #ifdef STRAY_USE_PARALLEL
    boost::random::mt19937 rng(omp_get_thread_num()+seed);
  #else 
    boost::random::mt19937 rng(seed);
  #endif 
  // storage for computation
  MatrixXd LambdaN(D-1, Q);
  MatrixXd XiN(D-1, D-1);
  MatrixXd LSigmaDraw(D-1, D-1);
  MatrixXd ELambda(D-1, Q);
  MatrixXd EEta(D-1, N);
  #pragma omp for 
  for (int i=0; i < iter; i++){
    //R_CheckUserInterrupt();
    const Map<const MatrixXd> Eta(&eta(i*N*(D-1)),D-1, N);
    LambdaN.noalias() = Eta*XTGammaN+ThetaGammaInvGammaN;
    ELambda = LambdaN-Theta;
    EEta.noalias() = Eta-LambdaN*X;
    XiN.noalias() = Xi+ EEta*EEta.transpose() + ELambda*GammaInv*ELambda.transpose();
    
    if (ret_mean){
      Map<VectorXd> LambdaNVec(LambdaN.data(), LambdaN.size());
      Map<VectorXd> XiNVec(XiN.data(), XiN.size());
      LambdaDraw0.col(i) = LambdaNVec;
      SigmaDraw0.col(i) = (upsilonN-D)*XiNVec; // mean of inverse wishart
    } else {
      // Draw Random Component
      rInvWishRevCholesky_thread_inplace(LSigmaDraw, upsilonN, XiN, rng);
      // Note: Below is valid even though LSigmaDraw is reverse cholesky factor
      Eigen::Ref<VectorXd> LambdaDraw_tmp = LambdaDraw0.col(i);
      Eigen::Map<MatrixXd> LambdaDraw(LambdaDraw_tmp.data(), D-1, Q);
      rMatNormalCholesky_thread_inplace(LambdaDraw, LambdaN, LSigmaDraw, 
                                        LGammaN.matrix(), rng);

      Eigen::Ref<VectorXd> SigmaDraw_tmp = SigmaDraw0.col(i);
      Eigen::Map<MatrixXd> SigmaDraw_tosquare(SigmaDraw_tmp.data(), D-1, D-1);
      SigmaDraw_tosquare.noalias() = LSigmaDraw*LSigmaDraw.transpose();
      
    }
  }
  }
  #ifdef STRAY_USE_PARALLEL
  if (ncores > 0){
    Eigen::setNbThreads(ncores);
  } else {
    Eigen::setNbThreads(omp_get_max_threads());  
  }
  #endif 

  IntegerVector dLambda = IntegerVector::create(D-1, Q, iter);
  IntegerVector dSigma = IntegerVector::create(D-1, D-1, iter);
  NumericVector nvLambda = wrap(LambdaDraw0);
  NumericVector nvSigma = wrap(SigmaDraw0);
  nvLambda.attr("dim") = dLambda;
  nvSigma.attr("dim") = dSigma;
  out[0] = nvLambda;
  out[1] = nvSigma;
  timer.step("Overall_stop");
  NumericVector t(timer);
  out[2] = timer;
  */
  return out;
}



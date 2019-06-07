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

// [[Rcpp::export]]
List uncollapseLabraduck(const Eigen::Map<Eigen::MatrixXd> eta, // note this is essentially eta
                    const Eigen::Map<Eigen::MatrixXd> F, 
                    const Eigen::Map<Eigen::MatrixXd> G, 
                    const Eigen::Map<Eigen::MatrixXd> W, 
                    const double W_scale,
                    const double gamma_scale,
                    const int upsilon,
                    const Eigen::Map<Eigen::MatrixXd> Xi, 
                    const Eigen::Map<Eigen::MatrixXd> M0, 
                    const Eigen::Map<Eigen::MatrixXd> C0, 
                    const Eigen::Map<Eigen::VectorXd> observations, 
                    long seed, 
                    bool ret_mean=false,
                    bool apply_smoother=false,
                    int ncores=-1){
  #ifdef STRAY_USE_PARALLEL
    Eigen::initParallel();
    if (ncores > 0) Eigen::setNbThreads(ncores);
    if (ncores > 0) {
      omp_set_num_threads(ncores);
    } else {
      omp_set_num_threads(omp_get_max_threads());
    }
  #endif 
  // we'll pull and return a single sample of Sigma for every sample of eta
  // alternatively we could use the mean of eta
  Timer timer;
  timer.step("Overall_start");
  List out(4);
  out.names() = CharacterVector::create("Sigma", "Timer", "Thetas_filtered_sample", "Thetas_smoothed_sample");
  int D = Xi.rows()+1;
  int N = observations.size();
  int T = observations.maxCoeff();
  int system_dim = G.rows();
  int iter = eta.size()/(N*(D-1)); // assumes result is an integer !!!
  MatrixXd SigmaDraw0((D-1)*(D-1), iter);
  MatrixXd ThetaFilteredDraw0(system_dim*(D-1)*T, iter);
  MatrixXd ThetaSmoothedDraw0(system_dim*(D-1)*T, iter);

  //iterate over all draws of eta - embarrassingly parallel with parallel rng
  #ifdef STRAY_USE_PARALLEL
    Eigen::setNbThreads(1);
    //Rcout << "thread: "<< omp_get_max_threads() << std::endl;
  #endif 
  #pragma omp parallel shared(D, N, SigmaDraw0, ThetaFilteredDraw0, ThetaSmoothedDraw0)
  {
  #ifdef STRAY_USE_PARALLEL
    boost::random::mt19937 rng(omp_get_thread_num()+seed);
  #else 
    boost::random::mt19937 rng(seed);
  #endif 
  MatrixXd LSigmaDraw(D-1, D-1);
  #pragma omp for 
  for (int i=0; i<iter; i++){
  //   //R_CheckUserInterrupt();
    const Map<const MatrixXd> Eta(&eta(i*N*(D-1)),D-1, N); // current sample
    TimeSeriesFit ts(F, G, W, W_scale, gamma_scale, upsilon, Xi, M0, C0, observations);
    ts.apply_Kalman_filter(Eta.transpose());
    if(ret_mean) {
      // return (1) mean of Sigma (2) sample of each t from marginal (3) sample from 1:T
      MatrixXd Sigma_mean = ts.XiT/(ts.upsilonT - D - 1);
      Eigen::Ref<VectorXd> SigmaDraw_tmp = SigmaDraw0.col(i);
      Eigen::Map<MatrixXd> SigmaDraw_tosquare(SigmaDraw_tmp.data(), D-1, D-1);
      SigmaDraw_tosquare.noalias() = Sigma_mean;
      // skip returning samples of filtered, smoothed Thetas if ret_mean=true (for now)
    } else {
      // need to reshape a 2 x D-1 x T tensor into a columns vector for inclusion as ThetaFilteredDraw0.col(i)
      // ts.Thetas_filtered is (2 x D-1, T)
      Eigen::Ref<VectorXd> ThetaFilteredDraw_tmp = ThetaFilteredDraw0.col(i); // this is 2 x D-1 x T
      Eigen::Map<MatrixXd> ThetaFDraw_tosquare(ThetaFilteredDraw_tmp.data(), system_dim*(D-1), T);
      ThetaFDraw_tosquare.noalias() = ts.Thetas_filtered;
      rInvWishRevCholesky_thread_inplace(LSigmaDraw, ts.upsilonT, ts.XiT, rng);
      Eigen::Ref<VectorXd> SigmaDraw_tmp = SigmaDraw0.col(i);
      Eigen::Map<MatrixXd> SigmaDraw_tosquare(SigmaDraw_tmp.data(), D-1, D-1);
      SigmaDraw_tosquare.noalias() = LSigmaDraw*LSigmaDraw.transpose();
      if(apply_smoother) {
        ts.apply_simulation_smoother();
        Eigen::Ref<VectorXd> ThetaSmoothedDraw_tmp = ThetaSmoothedDraw0.col(i); // this is 2 x D-1 x T
        Eigen::Map<MatrixXd> ThetaSDraw_tosquare(ThetaSmoothedDraw_tmp.data(), system_dim*(D-1), T);
        ThetaSDraw_tosquare.noalias() = ts.Thetas_simulation_smoothed;
      }
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

  IntegerVector dSigma = IntegerVector::create(D-1, D-1, iter);
  NumericVector nvSigma = wrap(SigmaDraw0);
  nvSigma.attr("dim") = dSigma;
  out[0] = nvSigma;
  timer.step("Overall_stop");
  NumericVector t(timer);
  out[1] = timer;
  out[2] = ThetaFilteredDraw0;
  out[3] = ThetaSmoothedDraw0;
  return out;
}



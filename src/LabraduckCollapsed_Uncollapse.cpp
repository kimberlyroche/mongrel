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
                    const double gamma,
                    const int upsilon,
                    const Eigen::Map<Eigen::MatrixXd> Xi, 
                    const Eigen::Map<Eigen::MatrixXd> M0, 
                    const Eigen::Map<Eigen::MatrixXd> C0, 
                    const Eigen::Map<Eigen::VectorXd> observations, 
                    long seed, 
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
  List out(3);
  out.names() = CharacterVector::create("Sigma", "Timer", "Thetas_filtered_sample");
  int D = Xi.rows()+1;
  int N = observations.size();
  int T = observations.maxCoeff();
  int iter = eta.size()/(N*(D-1)); // assumes result is an integer !!!
  MatrixXd SigmaDraw0((D-1)*(D-1), iter);
  MatrixXd ThetaFilteredDraw0;
  bool draw_taken = false;

  //iterate over all draws of eta - embarrassingly parallel with parallel rng
  #ifdef STRAY_USE_PARALLEL
    Eigen::setNbThreads(1);
    //Rcout << "thread: "<< omp_get_max_threads() << std::endl;
  #endif 
  #pragma omp parallel shared(D, N, SigmaDraw0, ThetaFilteredDraw0, draw_taken)
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
    TimeSeriesFit ts(F, G, W, gamma, upsilon, Xi, M0, C0, observations);
    ts.apply_Kalman_filter(Eta.transpose());
    if(!draw_taken) {
      // just grab the sampled marginal Thetas for the last iteration
      ThetaFilteredDraw0 = ts.Thetas_filtered;
      draw_taken = true;
    }
    rInvWishRevCholesky_thread_inplace(LSigmaDraw, ts.upsilonT, ts.XiT, rng);
    Eigen::Ref<VectorXd> SigmaDraw_tmp = SigmaDraw0.col(i);
    Eigen::Map<MatrixXd> SigmaDraw_tosquare(SigmaDraw_tmp.data(), D-1, D-1);
    SigmaDraw_tosquare.noalias() = LSigmaDraw*LSigmaDraw.transpose();
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
  return out;
}



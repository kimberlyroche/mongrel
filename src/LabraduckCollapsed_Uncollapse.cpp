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
List uncollapseLabraduck(const Eigen::Map<Eigen::VectorXd> eta, // note this is essentially eta
                    const Eigen::Map<Eigen::MatrixXd> X, 
                    const Eigen::Map<Eigen::MatrixXd> Theta,
                    const Eigen::Map<Eigen::MatrixXd> Gamma, 
                    const Eigen::Map<Eigen::MatrixXd> Xi, 
                    const double upsilon, 
                    long seed, 
                    bool ret_mean = false, 
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
  return out;
}


// #   obj <- eigen(obj)
// #   e_vec <- obj$vectors
// #   e_val <- diag(2)
// #   diag(e_val) <- obj$values
// #   if(it_begin != it_end) {
// #     if(descending) {
// #       if(it_begin > it_end) {
// #         power_it <- length(1:(it_begin-it_end+1))
// #         e_val <- e_val**power_it
// #       } else {
// #         # invalid case
// #         e_val <- matrix(0, 2, 2)
// #       }
// #     } else {
// #       if(it_begin < it_end) {
// #         power_it <- length(1:(it_end-it_begin+1))
// #         e_val <- e_val**power_it
// #       } else {
// #         # invalid case
// #         e_val <- matrix(0, 2, 2)
// #       }
// #     }
// #   }
// #   # explicitly only returning the real part of A
// #   # some tiny complex eigenvalues can be produced
// #   ret_val <- Re(e_vec%*%e_val%*%solve(e_vec)) # reconstitute
// #   return(ret_val)

// [[Rcpp::export]]
Eigen::MatrixXd power_G(Eigen::MatrixXd G, 
                        int it_begin,
                        int it_end){
  if(it_begin == it_end) {
    return(G);
  }
  int p = G.rows();
  Eigen::SelfAdjointEigenSolver<MatrixXd> Gdec(G);
  int power_it = std::abs(it_end - it_begin) + 1;
  VectorXd powered_eigenvalues = Gdec.eigenvalues().pow(power_it);
  MatrixXd Lambda = MatrixXd::Zero(p, p);
  Lambda.diagonal() = powered_eigenvalues;
  MatrixXd res = Gdec.eigenvectors()*Lambda*Gdec.eigenvectors().inverse();
  return res;
}













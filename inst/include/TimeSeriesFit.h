using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::Ref;

class TimeSeriesFit
{
  public:
    // filtering results
    MatrixXd Thetas_filtered;
    int upsilonT;
    MatrixXd XiT;
    // smoothing results to be added
    MatrixXd Thetas_simulation_smoothed;
  private:
    int N;
    int T; // T may be larger than N if there are gaps in observations
    int D;
    int system_dim;
    // these are passed in
    MatrixXd eta;
    const MatrixXd F;
    const MatrixXd G;
    const MatrixXd W;
    const double gamma;
    const int upsilon;
    const MatrixXd Xi;
    const MatrixXd M0;
    const MatrixXd C0;
    const VectorXd observations;
    // filtering results
    MatrixXd Rs;
    MatrixXd Ms;
    MatrixXd Cs;
    // indicate which applied
    bool filtered;
    bool smoothed;
    bool simulation_smoothed;

    
  public:
    TimeSeriesFit(const MatrixXd F_,
                  const MatrixXd G_,
                  const MatrixXd W_,
                  const double gamma_,
                  const int upsilon_, 
                  const MatrixXd Xi_,
                  const MatrixXd M0_,
                  const MatrixXd C0_,
                  const VectorXd observations_) :
    F(F_), G(G_), W(W_), gamma(gamma_), upsilon(upsilon_), Xi(Xi_), M0(M0_), C0(C0_), observations(observations_)
    {
      // note: eta is transposed in the DLM specification!
      filtered = false;
      smoothed = false;
      simulation_smoothed = false;
    }
    ~TimeSeriesFit(){}        
    
    // assuming time-invariant F, G, W, and gamma
    // need to adapt for time varying parameters (TODO)
    void apply_Kalman_filter(MatrixXd eta) {
      this->eta = eta;
      T = observations.maxCoeff(); // max observation assuming earliest observation is zero/baseline
      N = eta.rows();              // number of samples
      D = eta.cols()+1;            // number of multinomial categories minus 1
      system_dim = G.cols();
      Thetas_filtered = MatrixXd(system_dim*(D-1), T); // single marginal sample at each t
      Rs = MatrixXd(system_dim*system_dim, T);
      Ms = MatrixXd(system_dim*(D-1), T);
      Cs = MatrixXd(system_dim*system_dim, T);
      // init the objects we'll iteratively overwrite where they've got an initial value
      int upsilon_t = upsilon;
      MatrixXd Gt = G.transpose();
      MatrixXd Ft = F.transpose();
      MatrixXd Xi_t = Xi;
      MatrixXd M_t = M0;
      MatrixXd C_t = C0;
      // instantiate the others
      MatrixXd A_t(system_dim, D-1);
      MatrixXd R_t(system_dim, system_dim);
      MatrixXd ft_t(1, D-1);
      double q_t;
      MatrixXd et_t(1, D-1);
      MatrixXd S_t(system_dim, 1);
      MatrixXd LV(D-1, D-1);
      Eigen::LLT<Eigen::MatrixXd> lltOfCt;
      MatrixXd LU;
      MatrixXd Theta_t(system_dim, system_dim);
      VectorXd res(T);
      int curr_obs_idx = -1;
      bool found;
      for(int t=1; t<=T; t++) {
        // find the index of this time point in the observation vector (if it exists)
        // this will skip all but the first observation on a given date
        found = false;
        for(int k=curr_obs_idx+1; k<observations.size(); k++) {
          if(observations(k) == t) {
            curr_obs_idx = k;
            found = true;
          }
        }
        // system prior at t
        A_t = G*M_t;
        R_t = G*C_t*Gt + W;
        if(!found) {
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
        lltOfCt = Eigen::LLT<Eigen::MatrixXd>(C_t);
        LU = lltOfCt.matrixL();
        // sample Sigma(t); note: this returns Chol. upper triangular
        // swap this for thread_inplace version (TODO)
        LV = rInvWishRevCholesky(upsilon_t, Xi_t).matrix().transpose();
        // sample Theta(t)
        // swap this for thread_inplace version (TODO)
        Theta_t = rMatNormalCholesky(M_t, LU, LV);
        // pack up the samples to save
        for(int i=0; i<system_dim; i++) {
          for(int j=0; j<(D-1); j++) {
            Thetas_filtered(i+(j*system_dim),t-1) = Theta_t(i,j);
            Ms(i+(j*system_dim),t-1) = M_t(i,j);
          }
        }
        for(int i=0; i<system_dim; i++) {
          for(int j=0; j<system_dim; j++) {
            Rs(i+(j*system_dim),t-1) = R_t(i,j);
            Cs(i+(j*system_dim),t-1) = C_t(i,j);
          }
        }
      }
      upsilonT = upsilon_t;
      XiT = Xi_t;
      filtered = true;
    }

    // utility function
    MatrixXd unpack_sample(MatrixXd samples, int no_rows, int no_cols, int sample_idx) {
      // add check for bad dimensions here (TODO)
      MatrixXd sample(no_rows, no_cols);
      for(int i=0; i<no_rows; i++) {
        for(int j=0; j<no_cols; j++) {
          sample(i,j) = samples(i+(j*no_rows),sample_idx);
        }
      }
      return(sample);
    }

    // assuming time-invariant F, G, W, and gamma
    // need to adapt for time varying parameters (TODO)
    void apply_simulation_smoother(MatrixXd eta) {
      if(!filtered) {
        apply_Kalman_filter(eta);
      }
      Thetas_simulation_smoothed = MatrixXd(system_dim*(D-1), T); // single sample over the whole trajectory
      MatrixXd smoothed_Thetas(system_dim*(D-1), T);
      MatrixXd LV = rInvWishRevCholesky(upsilonT, XiT).matrix().transpose();
      // grab M_T
      MatrixXd M_t = unpack_sample(Ms, system_dim, (D-1), T-1);
      // grab C_T
      MatrixXd C_t = unpack_sample(Cs, system_dim, system_dim, T-1);
      Eigen::LLT<MatrixXd> lltOfCt(C_t);
      MatrixXd LU = lltOfCt.matrixL();
      // sample Theta_T
      MatrixXd smoothed_Theta_t = rMatNormalCholesky(M_t, LU, LV);
      for(int i=0; i<system_dim; i++) {
        for(int j=0; j<(D-1); j++) {
          smoothed_Thetas(i+(j*system_dim),T-1) = smoothed_Theta_t(i,j);
        }
      }
      MatrixXd R_t(system_dim, system_dim);
      MatrixXd R_t_inv(system_dim, system_dim);
      MatrixXd Z_t(system_dim, system_dim);
      MatrixXd A_t(system_dim, D-1);
      MatrixXd M_t_star(system_dim, D-1);
      MatrixXd C_t_star(system_dim, system_dim);
      for(int t=T-1; t>0; t--) {
        // note to self: 1-indexed loop to 0-indexed data structure
        R_t = unpack_sample(Rs, system_dim, system_dim, t); // R_t+1
        R_t_inv = R_t.inverse();
        M_t = unpack_sample(Ms, system_dim, D-1, t-1);
        C_t = unpack_sample(Cs, system_dim, system_dim, t-1);
        Z_t = C_t*(G.transpose())*R_t_inv;
        A_t = G*M_t;
        M_t_star = M_t + Z_t*(smoothed_Theta_t - A_t);
        C_t_star = C_t - Z_t*R_t*(Z_t.transpose());
        lltOfCt = Eigen::LLT<MatrixXd>(C_t_star); // reuse
        LU = lltOfCt.matrixL();
        smoothed_Theta_t = rMatNormalCholesky(M_t_star, LU, LV); // reuse
        for(int i=0; i<system_dim; i++) {
          for(int j=0; j<(D-1); j++) {
            Thetas_simulation_smoothed(i+(j*system_dim),t-1) = smoothed_Theta_t(i,j);
          }
        }
      }
      simulation_smoothed = true;
    }
};


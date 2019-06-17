# for reference, individuals passable as arguments are:
# "DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI"

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 6) {
  stop("Testing usage: Rscript rol_fit_labraduck.R ACA TRUE FALSE TRUE FALSE 100\nScaled up usage: Rscript rol_fit_labraduck.R ACA FALSE TRUE TRUE TRUE 200", call.=FALSE)
}

devtools::load_all("/data/mukherjeelab/labraduck")
source("rol_includes.R")

baboon <- args[1]
subset_time <- as.logical(args[2])
eval_MAP <- as.logical(args[3])
use_smooth <- as.logical(args[4])
save_fit <- as.logical(args[5])
n_samples <- as.numeric(args[6])

cat("Fitting",baboon,"over all taxa...\n")
load(paste0(data_path,"/",baboon,"_data.RData"))
  
W <- matrix(0, 3, 3)
F <- matrix(c(1, 0, 1), 3, 1)
  
diag(W) <- c(1, 1, 1/100)
  
# 1:1 signal:noise should be given by something like
# Tr(gamma_t * Sigma) = Tr(W_t * Sigma) = gamma_t = Tr(W_t)
# i.e. fixed_gamma_scale <- sum(diag(W)*fixed_W_scale)
# but this looks super fucked up in practice
fixed_W_scale <- 0
fixed_gamma_scale <- 0
cat("Using W_scale:",fixed_W_scale,"\n")
cat("Using gamma_scale:",fixed_gamma_scale,"\n\n")
  
# ALR prior covariance
upsilon <- D-1+10 # lesser certainty
# supsilon <- D-1+20 # greater certainty; this should tighten the distribution around this mean
GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference;
# take diag as covariance over log abundances
Xi <- GG%*%(diag(D)*1)%*%t(GG)
# mean-center
Xi <- Xi*(upsilon-D-1)

fit_obj <- fit_model(indiv_data, W, F, gamma_scale=fixed_gamma_scale, W_scale=fixed_W_scale,
                     upsilon, Xi,
                     n_samples=n_samples, ret_mean=FALSE,
                     apply_smoother=use_smooth, subset_time=subset_time)

if(save_fit) {
  cat("Saving fit object...\n")
  Sigma_samples <- fit_obj$fit$Sigma
  save(Sigma_samples, file=paste0(data_path,"/",baboon,"_fit.RData"))
}

fit <- fit_obj$fit
Y <- fit_obj$Y
observations <- fit_obj$observations
  
# add back later -- plot covariance associated with simulation smoother samples
# cov_Theta <- plot_cov_Theta(fit$D, fit$T, n_samples, F, fit$Thetas_smoothed, baboon, save_path=save_path, as_corr=FALSE)
# png(paste0(save_path,"/",baboon,"_covTheta_corr.png"))
# image(cov2cor(cov_Theta))
# dev.off()

# high abundance
plot_posterior(Y, fit, F, observations, baboon, plot_what=c("smoothed", "dlm_eta"),
               save_path=save_path, lr_idx=1)
plot_posterior(Y, fit, F, observations, baboon, plot_what=c("smoothed", "dlm_eta"),
                 save_path=save_path, lr_idx=2)
# low abundance
plot_posterior(Y, fit, F, observations, baboon, plot_what=c("smoothed", "dlm_eta"),
               save_path=save_path, lr_idx=21)
plot_posterior(Y, fit, F, observations, baboon, plot_what=c("smoothed", "dlm_eta"),
               save_path=save_path, lr_idx=25)

if(eval_MAP) {
  fit_obj <- fit_model(indiv_data, W, F, gamma_scale=fixed_gamma_scale, W_scale=fixed_W_scale,
                       upsilon, Xi,
                       n_samples=0, ret_mean=TRUE,
                       apply_smoother=FALSE, subset_time=subset_time)
      
  fit <- fit_obj$fit
  Y <- fit_obj$Y
  plot_Sigma(fit, Y, baboon, save_path=save_path, as_corr=FALSE)
  plot_Sigma(fit, Y, baboon, save_path=save_path, as_corr=TRUE)
}

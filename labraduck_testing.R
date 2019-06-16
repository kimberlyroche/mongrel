library(stray)
library(phyloseq)
library(ggplot2)
library(Rcpp)
library(RcppEigen)

#rm(list=ls())

sourceCpp("src/cov_viz_test.cpp")

plot_Sigma <- function(fit, Y, baboon, save_path="", append="", as_corr=FALSE) {
  if(as_corr) {
    png(paste0(save_path,"/",baboon,"_Sigma_corr",append,".png"))
    image(cov2cor(fit$Sigma[,,1]))
  } else {
    png(paste0(save_path,"/",baboon,"_Sigma",append,".png"))
    image(fit$Sigma[,,1])
  }
  dev.off()
  if(as_corr) {
    png(paste0(save_path,"/",baboon,"_empcov_corr",append,".png"))
    image(cov2cor(cov(driver::alr(t(Y + 0.5)))))
  } else {
    png(paste0(save_path,"/",baboon,"_empcov",append,".png"))
    image(cov(driver::alr(t(Y + 0.5))))
  }
  dev.off()
}

plot_cov_Theta <- function(D, T, n_samples, F, Thetas, baboon="", save_path="", as_corr=FALSE) {
  # dimensions of the smoothed sample Thetas are (Q x D-1 x T) x n_samples
  # n_samples as passed can be a subset of the samples
  Thetas <- Thetas[,1:n_samples]
  flattened_Theta <- matrix(NA, T*n_samples, D-1)
  for(i in 1:n_samples) {
    #cat(i,"\n")
    temp <- Thetas[,i]
    dim(temp) <- c(nrow(F), D-1, T)
    for(t in 1:T) {
      temp2 <- temp[,,t] # 3 x D-1 (etc.)
      flattened_Theta[T*(i-1) + t,] <- t(F)%*%temp2
    }
  }
  cov_obj <- t(flattened_Theta)%*%flattened_Theta/(T*n_samples)
  if(as_corr) {
    png(paste0(save_path,"/",baboon,"_covTheta_corr.png"))
    image(cov2cor(cov_obj))
  } else {
    png(paste0(save_path,"/",baboon,"_covTheta.png"))
    image(cov_obj)
  }
  dev.off()
  return(cov_obj)
}

get_Sigma_distances <- function(fit, method="fro", percent_sample=0.1, as_corr=FALSE) {
  # we'll use 10% of the samples
  n_samples <- dim(fit$Sigma)[3]
  D <- dim(fit$Sigma)[1]+1
  sample_idx <- sample(n_samples)[1:round(n_samples*percent_sample)]
  upper <- length(sample_idx)
  dist_distro <- numeric((upper^2 - upper)/2)
  all_sums <- numeric(upper*2) # calculate sum of elements of A, B
  all_sd <- numeric(upper*2) # calculate sd of elements of A, B
  idx <- 1
  for(i in 1:(upper-1)) {
    for(j in (i+1):upper) {
      A <- fit$Sigma[,,sample_idx[i]]
      B <- fit$Sigma[,,sample_idx[j]]
      if(as_corr) {
        A <- cov2cor(A)
        B <- cov2cor(B)
      }
      all_sums[((idx-1)*2)+1] <- sum(c(A))
      all_sums[((idx-1)*2)+2] <- sum(c(B))
      all_sd[((idx-1)*2)+1] <- sd(c(A))
      all_sd[((idx-1)*2)+2] <- sd(c(B))
      if(method == "fro") {
        d <- fro_dist(A, B)
      } else {
        d <- riemann_dist(A, B)
      }
      dist_distro[idx] <- d
      idx <- idx + 1
    }
  }
  return(list(distro=dist_distro, avg_sum=mean(all_sums), avg_dev=mean(all_sd)/(mean(all_sums)*(D-1)^2)))
}

plot_distance_distros <- function(fit_models, save_path="") {
  # plot the Sigma distance distributions for the models
  df <- data.frame(distance=c(), method=c(), W_scale=c())
  df_corr <- data.frame(distance=c(), method=c(), W_scale=c())
  for(i in 1:length(fit_models)) {
    cat("Calculating distances for model no.",i,"...\n")
    fit <- fit_models[[i]]
    temp <- get_Sigma_distances(fit, method="fro", percent_sample=0.1, as_corr=FALSE)
    df_sub <- data.frame(distance=temp$distro, method="Frobenius", W_scale=W_scales[i])
    df <- rbind(df, df_sub)
    temp <- get_Sigma_distances(fit, method="fro", percent_sample=0.1, as_corr=TRUE)
    df_sub <- data.frame(distance=temp$distro, method="Frobenius", W_scale=W_scales[i])
    df_corr <- rbind(df_corr, df_sub)
  }
  df$W_scale <- as.factor(df$W_scale)
  df_corr$W_scale <- as.factor(df$W_scale)
  # if we add in Riemannian distance might want: facet_wrap(~method, scales="free_x")
  p <- ggplot(data=df, aes(x=distance)) + geom_density(aes(fill=W_scale), alpha = 0.4) +
    theme_minimal()
  ggsave(filename=paste0(save_path,"/Frobenius_Sigma_distances_cov.png"),
         units="in", scale=2, width=4, height=3)
  p <- ggplot(data=df_corr, aes(x=distance)) + geom_density(aes(fill=W_scale), alpha = 0.4) +
    theme_minimal()
  ggsave(filename=paste0(save_path,"/Frobenius_Sigma_distances_corr.png"),
         units="in", scale=2, width=4, height=3)
}

plot_distant_samples <- function(fit_models, save_path="") {
  # cherry pick a couple of Sigma samples that are very difference and very similar and visually compare them!
  fit <- fit_models[[length(fit_models)]]
  n_samples <- dim(fit$Sigma)[3]
  i <- 1
  vdiff_idx_corr <- -1
  vsim_idx_corr <- -1
  vdiff_idx_cov <- -1
  vsim_idx_cov <- -1
  min_dist_corr <- Inf
  max_dist_corr <- -Inf
  min_dist_cov <- Inf
  max_dist_cov <- -Inf
  for(j in 2:n_samples) {
    A <- cov2cor(fit$Sigma[,,i])
    B <- cov2cor(fit$Sigma[,,j])
    d <- fro_dist(A, B)
    if(d < min_dist_corr) {
      min_dist_corr <- d
      vsim_idx_corr <- j
    }
    if(d > max_dist_corr) {
      max_dist_corr <- d
      vdiff_idx_corr <- j
    }
    A <- fit$Sigma[,,i]
    B <- fit$Sigma[,,j]
    d <- fro_dist(A, B)
    if(d < 6) {
      cat("HERE\n")
    }
    if(d < min_dist_cov) {
      min_dist_cov <- d
      vsim_idx_cov <- j
    }
    if(d > max_dist_cov) {
      max_dist_cov <- d
      vdiff_idx_cov <- j
    }
  }
  
  A <- cov2cor(fit$Sigma[,,1])
  B <- cov2cor(fit$Sigma[,,vdiff_idx_corr])
  d <- fro_dist(A, B)
  cat("Max distance (correlation) is:",d,"\n")
  png(paste0(save_path,"/diffA_corr.png"))
  image(A)
  dev.off()
  png(paste0(save_path,"/diffB_corr.png"))
  image(B)
  dev.off()
  B <- cov2cor(fit$Sigma[,,vsim_idx_corr])
  d <- fro_dist(A, B)
  cat("Min distance (correlation) is:",d,"\n")
  png(paste0(save_path,"/simB_corr.png"))
  image(B)
  dev.off()
  
  A <- fit$Sigma[,,1]
  B <- fit$Sigma[,,vdiff_idx_cov]
  d <- fro_dist(A, B)
  cat("Max distance (covariance) is:",d,"\n")
  png(paste0(save_path,"/diffA_cov.png"))
  image(A)
  dev.off()
  png(paste0(save_path,"/diffB_cov.png"))
  image(B)
  dev.off()
  B <- fit$Sigma[,,vsim_idx_cov]
  d <- fro_dist(A, B)
  cat("Min distance (covariance) is:",d,"\n")
  png(paste0(save_path,"/simB_cov.png"))
  image(B)
  dev.off()
  
  j <- vdiff_idx_cov # very similar
  B <- fit$Sigma[,,j]
  
  df <- data.frame(value=c(A), matrix="reference sample")
  df <- rbind(df, data.frame(value=c(B), matrix="different sample"))
  
  j <- vsim_idx_cov # very similar
  B <- fit$Sigma[,,j]
  
  df <- rbind(df, data.frame(value=c(B), matrix="similar sample"))
  
  p <- ggplot(data=df, aes(x=value)) + geom_density(aes(fill=matrix), alpha = 0.4) +
    theme_minimal()
  ggsave(filename=paste0(save_path,"/element_distributions_covariance.png"),
         units="in", scale=2, width=4, height=3)
}

plot_element_distros <- function(fit_models, W_scales, percent_sample=0.1, save_path="") {
  df <- data.frame(value=c(), W_scale=c())
  n_samples <- dim(fit$Sigma)[3]
  sample_idx <- sample(n_samples)[1:round(n_samples*percent_sample)]
  for(i in 1:length(fit_models)) {
    fit <- fit_models[[i]]
    # there are a shitton of these, just sample
    elements <- c(fit$Sigma[,,sample_idx])
    df <- rbind(df, data.frame(value=elements, W_scale=W_scales[i]))
  }
  df$W_scale <- as.factor(df$W_scale)
  p <- ggplot(data=df, aes(x=value)) + geom_density(aes(fill=W_scale), alpha = 0.4) +
    theme_minimal()
  ggsave(filename=paste0(save_path,"/element_distributions_W_sweep.png"),
         units="in", scale=2, width=4, height=3)
}

# assumes both smoother and filter have been run!
plot_posterior <- function(Y, fit, F, observations, baboon, save_path="", 
                           plot_what=c("filtered", "smoothed_mean", "smoothed", "dlm_eta", "optimized_eta"),
                           append_str="", lr_idx=1, ylim=NULL) {
  # plot the first sample of (1) log-transformed Y (2) approximated eta (3) filtered Theta (with F applied)
  # these should be similar
  T <- max(observations)
  N <- length(observations)
  etaopt_ymin <- rep(Inf, N)
  etaopt_ymax <- rep(-Inf, N)
  filtered_ymin <- rep(Inf, T)
  filtered_ymax <- rep(-Inf, T)
  smoothed_ymin_Theta <- rep(Inf, T)
  smoothed_ymax_Theta <- rep(-Inf, T)
  smoothed_ymin_MStar <- rep(Inf, T)
  smoothed_ymax_MStar <- rep(-Inf, T)
  smoothed_ymin_eta <- rep(Inf, T)
  smoothed_ymax_eta <- rep(-Inf, T)
  Ft <- t(F)
  Q <- length(F)
  lr <- driver::alr(t(Y + 0.5))
  n_samples <- dim(fit$Eta)[3]
  D <- dim(fit$Eta)[1]+1
  for(s in 1:n_samples) {
    if("filtered" %in% plot_what) {
      ThetasF_1T <- fit$Thetas_filtered[,s]
      dim(ThetasF_1T) <- c(Q, D-1, T)
    }
    if("smoothed_mean" %in% plot_what) {
      MStar_1T <- fit$M_star[,s]
      dim(MStar_1T) <- c(Q, D-1, T)
    }
    if("smoothed" %in% plot_what) {
      ThetasS_1T <- fit$Thetas_smoothed[,s]
      dim(ThetasS_1T) <- c(Q, D-1, T)
    }
    if("dlm_eta" %in% plot_what) {
      EtaDLM_1T <- fit$Eta_DLM[,s]
      dim(EtaDLM_1T) <- c(D-1, T)
    }
    if("optimized_eta" %in% plot_what) {
      Eta_1T <- fit$Eta[,,s]
      for(t in 1:N) {
        Etaopt_lr <- Eta_1T[lr_idx,t]
        if(Etaopt_lr < etaopt_ymin[t]) {
          etaopt_ymin[t] <- Etaopt_lr
        }
        if(Etaopt_lr > etaopt_ymax[t]) {
          etaopt_ymax[t] <- Etaopt_lr
        }
      }
    }
    for(t in 1:T) {
      if("filtered" %in% plot_what) {
        Theta_f_t <- ThetasF_1T[,,t]
        filtered_pt <- (Ft%*%Theta_f_t)[lr_idx]
        if(filtered_pt < filtered_ymin[t]) {
          filtered_ymin[t] <- filtered_pt
        }
        if(filtered_pt > filtered_ymax[t]) {
          filtered_ymax[t] <- filtered_pt
        }
      }
      if("smoothed_mean" %in% plot_what) {
        MStar_t <- MStar_1T[,,t]
        smoothed_mean <- (Ft%*%MStar_t)[lr_idx]
        if(smoothed_mean < smoothed_ymin_MStar[t]) {
          smoothed_ymin_MStar[t] <- smoothed_mean
        }
        if(smoothed_mean > smoothed_ymax_MStar[t]) {
          smoothed_ymax_MStar[t] <- smoothed_mean
        }
      }
      if("smoothed" %in% plot_what) {
        Theta_s_t <- ThetasS_1T[,,t]
        smoothed_pt <- (Ft%*%Theta_s_t)[lr_idx]
        if(smoothed_pt < smoothed_ymin_Theta[t]) {
          smoothed_ymin_Theta[t] <- smoothed_pt
        }
        if(smoothed_pt > smoothed_ymax_Theta[t]) {
          smoothed_ymax_Theta[t] <- smoothed_pt
        }
      }
      if("dlm_eta" %in% plot_what) {
        Eta_lr <- EtaDLM_1T[lr_idx,t]
        if(Eta_lr < smoothed_ymin_eta[t]) {
          smoothed_ymin_eta[t] <- Eta_lr
        }
        if(Eta_lr > smoothed_ymax_eta[t]) {
          smoothed_ymax_eta[t] <- Eta_lr
        }
      }
    }
  }
  df <- data.frame(timepoint=as(observations, "vector"), alrY=lr[,lr_idx])
  if("filtered" %in% plot_what) {
    df2 <- data.frame(timepoint=1:T, ymin_f=filtered_ymin[1:T], ymax_f=filtered_ymax[1:T])
    df <- merge(df, df2, by='timepoint', all=TRUE)
  }
  if("smoothed_mean" %in% plot_what) {
    df2 <- data.frame(timepoint=1:T, ymin_ms=smoothed_ymin_MStar[1:T], ymax_ms=smoothed_ymax_MStar[1:T])
    df <- merge(df, df2, by='timepoint', all=TRUE)
  }
  if("smoothed" %in% plot_what) {
    df2 <- data.frame(timepoint=1:T, ymin_s=smoothed_ymin_Theta[1:T], ymax_s=smoothed_ymax_Theta[1:T])
    df <- merge(df, df2, by='timepoint', all=TRUE)
  }
  if("dlm_eta" %in% plot_what) {
    df2 <- data.frame(timepoint=1:T, ymin_etadlm=smoothed_ymin_eta[1:T], ymax_etadlm=smoothed_ymax_eta[1:T])
    df <- merge(df, df2, by='timepoint', all=TRUE)
  }
  if("optimized_eta" %in% plot_what) {
    df2 <- data.frame(timepoint=as(observations, "vector"), ymin_etaopt=etaopt_ymin, ymax_etaopt=etaopt_ymax)
    df <- merge(df, df2, by='timepoint', all=TRUE)
  }
  p <- ggplot(df, aes(timepoint))
  gray_idx <- 1
  grays <- c("grey90", "grey80", "grey70", "grey60")
  if("filtered" %in% plot_what) {
    p <- p + geom_ribbon(aes(ymin=ymin_f, ymax=ymax_f), fill=grays[gray_idx])
    gray_idx <- gray_idx + 1
  }
  if("dlm_eta" %in% plot_what) {
    p <- p + geom_ribbon(aes(ymin=ymin_etadlm, ymax=ymax_etadlm), fill=grays[gray_idx])
    gray_idx <- gray_idx + 1
  }
  if("smoothed" %in% plot_what) {
    p <- p + geom_ribbon(aes(ymin=ymin_s, ymax=ymax_s), fill=grays[gray_idx])
    gray_idx <- gray_idx + 1
  }
  if("smoothed_mean" %in% plot_what) {
    p <- p + geom_ribbon(aes(ymin=ymin_ms, ymax=ymax_ms), fill=grays[gray_idx])
    gray_idx <- gray_idx + 1
  }
  if("optimized_eta" %in% plot_what) {
    p <- p + geom_point(aes(x=timepoint, y=ymin_etaopt), color="red", size=1) + geom_point(aes(x=timepoint, y=ymax_etaopt), color="red", size=1)
  }
  p <- p + geom_point(aes(x=timepoint, y=alrY), color="blue", size=1) + theme_minimal()
  if(lr_idx == 1) {
    p <- p + ylab("ALR(Bifidobacteriaceae/Helicobacteraceae)")
  } else if(lr_idx == 2) {
    p <- p + ylab("ALR(Prevotellaceae/Helicobacteraceae)")
  } else if(lr_idx == 7) {
    p <- p + ylab("ALR(Muribaculaceae/Helicobacteraceae)")
  } else if(lr_idx == 21) {
    p <- p + ylab("ALR(Christensenellaceae/Helicobacteraceae)")
  } else {
    # index 25
    p <- p + ylab("ALR(Lactobacillaceae/Helicobacteraceae)")
  }
  #show(p)
  if(!is.null(ylim)) {
    p <- p + ylim(ylim[1], ylim[2])
  }
  width <- round(log(T/100)*3)
  ggsave(filename=paste0(save_path,"/",baboon,"_Theta_smoothed",append_str,"_LR",lr_idx,".png"),
         units="in", scale=2, width=width, height=2)
}

record_time <- function(fit, baboon, save_path="", append_str="") {
  sink(paste0(save_path,"/",baboon,append_str,".txt"), append=FALSE, split=FALSE)
  for(i in 1:length(fit$Timer)) {
    cat(names(fit$Timer)[i]," -- ",round(fit$Timer[i],1),"\n")
  }
  sink()
}

fit_model <- function(indiv_data, W, F, gamma_scale=0, W_scale=0, upsilon, Xi,
                      n_samples=2000, ret_mean=FALSE,
                      apply_smoother=FALSE, subset_time=TRUE, useSylv=FALSE) {
  Y_full <- indiv_data$ys
  alr_Y_full <- driver::alr(Y_full + 0.5)
  alr_means <- colMeans(alr_Y_full)
  observations_full <- indiv_data$observation_vec
  
  if(subset_time) {
    date_lower_limit <- "2001-10-01"
    date_upper_limit <- "2002-11-30"
    date_upper_limit <- "2003-11-30"

    min_idx <- min(which(names(observations_full) >= date_lower_limit))
    max_idx <- max(which(names(observations_full) <= date_upper_limit))

    # subset dates
    Y <- t(Y_full[min_idx:max_idx,])
    colnames(Y) <- NULL
    rownames(Y) <- NULL
    # observations currently requires baseline observation at t=1
    observations <- matrix(observations_full[min_idx:max_idx], nrow=1) - observations_full[min_idx] + 1
  } else {
    Y <- t(Y_full)
    colnames(Y) <- NULL
    rownames(Y) <- NULL
    observations <- matrix(observations_full, nrow=1)
  }

  D <- nrow(Y)
  N <- length(observations)
  if(N < (D-1)) {
    cat("Too few observations (",N,")...\n")
    #return(NULL)
  }
  omega <- 2*pi/365
  G <- matrix(c(cos(omega), -sin(omega), 0, sin(omega), cos(omega), 0, 0, 0, 1), 3, 3)
  C0 <- W
  C0[3,3] <- C0[3,3]*10

  # set Fourier coefficients uniformly at 1; set offset to mean for this logratio over all timepoints
  M0 <- matrix(1, 3, D-1)
  #M0[3,] <- alr_means

  fit <- labraduck(Y=Y, upsilon=upsilon, Xi=Xi,
                   gamma=gamma, F=F, G=G, W=W, M0=M0, C0=C0, observations=observations,
                   gamma_scale=gamma_scale, W_scale=W_scale,
                   max_iter=100000, b1=0.95, step_size=0.002, eps_f=1e-11, decomp_method="eigen",
                   n_samples=n_samples, ret_mean=ret_mean,
                   apply_smoother=apply_smoother, useSylv=useSylv, verbose=FALSE)
  return(list(fit=fit, Y=Y, observations=observations))
}

best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")
best_sampled <- c("ACA", "DUI", "LOG", "VET")

subset_time <- TRUE
eval_MAP <- TRUE

n_samples <- 100

save_path <- "C:/Users/kim/Documents/rules_of_life/plots/stray"
#save_path <- "/Users/ladlab/Desktop/temp"

data_path <- "C:/Users/kim/Documents/rules_of_life/subsetted_indiv_data"
#data_path <- "/Users/ladlab/Desktop/indiv_baboons"

# very crude flag to me to indicate subsequent runs converged to different places
W_scales_fit <- list()
flagged <- c()

for(baboon in best_sampled) {
  cat("Fitting",baboon,"over all taxa...\n")
  load(paste0(data_path,"/",baboon,"_data.RData"))
  
  D <- ncol(indiv_data$ys)
  W <- matrix(0, 3, 3)
  F <- matrix(c(1, 0, 1), 3, 1)
  
  fixed_gamma_scale <- 0
  fixed_W_scale <- 0
  
  diag(W) <- c(1, 1, 1/100)
  
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
                       apply_smoother=TRUE, subset_time=subset_time)
  
  W_scales_fit[[baboon]] <- fit_obj$fit$W_scale
  
  cat("Saving fit object...\n")
  Sigma_samples <- fit_obj$fit$Sigma
  save(Sigma_samples, file=paste0(data_path,"/",baboon,"_fit.RData"))

  if(!is.null(fit_obj)) {
    fit <- fit_obj$fit
    Y <- fit_obj$Y
    observations <- fit_obj$observations
  
    # plot covariance associated with simulation smoother samples
    #cov_Theta <- plot_cov_Theta(fit$D, fit$T, n_samples, F, fit$Thetas_smoothed, baboon, save_path=save_path, as_corr=FALSE)
    #png(paste0(save_path,"/",baboon,"_covTheta_corr.png"))
    #image(cov2cor(cov_Theta))
    #dev.off()

    # high abundance
    plot_posterior(Y, fit, F, observations, baboon, plot_what=c("smoothed", "dlm_eta"),
                   save_path=save_path, lr_idx=1)
    #plot_posterior(Y, fit, F, observations, baboon, plot_what=c("smoothed", "dlm_eta"),
    #                 save_path=save_path, lr_idx=2)
    # low abundance
    #plot_posterior(Y, fit, F, observations, baboon, plot_what=c("smoothed", "dlm_eta"),
    #               save_path=save_path, lr_idx=21)
    plot_posterior(Y, fit, F, observations, baboon, plot_what=c("smoothed", "dlm_eta"),
                   save_path=save_path, lr_idx=25)

    if(eval_MAP) {
      fit_obj <- fit_model(indiv_data, W, F, gamma_scale=fixed_gamma_scale, W_scale=fixed_W_scale,
                           upsilon, Xi,
                           n_samples=0, ret_mean=TRUE,
                           apply_smoother=FALSE, subset_time=subset_time)
      
      if(abs(fit_obj$fit$W_scale - W_scales_fit[[baboon]]) > 0.1) {
        flagged <- c(flagged, baboon)
      }
      
      fit <- fit_obj$fit
      Y <- fit_obj$Y
      plot_Sigma(fit, Y, baboon, save_path=save_path, as_corr=FALSE)
      plot_Sigma(fit, Y, baboon, save_path=save_path, as_corr=TRUE)
    }
  }
}

if(length(flagged) > 0) {
  cat("Local minima flagged for:",flagged,"\n")
}

n_indiv <- length(best_sampled)
all_samples <- matrix(NA, D-1, (D-1)*n_samples*n_indiv)
labels <- c()
for(i in 1:n_indiv) {
  load(paste0(data_path,"/",best_sampled[i],"_fit.RData"))
  dim(Sigma_samples) <- c((D-1), (D-1)*n_samples)
  all_samples[,((i-1)*(D-1)*n_samples+1):(i*(D-1)*n_samples)] <- Sigma_samples
  labels <- c(labels, rep(best_sampled[i], n_samples))
}

distance_mat <- matrix(NA, n_samples*n_indiv, n_samples*n_indiv)
for(i in 1:(n_indiv*n_samples)) {
  for(j in i:(n_indiv*n_samples)) {
    i_idx <- (i-1)*(D-1)
    A <- all_samples[,(i_idx+1):(i_idx+(D-1))]
    j_idx <- (j-1)*(D-1)
    B <- all_samples[,(j_idx+1):(j_idx+(D-1))]
    distance_mat[i,j] <- mat_dist(A, B, use_Riemann=TRUE)
    distance_mat[j,i] <- distance_mat[i,j]
  }
}

fit <- cmdscale(distance_mat, eig=TRUE, k=2) # k is the number of dim
cat("Lambda 1:",fit$eig[1],"\n")
cat("Lambda 2:",fit$eig[2],"\n")
cat("Lambda 3:",fit$eig[3],"\n")

df <- data.frame(x=fit$points[,1], y=fit$points[,2], labels=labels)
p <- ggplot(df, aes(x=x, y=y, color=labels)) +
  geom_point() +
  ggtitle("Riemannian distance (4 indiv)")
p
ggsave(paste0(save_path,"/posterior_ordination_test.png"), scale=2,
       width=4, height=4, units="in", dpi=100)

# TODO: how to visualize distance: find most dissimiliar pair within and most dissimilar pair between
# max_dist <- which(distance_mat == max(distance_mat), arr.ind = TRUE)
# max_id1 <- max_dist[1,"row"][[1]]
# max_id2 <- max_dist[1,"col"][[1]]

# first 2 samples for individual 1
image(all_samples[,1:(D-1)])
image(all_samples[,((D-1)+1):(2*(D-1))])

image(all_samples[,1:(D-1)])
image(all_samples[,(((D-1)*n_samples)+1):(((D-1)*n_samples)+(D-1))])










library(stray)
library(driver)
library(phyloseq)
library(ggplot2)
library(Rcpp)
library(RcppEigen)

sourceCpp("src/cov_viz_test.cpp")

best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")

save_path <- "/data/mukherjeelab/rulesoflife/plots/stray"
#save_path <- "C:/Users/kim/Documents/rules_of_life/plots/stray"
#save_path <- "/Users/ladlab/Desktop/temp"

data_path <- "/data/mukherjeelab/rulesoflife/subsetted_indiv_data"
#data_path <- "C:/Users/kim/Documents/rules_of_life/subsetted_indiv_data"
#data_path <- "/Users/ladlab/Desktop/indiv_baboons"

D <- 26

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
         units="in", scale=2, dpi=100, width=width, height=2)
}

record_time <- function(fit, baboon, save_path="", append_str="") {
  sink(paste0(save_path,"/",baboon,append_str,".txt"), append=FALSE, split=FALSE)
  for(i in 1:length(fit$Timer)) {
    cat(names(fit$Timer)[i]," -- ",round(fit$Timer[i],1),"\n")
  }
  sink()
}

fit_model <- function(indiv_data, W, F, gamma_scale=0, W_scale=0, upsilon, Xi, M0, C0,
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

  fit <- labraduck(Y=Y, upsilon=upsilon, Xi=Xi,
                   gamma=gamma, F=F, G=G, W=W, M0=M0, C0=C0, observations=observations,
                   gamma_scale=gamma_scale, W_scale=W_scale,
                   max_iter=250000, b1=0.97, step_size=0.002, eps_f=1e-12, decomp_method="eigen",
                   n_samples=n_samples, ret_mean=ret_mean,
                   apply_smoother=apply_smoother, useSylv=useSylv, verbose=FALSE)
  return(list(fit=fit, Y=Y, observations=observations))
}

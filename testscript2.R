library(stray)
library(phyloseq)
library(ggplot2)

rm(list=ls())

fro_dist <- function(A, B) {
  D <- A - B
  sum_all <- sum(c(D^2))
  return(sqrt(sum_all))
}

# TODO: this needs to be debugged; it gives a different ordering over distance distributions
# than the (dead simple) Frobenius norm distance OR is this scale-invariant?
riemann_dist <- function(A, B) {
  Ainv <- solve(A)
  U <- chol(A) # upper triangular s.t. A <- t(U)%*%U
  X <- t(U)%*%B%*%U
  eigX <- eigen(X)
  # if V are eigenvectors of A and A' is a diagonal matrix containing eigenvalues of A
  # then log(A) = V%*%log(A')%*%V^(-1)
  logX <- eigX$vectors%*%diag(log(eigX$values))%*%solve(eigX$vectors)
  return(sum(diag(logX^2)))
}

plot_Sigma <- function(fit, Y, baboon, save_path="", append="", as_corr=FALSE) {
  png(paste0(save_path,"/",baboon,"_Sigma",append,".png"))
  if(as_corr) {
    image(cov2cor(fit$Sigma[,,1]))
  } else {
    image(fit$Sigma[,,1])
  }
  dev.off()
  png(paste0(save_path,"/",baboon,"_empcov",append,".png"))
  if(as_corr) {
    image(cov2cor(cov(driver::alr(t(Y + 0.5)))))
  } else {
    image(cov(driver::alr(t(Y + 0.5))))
  }
  dev.off()
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
plot_posterior <- function(Y, fit, F, observations, baboon, show_filtered=FALSE, save_path="", append_str="", lr_idx=1) {
  # plot the first sample of (1) log-transformed Y (2) approximated eta (3) filtered Theta (with F applied)
  # these should be similar
  show_Eta <- FALSE # can add back at will
  T <- max(observations)
  N <- length(observations)
  filtered_ymin <- rep(Inf, T)
  filtered_ymax <- rep(-Inf, T)
  smoothed_ymin <- rep(Inf, T)
  smoothed_ymax <- rep(-Inf, T)
  Ft <- t(F)
  Q <- length(F)
  lr <- driver::alr(t(Y + 0.5))
  n_samples <- dim(fit$Eta)[3]
  D <- dim(fit$Eta)[1]+1
  for(s in 1:n_samples) {
    ThetasF_1T <- fit$Thetas_filtered[,s]
    dim(ThetasF_1T) <- c(Q, D-1, T)
    ThetasS_1T <- fit$Thetas_smoothed[,s]
    dim(ThetasS_1T) <- c(Q, D-1, T)
    for(t in 1:(T-1)) {
      Theta_f_t <- ThetasF_1T[,,t]
      filtered_pt <- (Ft%*%Theta_f_t)[lr_idx]
      if(filtered_pt < filtered_ymin[t]) {
        filtered_ymin[t] <- filtered_pt
      }
      if(filtered_pt > filtered_ymax[t]) {
        filtered_ymax[t] <- filtered_pt
      }
      Theta_s_t <- ThetasS_1T[,,t]
      smoothed_pt <- (Ft%*%Theta_s_t)[lr_idx]
      if(smoothed_pt < smoothed_ymin[t]) {
        smoothed_ymin[t] <- smoothed_pt
        if(t == T) {
          cat("Smoothed min at s =",s,":",smoothed_pt,"\n")
        }
      }
      if(smoothed_pt > smoothed_ymax[t]) {
        smoothed_ymax[t] <- smoothed_pt
      }
    }
  }
  # don't show eta for now since there are multiple samples for every time point and logratio for eta
  # how best to represent?
  df <- data.frame(timepoint=as(observations, "vector"), alrY=lr[,lr_idx])
  if(show_Eta) {
    df2 <- data.frame(timepoint=as(observations, "vector"), eta_hat=fit$Eta[lr_idx,,n_samples]) # just use the last sample
    df <- merge(df, df2, by='timepoint', all=TRUE)
  }
  if(show_filtered) {
    df2 <- data.frame(timepoint=1:(T-1), ymin_f=filtered_ymin[1:(T-1)], ymax_f=filtered_ymax[1:(T-1)])
    df <- merge(df, df2, by='timepoint', all=TRUE)
  }
  df2 <- data.frame(timepoint=1:(T-1), ymin_s=smoothed_ymin[1:(T-1)], ymax_s=smoothed_ymax[1:(T-1)])
  df <- merge(df, df2, all=TRUE)
  p <- ggplot(df, aes(timepoint))
  if(show_filtered) {
    p <- p + geom_ribbon(aes(ymin=ymin_f, ymax=ymax_f), fill = "grey90")
  }
  p <- p + geom_ribbon(aes(ymin=ymin_s, ymax=ymax_s), fill = "grey80") +
    geom_point(aes(x=timepoint, y=alrY), color="blue")
  if(show_Eta) {
    p <- p + geom_point(aes(x=timepoint, y=eta_hat), color="red")
  }
  p <- p + theme_minimal() +
    ylab("ALR(Bifidobacteriaceae/Helicobacteraceae)")
  width <- round(log(T/100)*3)
  ggsave(filename=paste0(save_path,"/",baboon,"_Theta_smoothed",append_str,".png"),
         units="in", scale=2, width=width, height=2)
}

record_time <- function(fit, baboon, save_path="", append_str="") {
  sink(paste0(save_path,"/",baboon,append_str,".txt"), append=FALSE, split=FALSE)
  for(i in 1:length(fit$Timer)) {
    cat(names(fit$Timer)[i]," -- ",round(fit$Timer[i],1),"\n")
  }
  sink()
}

fit_model <- function(indiv_data, W, W_scale_init, F, n_samples=2000, ret_mean=FALSE, apply_smoother=FALSE, subset_time=TRUE) {
  Y_full <- indiv_data$ys
  alr_Y_full <- driver::alr(Y_full + 0.5)
  alr_means <- colMeans(alr_Y_full)
  observations_full <- indiv_data$observation_vec
  
  if(subset_time) {
    date_lower_limit <- "2001-10-01"
    date_upper_limit <- "2002-11-30"

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

  omega <- 2*pi/365
  G <- matrix(c(cos(omega), -sin(omega), 0, sin(omega), cos(omega), 0, 0, 0, 1), 3, 3)
  C0 <- W*10
  gamma <- 1

  D <- nrow(Y)
  # set Fourier coefficients uniformly at 1; set offset to mean for this logratio over all timepoints
  M0 <- matrix(1, 3, D-1)
  M0[3,] <- alr_means
  upsilon <- D+10
  Xi <- diag(D-1)
  fit <- labraduck(Y=Y, upsilon=upsilon, Xi=Xi,
                   gamma=gamma, F=F, G=G, W=W, W_scale_init=log(W_scale_init), M0=M0, C0=C0, observations=observations,
                   max_iter=100000, decomp_method="eigen", n_samples=n_samples, ret_mean=ret_mean, apply_smoother=apply_smoother)
  return(list(fit=fit, Y=Y, observations=observations))
}

best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")
best_sampled <- c("DUI")

W <- matrix(0, 3, 3)
F <- matrix(c(1, 0, 1), 3, 1)

W_scale_init <- 1
fit_models <- list()

ret_mean <- FALSE # sample from Sigma | eta
n_samples <- 2000

verbose <- TRUE
subset_time <- TRUE
apply_smoother <- TRUE

save_path <- "C:/Users/kim/Desktop/rules_of_life_stray_run1"

baboon <- "DUI"
w <- 1

for(baboon in best_sampled) {
  cat("Fitting",baboon,"over all taxa...\n")
  load(paste0("C:/Users/kim/Documents/rules_of_life/subsetted_indiv_data/",baboon,"_data.RData"))
  
  diag(W) <- c(1, 1, 1/100)
  append_str <- paste0("_W",w)
  
  fit_all <- fit_model(indiv_data, W, W_scale_init, F, n_samples=n_samples, ret_mean=ret_mean, apply_smoother=apply_smoother, subset_time=TRUE)
  fit_models[[length(fit_models)+1]] <- fit_all$fit
  fit <- fit_all$fit
  Y <- fit_all$Y
  observations <- fit_all$observations

  if(verbose & ret_mean) {
    plot_Sigma(fit, Y, baboon, save_path=paste0(save_path,"/plots"), append=append_str, as_corr=FALSE)
    #plot_Sigma(fit, Y, baboon, save_path=paste0(save_path,"/plots"), append=paste0(append_str,"_corr"), as_corr=TRUE)
  }
  if(verbose & !ret_mean) {
    plot_posterior(Y, fit, F, observations, baboon, show_filtered=FALSE, save_path=paste0(save_path,"/plots"), append_str=append_str, lr_idx=1)
    #plot_posterior(Y, fit, F, observations, baboon, show_filtered=TRUE, save_path=paste0(save_path,"/plots"), append_str=paste0(append_str,"_filter"))
  }
  if(verbose & !ret_mean) {
    #record_time(fit, baboon, save_path=paste0(save_path,"/time"), append=append_str)
  }
}

if(verbose & n_samples > 0) {
#  plot_distance_distros(fit_models, save_path)
#  plot_distant_samples(fit_models, save_path)
#  plot_element_distros(fit_models, W_scales=W_scales, percent_sample=0.1, save_path=save_path)
}

















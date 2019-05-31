library(stray)
library(phyloseq)
library(ggplot2)

rm(list=ls())

best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")
best_sampled <- c("ACA")

Q <- 3
omega <- 2*pi/365
G <- matrix(c(cos(omega), -sin(omega), 0, sin(omega), cos(omega), 0, 0, 0, 1), 3, 3)
F <- matrix(c(1, 0, 1), Q, 1)
W <- matrix(0, Q, Q)
diag(W) <- c(0.1, 0.1, 0.1)
C0 <- W*10
gamma <- 1
smooth <- TRUE
plot_filter <- FALSE
subset_time <- FALSE
n_samples <- 2000
save_img <- TRUE

for(baboon in best_sampled) {
  for(full in c(TRUE)) {
    cat("Fitting",baboon,"over all taxa =",full,"...\n")
    if(full) {
      load(paste0("C:/Users/kim/Documents/rules_of_life/subsetted_indiv_data/",baboon,"_data.RData"))
    } else {
      load(paste0("C:/Users/kim/Documents/rules_of_life/subsetted_indiv_data/",baboon,"_subset_data.RData"))
    }
  
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
  
    D <- nrow(Y)
    # set Fourier coefficients uniformly at 1; set offset to mean for this logratio over all timepoints
    M0 <- matrix(1, Q, D-1)
    M0[3,] <- alr_means
    upsilon <- D+10
    Xi <- diag(D-1)
    
    #fit <- labraduck(Y=Y, upsilon=upsilon, Xi=Xi,
    #                 gamma=gamma, F=F, G=G, W=W, M0=M0, C0=C0, observation=observations,
    #                 max_iter=100000, step_size=0.002, eps_f=1e-11, b1=0.99,
    #                 optim_method="adam", decomp_method="eigen", useSylv=FALSE, smooth=smooth)
    fit <- labraduck(Y=Y, upsilon=upsilon, Xi=Xi,
                     gamma=gamma, F=F, G=G, W=W, M0=M0, C0=C0, observation=observations,
                     max_iter=100000, decomp_method="eigen", n_samples=n_samples, ret_mean=FALSE, smooth=smooth)
    sample_no <- 1
  
    # does a sample from Sigma look ok in scale?
    if(full) {
      if(save_img) {
      png(paste0("C:/Users/kim/Desktop/rules_of_life_stray_run1/plots_full/",baboon,"_Sigma.png"))
      image(fit$Sigma[,,sample_no])
      dev.off()
      # compare to empirical covariance over taxa
      png(paste0("C:/Users/kim/Desktop/rules_of_life_stray_run1/plots_full/",baboon,"_empcov.png"))
      image(cov(driver::alr(t(Y + 0.5))))
      dev.off()
      } else {
        image(fit$Sigma[,,sample_no])
      }
    }
    
    # plot the first sample of (1) log-transformed Y (2) approximated eta (3) filtered Theta (with F applied)
    # these should be similar
    T <- max(observations)
    N <- length(observations)
    filtered_ymin <- rep(Inf, T)
    filtered_ymax <- rep(-Inf, T)
    smoothed_ymin <- rep(Inf, T)
    smoothed_ymax <- rep(-Inf, T)
    Ft <- t(F)
    lr <- driver::alr(t(Y + 0.5))
    lr_idx <- 1
    for(s in 1:n_samples) {
      ThetasF_1T <- fit$Thetas_filtered[,s]
      dim(ThetasF_1T) <- c(Q, D-1, T)
      if(smooth) {
        ThetasS_1T <- fit$Thetas_smoothed[,s]
        dim(ThetasS_1T) <- c(Q, D-1, T)
      }
      for(t in 1:T) {
        Theta_f_t <- ThetasF_1T[,,t]
        filtered_pt <- (Ft%*%Theta_f_t)[lr_idx]
        if(filtered_pt < filtered_ymin[t]) {
          filtered_ymin[t] <- filtered_pt
        }
        if(filtered_pt > filtered_ymax[t]) {
          filtered_ymax[t] <- filtered_pt
        }
        if(smooth) {
          Theta_s_t <- ThetasS_1T[,,t]
          smoothed_pt <- (Ft%*%Theta_s_t)[lr_idx]
          if(smoothed_pt < smoothed_ymin[t]) {
            smoothed_ymin[t] <- smoothed_pt
          }
          if(smoothed_pt > smoothed_ymax[t]) {
            smoothed_ymax[t] <- smoothed_pt
          }
        }
      }
    }
    df <- data.frame(timepoint=as(observations, "vector"), alrY=lr[,lr_idx])
    df2 <- data.frame(timepoint=as(observations, "vector"), eta_hat=fit$Eta[lr_idx,,sample_no])
    df <- merge(df, df2, by='timepoint', all=TRUE)
    if(plot_filter) {
      df2 <- data.frame(timepoint=1:T, ymin_f=filtered_ymin, ymax_f=filtered_ymax)
      df <- merge(df, df2, by='timepoint', all=TRUE)
    }
    if(smooth) {
      df2 <- data.frame(timepoint=1:T, ymin_s=smoothed_ymin, ymax_s=smoothed_ymax)
      df <- merge(df, df2, all=TRUE)
    }
    p <- ggplot(df, aes(timepoint))
    if(plot_filter) {
      p <- p + geom_ribbon(aes(ymin=ymin_f, ymax=ymax_f), fill = "grey90")
    }
    if(smooth) {
      p <- p + geom_ribbon(aes(ymin=ymin_s, ymax=ymax_s), fill = "grey80")
    }
    p <- p + geom_point(aes(x=timepoint, y=alrY), color="red") +
      geom_point(aes(x=timepoint, y=eta_hat), color="blue") +
      theme_minimal()
    # FILTERING DISTRIBUTION FALLS TOWARD NEGATIVE FOR SOME LOG RATIOS, SOME INDIVIDUALS ONLY???
    # e.g. THR x 1 versus THR x 25
    # this definitely has something to do with MORE abundant vs less!!!
    p <- p + ylab("ALR(Bifidobacteriaceae/Helicobacteraceae)")
    if(save_img) {
      if(full) {
        width <- 8
        if(subset_time) {
          width <- 5
        }
        ggsave(filename=paste0("C:/Users/kim/Desktop/rules_of_life_stray_run1/plots_full/",baboon,"_Theta_smoothed.png"),
               units="in", scale=2, width=width, height=2)
      } else {
        ggsave(filename=paste0("C:/Users/kim/Desktop/rules_of_life_stray_run1/plots_subset/",baboon,"_Theta_smoothed.png"),
                               units="in", scale=2, width=5, height=2)
      }
      if(full) {
        sink(paste0("C:/Users/kim/Desktop/rules_of_life_stray_run1/time_breakdowns_full/",baboon,".txt"), append=FALSE, split=FALSE)
        for(i in 1:length(fit$Timer)) {
          cat(names(fit$Timer)[i]," -- ",round(fit$Timer[i],1),"\n")
        }
        sink()
      }
    }
  }
}
















library(stray)
library(phyloseq)
library(ggplot2)

rm(list=ls())

best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")

Q <- 2
omega <- 2*pi/365
G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)
F <- matrix(c(1, 0), 2, 1)
W <- diag(Q)*0.5
C0 <- W*10
gamma <- 1
smooth <- TRUE
reference_taxon <- 9 # low count group

for(baboon in best_sampled) {
  for(full in c(FALSE, TRUE)) {
    cat("Fitting",baboon,"over all taxa =",full,"...\n")
    if(full) {
      load(paste0("C:/Users/kim/Desktop/temp/",baboon,"_data.RData"))
    } else {
      load(paste0("C:/Users/kim/Desktop/temp/",baboon,"_subset_data.RData"))
    }
  
    Y_full <- indiv_data$ys
    observations_full <- indiv_data$observation_vec
  
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
  
    D <- nrow(Y)
    M0 <- matrix(rnorm(Q*(D-1)), Q, D-1)
    upsilon <- D+10
    Xi <- diag(D-1)
    
    #fit <- labraduck(Y=Y, upsilon=upsilon, Xi=Xi,
    #                 gamma=gamma, F=F, G=G, W=W, M0=M0, C0=C0, observation=observations,
    #                 max_iter=100000, step_size=0.002, eps_f=1e-11, b1=0.99,
    #                 optim_method="adam", decomp_method="eigen", useSylv=FALSE, smooth=smooth)
    fit <- labraduck(Y=Y, upsilon=upsilon, Xi=Xi,
                     gamma=gamma, F=F, G=G, W=W, M0=M0, C0=C0, observation=observations,
                     max_iter=100000, decomp_method="eigen", smooth=smooth)
    sample_no <- 1
  
    # does a sample from Sigma look ok in scale?
    # fit$Sigma[,,1]
    if(full) {
      png(paste0("C:/Users/kim/Desktop/rules_of_life_stray_run1/plots_full/",baboon,"_Sigma.png"))
      image(fit$Sigma[,,1])
      dev.off()
      # compare to empirical covariance over taxa
      png(paste0("C:/Users/kim/Desktop/rules_of_life_stray_run1/plots_full/",baboon,"_empcov.png"))
      image(cov(driver::alr(t(Y + 0.5))))
      dev.off()
    }
    
    # plot the first sample of (1) log-transformed Y (2) approximated eta (3) filtered Theta (with F applied)
    # these should be similar
    T <- max(observations)
    N <- length(observations)
    filtered_pts <- numeric(T)
    smoothed_pts <- numeric(T)
    Ft <- t(F)
    lr <- driver::alr(t(Y + 0.5))
    if(full) {
      lr_idx <- 20
    } else {
      lr_idx <- 1
    }
    for(t in 1:T) {
      Theta_t <- fit$Thetas_filtered[,t] # always returns first sample
      dim(Theta_t) <- c(Q, D-1) # system_dim x D-1
      filtered_pts[t] <- (Ft%*%Theta_t)[lr_idx]
      if(smooth) {
        Theta_t <- fit$Thetas_smoothed[,t] # always returns first sample
        dim(Theta_t) <- c(Q, D-1) # system_dim x D-1
        smoothed_pts[t] <- (Ft%*%Theta_t)[lr_idx]
      }
    }
    # df <- data.frame(timepoint=1:T, which="filtered", value=filtered_pts)
    df <- data.frame(timepoint=as(observations, "vector"), estimator=rep("alr(Y)", N), value=lr[,lr_idx])
    df <- rbind(df, data.frame(timepoint=as(observations, "vector"), estimator=rep("eta_hat", N), value=fit$Eta[lr_idx,,sample_no]))
    if(smooth) {
      df <- rbind(df, data.frame(timepoint=1:T, estimator=rep("smoothed", T), value=smoothed_pts))
    }
    # alternative line (smoothed) and points (alr(Y) and eta_hat) -- not so pretty but might want to return to
    # p <- ggplot(df, aes(x=timepoint, y=value, which="smoothed")) + geom_line()
    # p <- p + geom_point(aes(x=timepoint, y=value, color="red"), subset(df, df$which=="alr(Y)"))
    # p <- p + geom_point(aes(x=timepoint, y=value, color="blue"), subset(df, df$which=="eta_hat"))
    p <- ggplot(df, aes(x=timepoint, y=value, color=estimator)) + geom_line()
    if(full) {
      p <- p + ylab("ALR(Atopobiaceae/Helicobacteraceae)")
    } else {
      p <- p + ylab("ALR(Bifidobacteriaceae/Helicobacteraceae)")
    }
    if(full) {
      ggsave(filename=paste0("C:/Users/kim/Desktop/rules_of_life_stray_run1/plots_full/",baboon,"_Theta_smoothed.png"),
                             units="in", scale=2, width=5, height=2)
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
















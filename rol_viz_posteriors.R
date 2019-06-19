args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  stop("Testing usage: Rscript rol_viz_posteriors.R TRUE", call.=FALSE)
}

use_Riemann <- as.logical(args[1])

devtools::load_all("/data/mukherjeelab/labraduck")
source("rol_includes.R")

n_samples <- 100 # eta samples to draw; Kalman filter, simulation smoother run on these

n_samples_subset <- 100

n_indiv <- length(best_sampled)
all_samples <- matrix(NA, D-1, (D-1)*n_samples_subset*n_indiv)
labels <- c()
for(i in 1:n_indiv) {
  load(paste0(data_path,"/",best_sampled[i],"_fit.RData"))
  dim(Sigma_samples) <- c((D-1), (D-1), n_samples)
  Sigma_samples <- Sigma_samples[,,1:n_samples_subset]
  all_samples[,((i-1)*(D-1)*n_samples_subset+1):(i*(D-1)*n_samples_subset)] <- Sigma_samples
  labels <- c(labels, rep(best_sampled[i], n_samples_subset))
}

distance_mat <- matrix(NA, n_samples_subset*n_indiv, n_samples_subset*n_indiv)
for(i in 1:(n_indiv*n_samples_subset)) {
  for(j in i:(n_indiv*n_samples_subset)) {
    i_idx <- (i-1)*(D-1)
    A <- all_samples[,(i_idx+1):(i_idx+(D-1))]
    j_idx <- (j-1)*(D-1)
    B <- all_samples[,(j_idx+1):(j_idx+(D-1))]
    distance_mat[i,j] <- mat_dist(A, B, use_Riemann=use_Riemann)
    distance_mat[j,i] <- distance_mat[i,j]
  }
}

fit <- cmdscale(distance_mat, eig=TRUE, k=2) # k is the number of dim
cat("Lambda 1:",fit$eig[1],"\n")
cat("Lambda 2:",fit$eig[2],"\n")
cat("Lambda 3:",fit$eig[3],"\n")

df <- data.frame(x=fit$points[,1], y=fit$points[,2], labels=labels)
p <- ggplot(df, aes(x=x, y=y, color=labels)) + geom_point()
if(use_Riemann) {
  p <- p + ggtitle("Riemannian distance")
} else {
  p <- p + ggtitle("Frobenius norm of difference")
}
plot_save_name <- "Sigma_posterior_ordination_"
if(use_Riemann) {
  plot_save_name <- paste0(plot_save_name,"Riemann.png")
} else {
  plot_save_name <- paste0(plot_save_name,"Frobenius.png")
}
ggsave(paste0(save_path,"/",plot_save_name), scale=2,
         width=4, height=4, units="in", dpi=100)






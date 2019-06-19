library(matrixsampling)

# compare behavior of Frobenius norm of difference and Riemannian distance on some contrived examples
# this should give some intuition about behavior

# we'll test 6 combinations
# diagonal matrix:            loosely centered, tightly centered, tightly centered with large scale
# positive off-diagonal mass: loosely centered, tightly centered, tightly centered with large scale

draw_samples <- function(upsilon, Xi, N) {
  samples <- rinvwishart(N, upsilon, Xi)
  dim(samples) <- c(dim(Xi)[1], dim(Xi)[1]*N)
  return(samples)
}

calc_dist_mat <- function(n_tests, D, N, all_samples, use_Riemann) {
  distance_mat <- matrix(NA, N*n_tests, N*n_tests)
  for(i in 1:(n_tests*N)) {
    for(j in i:(n_tests*N)) {
      i_idx <- (i-1)*D
      A <- all_samples[,(i_idx+1):(i_idx+D)]
      j_idx <- (j-1)*D
      B <- all_samples[,(j_idx+1):(j_idx+D)]
      distance_mat[i,j] <- mat_dist(A, B, use_Riemann=use_Riemann)
      distance_mat[j,i] <- distance_mat[i,j]
    }
  }
  return(distance_mat)
}

ordinate_plot <- function(n_tests, D, N, all_samples, all_labels, use_Riemann, save_tag) {
  distance_mat <- calc_dist_mat(n_tests, D, N, all_samples, use_Riemann)

  fit <- cmdscale(distance_mat, eig=TRUE, k=2) # k is the number of dim
  cat("Top three eigenvalues:",fit$eig[1],",",fit$eig[2],",",fit$eig[3],"\n")

  df <- data.frame(x=fit$points[,1], y=fit$points[,2], labels=all_labels)
  p <- ggplot(df, aes(x=x, y=y, color=all_labels)) + geom_point()
  if(use_Riemann) {
    p <- p + ggtitle("Riemannian distance")
    ggsave(paste0("test_ordination_",save_tag,"_Rie.png"), scale=2, width=4, height=4, units="in", dpi=100)
  } else {
    p <- p + ggtitle("Frobenius norm of difference")
    ggsave(paste0("test_ordination_",save_tag,"_Fro.png"), scale=2, width=4, height=4, units="in", dpi=100)
  }
}

devtools::load_all("/data/mukherjeelab/labraduck")
source("rol_includes.R")

D <- 20
N <- 100

cat("Plotting diagonal visualization...\n")

diag_labels <- c(rep("diag_loose", N), rep("diag_tight", N), rep("diag_tight_large", N))
diag_samples <- matrix(NA, D, D*N*3)
diag_samples[,1:(1*D*N)] <- draw_samples(D+10, diag(D)*(10-1), N)
diag_samples[,((D*N*1)+1):(2*D*N)] <- draw_samples(D+100, diag(D)*(100-1), N)
diag_samples[,((D*N*2)+1):(3*D*N)] <- draw_samples(D+100, diag(D)*(100-1)*10, N)

ordinate_plot(n_tests=3, D, N, diag_samples, diag_labels, use_Riemann=TRUE, "diag")
ordinate_plot(n_tests=3, D, N, diag_samples, diag_labels, use_Riemann=FALSE, "diag")

cat("Plotting positive dense visualization...\n")

pos_labels <- c(rep("pos_loose", N), rep("pos_tight", N), rep("pos_tight_large", N))
pos_samples <- matrix(NA, D, D*N*3)
off_diag <- matrix(0.5, D, D)
diag(off_diag) <- 0
pos_samples[,1:(1*D*N)] <- draw_samples(D+10, (diag(D)+off_diag)*(10-1), N)
pos_samples[,((D*N*1)+1):(2*D*N)] <- draw_samples(D+100, (diag(D)+off_diag)*(100-1), N)
pos_samples[,((D*N*2)+1):(3*D*N)] <- draw_samples(D+100, (diag(D)+off_diag)*(100-1)*10, N)

ordinate_plot(n_tests=3, D, N, pos_samples, pos_labels, use_Riemann=TRUE, "pos")
ordinate_plot(n_tests=3, D, N, pos_samples, pos_labels, use_Riemann=FALSE, "pos")

cat("Plotting ALL visualization...\n")

all_samples <- matrix(NA, D, D*N*6)
all_samples[,1:(D*N*3)] <- diag_samples
all_samples[,(D*N*3+1):(D*N*6)] <- pos_samples
all_labels <- c(diag_labels, pos_labels)
ordinate_plot(n_tests=6, D, N, all_samples, all_labels, use_Riemann=TRUE, "all")
ordinate_plot(n_tests=6, D, N, all_samples, all_labels, use_Riemann=FALSE, "all")

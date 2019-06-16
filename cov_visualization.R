library(matrixsampling)
library(ggplot2)
library(Rcpp)
library(RcppEigen)

sourceCpp("src/cov_viz_test.cpp")

# note: ln(AB) = ln(A) + ln(B) only true for matrices that commute
# this is a formula for matrices that don't commute but it's complex
# (in other words, we can simplify the Riemannian distance)
mat_log <- function(A) {
  eigen.A <- eigen(A)
  log.D.mat <- diag(length(eigen.A$values))
  diag(log.D.mat) <- log(eigen.A$values)
  return(eigen.A$vectors%*%log.D.mat%*%solve(eigen.A$vectors))
}

mat_exp <- function(A) {
  eigen.A <- eigen(A)
  log.D.mat <- diag(length(eigen.A$values))
  diag(log.D.mat) <- exp(eigen.A$values)
  return(eigen.A$vectors%*%log.D.mat%*%solve(eigen.A$vectors))
}

fro_dist <- function(A, B) {
  mat_diff <- A - B
  return(norm(mat_diff, type="F"))
}

riemannian_dist <- function(A, B) {
  U <- chol(A)
  U.inv <- solve(U)
  X <- t(U.inv)%*%B%*%U.inv
  eig.X <- eigen(X)
  return(sqrt(sum((log(eig.X$values))^2)))
}

N <- 100 # samples
D <- 10 # taxa

all_samples <- matrix(NA, D, D*N*4)

Xi.base <- diag(D)

# loosely centered diagonal (W)
upsilon <- D+10
Xi <- Xi.base*(upsilon - D - 1)
temp <- rinvwishart(N, upsilon, Xi)
dim(temp) <- c(D, D*N)
all_samples[,1:(D*N)] <- temp

# tightly centered diagonal (X)
upsilon <- D+100
temp <- rinvwishart(N, upsilon, Xi)
dim(temp) <- c(D, D*N)
all_samples[,(D*N+1):(2*D*N)] <- temp

Xi.base <- diag(D)*0.5 + 0.5

# loosely centered dense (Y)
upsilon <- D+10
Xi <- Xi.base*(upsilon - D - 1)
temp <- rinvwishart(N, upsilon, Xi)
dim(temp) <- c(D, D*N)
all_samples[,(2*D*N+1):(3*D*N)] <- temp

# tightly centered dense (Z)
upsilon <- D+100
temp <- rinvwishart(N, upsilon, Xi)
dim(temp) <- c(D, D*N)
all_samples[,(3*D*N+1):(4*D*N)] <- temp

# big cross-distance matrix
all_d <- matrix(0, N*4, N*4)
for(i in 1:(4*N-1)) {
  for(j in i:(4*N)) {
    i_idx <- (i-1)*D
    A <- all_samples[,(i_idx+1):(i_idx+D)]
    j_idx <- (j-1)*D
    B <- all_samples[,(j_idx+1):(j_idx+D)]
    #all_d[i,j] <- mat_dist(A, B, use_Riemann=FALSE)
    all_d[i,j] <- mat_dist(A, B, use_Riemann=TRUE)
    all_d[j,i] <- all_d[i,j]
  }
}

fit <- cmdscale(all_d, eig=TRUE, k=2) # k is the number of dim
cat("Lambda 1:",fit$eig[1],"\n")
cat("Lambda 2:",fit$eig[2],"\n")
cat("Lambda 3:",fit$eig[3],"\n")

df <- data.frame(x=fit$points[,1], y=fit$points[,2],
                 labels=c(rep("W",N), rep("X",N), rep("Y",N), rep("Z",N)))
p <- ggplot(df, aes(x=x, y=y, color=labels)) +
  geom_point() +
  ggtitle("Frobenius diff distance test case") +
  theme_minimal()
#ggsave("C:/Users/kim/Desktop/Frobeniusdiff_test.png", scale=2,
#       width=4, height=4, units="in", dpi=100)
p







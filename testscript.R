library(LaplacesDemon)
library(driver)
library(stray)
library(ggplot2)

# need to look at the simulation smoother -- bugs? fits look suprisingly shitty

rm(list=ls()) # really seems necessary to flush; should check on this

D <- 5
N <- 20

add_gaps <- TRUE
periodic <- TRUE
smooth <- FALSE

# simulate data
Y <- matrix(0, D, N)
for(i in 1:N) {
  props <- LaplacesDemon::rdirichlet(1, rep(0.2, D))
  Y[,i] <- rmultinom(1, 5000, props)
}

if(periodic)  {
  Q <- 2
  omega <- 2*pi/5
  G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)
  F <- matrix(c(1, 0), 2, 1)
  W <- diag(Q)*0.1
  M0 <- matrix(rnorm(Q*(D-1)), Q, D-1)
  C0 <- W*10
  upsilon <- D+10
  Xi <- diag(D-1)
  gamma <- 1
} else {
  Q <- 1
  G <- diag(Q)
  F <- matrix(1, Q, 1)
  W <- diag(Q)*0.5
  M0 <- matrix(rnorm(Q*(D-1)), Q, D-1)
  C0 <- W
  upsilon <- D+10
  Xi <- diag(D-1)
  gamma <- 1
}

if(add_gaps) {
  observations <- matrix(1:N, nrow=1)
  halfway <- round(N/2)
  observations[(halfway+1):N] <- observations[(halfway+1):N] + 10 # gap of 10 observations
} else {
  observations <- matrix(as.numeric(1:N), nrow=1)
}

fit <- labraduck(Y=Y, upsilon=upsilon, Xi=Xi,
          gamma=gamma, F=F, G=G, W=W, M0=M0, C0=C0, observation=observations,
          max_iter=100000, step_size=0.002, eps_f=1e-11, b1=0.99,
          optim_method="adam", decomp_method="eigen", useSylv=FALSE, smooth=smooth)

sample_no <- 1

# does a sample from Sigma look ok in scale?
fit$Sigma[,,1]
image(fit$Sigma[,,1])

# compare to empirical covariance over taxa
image(cov(driver::alr(t(Y + 0.5))))

# plot the first sample of (1) log-transformed Y (2) approximated eta (3) filtered Theta (with F applied)
# these should be similar
T <- max(observations)
filtered_pts <- numeric(T)
smoothed_pts <- numeric(T)
Ft <- t(F)
lr <- driver::alr(t(Y + 0.5))
lr_idx <- 1
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
df <- data.frame(timepoint=1:T, which="filtered", value=filtered_pts)
df <- rbind(df, data.frame(timepoint=as(observations, "vector"), which=rep("alr(Y)", N), value=lr[,lr_idx]))
df <- rbind(df, data.frame(timepoint=as(observations, "vector"), which=rep("eta_hat", N), value=fit$Eta[lr_idx,,sample_no]))
if(smooth) {
  df <- rbind(df, data.frame(timepoint=1:T, which=rep("smoothed", T), value=smoothed_pts))
}
p <- ggplot(df, aes(x=timepoint, y=value, color=which)) + geom_line()
p


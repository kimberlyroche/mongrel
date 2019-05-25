library(Rcpp)
library(RcppEigen)
library(matrixsampling)
library(mvtnorm)

# powers G through eigenvalue decomposition
stack_G <- function(G, it_begin, it_end, descending=TRUE, transpose=FALSE) {
  obj <- G
  if(transpose) {
    obj <- t(G)
  }
  obj <- eigen(obj)
  e_vec <- obj$vectors
  e_val <- diag(2)
  diag(e_val) <- obj$values
  if(it_begin != it_end) {
    if(descending) {
      if(it_begin > it_end) {
        power_it <- length(1:(it_begin-it_end+1))
        e_val <- e_val**power_it
      } else {
        # invalid case
        e_val <- matrix(0, 2, 2)
      }
    } else {
      if(it_begin < it_end) {
        power_it <- length(1:(it_end-it_begin+1))
        e_val <- e_val**power_it
      } else {
        # invalid case
        e_val <- matrix(0, 2, 2)
      }
    }
  }
  # explicitly only returning the real part of A
  # some tiny complex eigenvalues can be produced -- cool to truncate in this way?
  ret_val <- Re(e_vec%*%e_val%*%solve(e_vec))
  return(ret_val)
}

build_A <- function(T, gamma, F, G, W, C0) {
  # calculate A (covariance matrix over states) for this simulation
  # this is the expression exactly as in the manuscript, calculated from Cov(eta_t, eta_{t-k})
  A <- matrix(0, T, T)
  for(i in 1:T) {
    for(j in 1:T) {
      if(i == j) {
        # diagonal
        t <- j
        first_sum <- matrix(0, 2, 2)
        if(t >= 2) {
          for(ell in t:2) {
            G_left <- stack_G(G, t, ell)
            G_right <- stack_G(G, ell, t, descending=FALSE, transpose=TRUE)
            addend <- G_left%*%W%*%G_right
            first_sum <- first_sum + addend
          }
        }
        # second sum
        G_left <- stack_G(G, t, 1)
        G_right <- stack_G(G, 1, t, descending=FALSE, transpose=TRUE)
        second_sum <- G_left%*%C0%*%G_right
        A[t,t] <- gamma + F%*%(W + first_sum + second_sum)%*%t(F)
      } else {
        tk <- i
        t <- j
        if(j < i) {
          tk <- j
          t <- i
        }
        # off-diagonal
        first_sum <- matrix(0, 2, 2)
        if(tk >= 2) {
          for(ell in tk:2) {
            G_left <- stack_G(G, t, ell)
            G_right <- stack_G(G, ell, tk, descending=FALSE, transpose=TRUE)
            first_sum <- first_sum + G_left%*%W%*%G_right
          }
        }
        G_left <- stack_G(G, t, 1)
        G_right <- stack_G(G, 1, tk, descending=FALSE, transpose=TRUE)
        second_sum <- G_left%*%C0%*%G_right
        G_left <- stack_G(G, t, tk+1)
        A[i,j] <- F%*%(G_left%*%W + first_sum + second_sum)%*%t(F)
      }
    }
  }
  if(min(eigen(A)$values) < 0) {
    cat("Matrix A has negative eigenvalue(s)!\n")
  }
  return(A)
}

# test R vs. C++ implementations of matrix A construction

sourceCpp("test.cpp")

T <- 10
omega <- 2*pi/5
G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)
F <- matrix(c(1, 0), 2, 1)
W <- diag(2)
M0 <- matrix(rnorm(10), 2, 5)
C0 <- W
A <- dlm_A(T, 1, F, G, W, C0)
A_R <- build_A(T, 1, t(F), G, W, C0) # expects transposed F!

B <- dlm_B(F, G, M0, T)

test <- rmatrixnormal(1, B, diag(5), res)[,,1]

# Kalman filter
# data_obj is a list containing ys, F, W, G, upsilon, Xi, gamma, Sigma, M.0, C.0
# observation_vec (if present) indicates the spacing of observations, e.g. c(1, 3, 4, 7)
#   indicates observations 2, 5, & 6 are missing and should be imputed in the usual way
fit_filter <- function(data_obj, observation_vec=NULL) {
  D <- ncol(data_obj$ys)
  Theta.dim <- ncol(data_obj$G)
  T <- nrow(data_obj$ys)
  if(!is.null(observation_vec)) {
    T <- max(observation_vec)
  }
  upsilon.t <- data_obj$upsilon
  Xi.t <- data_obj$Xi
  M.t <- data_obj$M.0
  C.t <- data_obj$C.0
  Thetas.t <- array(0, dim=c(Theta.dim, D, T)) # sample at each t
  Cs.t <- array(0, dim=c(Theta.dim, Theta.dim, T))
  Ms.t <- array(0, dim=c(Theta.dim, D, T))
  Rs.t <- array(0, dim=c(Theta.dim, Theta.dim, T))
  W.t <- NULL # last observed W.t (when using discount)
  for(t in 1:T) {
    R.t <- data_obj$G%*%C.t%*%t(data_obj$G) + data_obj$W
    R.t <- round(R.t, 10)
    Rs.t[,,t] <- R.t
    A.t <- data_obj$G%*%M.t
    if(!is.null(observation_vec) && !(t %in% observation_vec)) {
      # impute
      # carry forward without update
      M.t <- A.t
      Ms.t[,,t] <- M.t
      C.t <- R.t
      Cs.t[,,t] <- C.t
      # no change to Sigma.t parameters (zero change in uncertainty)
    } else {
      f.t.T <- data_obj$F%*%A.t
      q.t <- data_obj$gamma + (data_obj$F%*%R.t%*%t(data_obj$F))[1,1]
      if(is.null(observation_vec)) {
        e.t.T <- data_obj$ys[t,] - f.t.T
      } else {
        # there are missing observations but this is not one of them
        # index[1]: if there are more than one observations on this date (VERY rare) just take the first
        # need a better way of handling this
        e.t.T <- data_obj$ys[as(which(observation_vec == t)[1], "numeric"),] - f.t.T
      }
      S.t <- R.t%*%t(data_obj$F)/q.t
      M.t <- A.t + S.t%*%e.t.T
      Ms.t[,,t] <- M.t
      C.t <- R.t - q.t*S.t%*%t(S.t)
      C.t <- round(C.t, 10)
      Cs.t[,,t] <- C.t
      upsilon.t <- upsilon.t + 1
      Xi.t <- Xi.t + t(e.t.T)%*%e.t.T/q.t
    }
    Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
    LU <- t(chol(C.t))
    LV <- chol(Sigma.t)
    #Thetas.t[,,t] <- rmatrixnormal(1, M.t, C.t, Sigma.t)[,,1]
    Thetas.t[,,t] <- M.t + LU%*%matrix(rnorm(10), 2, 5)%*%LV
  }
  return(list(Thetas.t=Thetas.t, upsilon=upsilon.t, Xi=Xi.t, Ms.t=Ms.t, Cs.t=Cs.t, Rs.t=Rs.t))
}

# simulation smoother
# fit_obj is the output of fit_filter()
fit_smoother <- function(data_obj, fit_obj) {
  D <- ncol(data_obj$ys)
  Theta.dim <- ncol(data_obj$G)
  T <- dim(fit_obj$Thetas.t)[3]
  Sigma.t <- rinvwishart(1, fit_obj$upsilon, fit_obj$Xi)[,,1]
  Thetas.t.smoothed <- array(0, dim=c(Theta.dim, D, T))
  Ms.t <- array(0, dim=c(Theta.dim, D, T))
  etas.t <- matrix(0, D, T)
  rmatnorm_mean <- fit_obj$Ms.t[,,T]
  dim(rmatnorm_mean) <- c(Theta.dim, D) # fix if the state has one dimension only
  # can't do drop=FALSE here
  Thetas.t.smoothed[,,T] <- rmatrixnormal(1, rmatnorm_mean, fit_obj$Cs.t[,,T], Sigma.t)[,,1]
  etas.t[,T] <- rmatrixnormal(1, data_obj$F%*%Thetas.t.smoothed[,,T], 1, Sigma.t)[,,1]
  Ms.t[,,T] <- rmatnorm_mean
  for(t in (T-1):1) {
    Z.t <- fit_obj$Cs.t[,,t]%*%t(data_obj$G)%*%solve(fit_obj$Rs.t[,,(t+1)])
    M.t.star <- fit_obj$Ms.t[,,t] + Z.t%*%(Thetas.t.smoothed[,,(t+1)] - data_obj$G%*%fit_obj$Ms.t[,,t])
    Ms.t[,,t] <- M.t.star
    C.t.star <- round(fit_obj$Cs.t[,,t] - Z.t%*%fit_obj$Rs.t[,,(t+1)]%*%t(Z.t), 10)
    Thetas.t.smoothed[,,t] <- rmatrixnormal(1, M.t.star, C.t.star, Sigma.t)[,,1]
    # draw an eta, this is how we'll estimate modeled covariance
    etas.t[,t] <- rmatrixnormal(1, data_obj$F%*%Thetas.t.smoothed[,,t], 1, Sigma.t)[,,1]
  }
  return(list(Thetas.t=Thetas.t.smoothed, etas.t=etas.t, Ms.t=Ms.t))
}

# test random-walk ungapped Kalman filter

sourceCpp("test.cpp")

T <- 50
M0 <- matrix(rnorm(5), 2, 5)
eta <- rmvnorm(T, matrix(rnorm(5), 1, 5), diag(5)) # rows are samples

F <- matrix(c(1, 1), 2, 1)
G <- diag(2)
W <- diag(2)*0.5
V <- diag(5)*0.05
C0 <- diag(2)
upsilon <- 7
Xi <- diag(5)*0.1
observation_vec <- 1:T
observation_vec <- numeric(50)
observation_vec[1:25] <- 1:25
observation_vec[26:50] <- 51:75 # gap in the middle!

data_obj <- list(ys=eta,
                 F=t(F),
                 W=W,
                 G=G,
                 upsilon=upsilon,
                 Xi=Xi, # center on empirical covariance
                 gamma=1,
                 Sigma=NULL,
                 M.0=M0, # could randomly init too
                 C.0=C0)



# C++
res <- filter(eta, F, G, W, 1, upsilon, Xi, M0, C0, observations=observation_vec)
smooth <- simulation_smooth(eta, G, res$upsilon.T, res$Xi.T, res$Rs, res$Ms, res$Cs)

# R
resR <- fit_filter(data_obj, observation_vec=observation_vec)
smoothR <- fit_smoother(data_obj, resR)



T_actual <- max(observation_vec)

filtered_obs1 <- numeric(T_actual)
filtered_smooth1 <- numeric(T_actual)
for(t in 1:T_actual) {
  filtered_obs1[t] <- t(F)%*%res$Thetas[1:2,t]
  filtered_smooth1[t] <- t(F)%*%smooth$Thetas[1:2,t]
}
plot(observation_vec, eta[,1], type="l", ylim=c(-5,5))
lines(1:T_actual, filtered_obs1, col="blue", lty=2)
lines(1:T_actual, filtered_smooth1, col="blue")


filtered_obs1R <- numeric(T_actual)
filtered_smooth1R <- numeric(T_actual)
for(t in 1:T_actual) {
  filtered_obs1R[t] <- t(F)%*%resR$Thetas.t[1:2,1,t]
  filtered_smooth1R[t] <- t(F)%*%smoothR$Thetas.t[1:2,1,t]
}
#plot(observation_vec, eta[,1], type="l", ylim=c(-5,5))
lines(1:T_actual, filtered_obs1R, col="red", lty=2)
lines(1:T_actual, filtered_smooth1R, col="red")













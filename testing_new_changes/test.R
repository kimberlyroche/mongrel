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
C0 <- W
res <- dlm_cov(T, 1, F, G, W, C0)
res_R <- build_A(T, 1, t(F), G, W, C0) # expects transposed F!

res
res_R

image(res)
image(res_R)

# test random-walk ungapped Kalman filter

T <- 100
M0 <- matrix(rnorm(5), 2, 5)
eta <- rmvnorm(T, matrix(rnorm(5), 1, 5), diag(5)) # rows are samples

F <- matrix(c(1, 1), 2, 1)
G <- diag(2)
W <- diag(2)*0.5
V <- diag(5)*0.5
C0 <- diag(2)
upsilon <- 7
Xi <- diag(5)

Sigma_sample <- filter(eta, F, G, W, 1, upsilon, Xi, M0, C0, c(1:T))
Sigma_sample
image(Sigma_sample)













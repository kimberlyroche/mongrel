rm(list=ls())
   
library(stray)
library(matrixsampling)
library(driver)

#source("C:/Users/kim/Documents/batch16s/R/helper_functions.R")

m <- 10 # individuals
q <- 3 # batches
n <- m*q # total no. samples
d <- 20 # taxa

X <- matrix(0, m+(q-1), n)
X[1,] <- 1 # intercept
X[2:m,2:m] <- diag(m-1) # individual intercepts (offsets from #1)
X[2:m,(m+2):(2*m)] <- diag(m-1)
X[2:m,(2*m+2):(3*m)] <- diag(m-1)
X[(m+1),(m+2):(2*m)] <- 1 # batch 2 offset (from batch #1)
X[(m+2),(2*m+2):(3*m)] <- 1 # batch 3 offset (from batch #2)
image(X)

U1 <- diag(m+q-1) # what's the effect of extra variance at 1,1?
U1[1,1] <- 2
U1[m+1,m+1] <- 0
U1[m+2,m+2] <- 0
U2 <- matrix(0, m+q-1, m+q-1)
U2[m+1,m+1] <- 1
U3 <- matrix(0, m+q-1, m+q-1)
U3[m+2,m+2] <- 1

Gamma <- U1 + U2 + U3

upsilon <- (d-1)+10
Xi <- diag(d-1)*(upsilon-(d-1)-1)
Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
Sigma.sq <- chol(Sigma) # returns U of t(U)%*%U = Sigma
                        # L = t(U)

Theta <- matrix(0, d-1, nrow(X))
Lambda <- rmatrixnormal(1, Theta, Sigma, Gamma)[,,1]
Eta <- Lambda%*%X + t(Sigma.sq)%*%rmatrixnormal(1, matrix(0, d-1, n), diag(d-1), diag(n))[,,1]

Pi <- alrInv(t(Eta))
Y <- apply(Pi, 1, function(x) rmultinom(1, rpois(1, 10000), x)) # D x N

# visualize truth
#image(log(Y + 0.5))

U <- rbind(U1, U2, U3)

for(i in 1:5) {
  #Theta_test <- Theta + matrix(rnorm((d-1)*5, 0, 3), d-1, 5)
  maltipoofit <- maltipoo(Y=Y, X=X, upsilon=upsilon, Theta=Theta, U=U, Xi=Xi, n_samples=0)
  cat(maltipoofit$VCScale,"\n")
}

# visualize estimate
#image(maltipoofit$Eta[,,1])

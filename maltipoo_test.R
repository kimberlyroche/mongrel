library(matrixsampling)
library(driver)
devtools::load_all("/Users/ladlab/stray")
# simulate data from maltipoo model

D <- 20
N <- 100
# covariates will be an intercept and a batch label

X <- matrix(0, 2, N)
X[1,1:(N/2)] <- 1
X[2,((N/2)+1):N] <- 1

U1 <- matrix(0, 2, 2)
U1[1,1] <- 1
U2 <- matrix(0, 2, 2)
U2[2,2] <- 1

Gamma <- 0.1*U1 + 5*U2

upsilon <- (D-1)+100
Xi <- diag(D-1)*(upsilon-(D-1)-1)
Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
Sigma.sq <- chol(Sigma) # returns U of t(U)%*%U = Sigma
                        # L = t(U)

Theta <- matrix(0, D-1, 2)
Lambda <- rmatrixnormal(1, Theta, Sigma, Gamma)[,,1]
Eta <- Lambda%*%X + t(Sigma.sq)%*%rmatrixnormal(1, matrix(0, D-1, N), diag(D-1), diag(N))[,,1]
Pi <- alrInv(t(Eta))
Y <- apply(Pi, 1, function(x) rmultinom(1, rpois(1, 10000), x)) # D x N
logB1 <- log(Y[,1:(N/2)] + 0.5)
logB2 <- log(Y[,(N/2+1):N] + 0.5)
plot(density(logB1)) # distribution of counts for batch 1
lines(density(logB2), col="red") # distribution of counts for batch 2
cat("SD batch 1:",sd(c(logB1))," SD batch 2:",sd(c(logB2)),"\n")

# does this make sense?

U <- matrix(0, 4, 2)
U[1:2,1:2] <- U1
U[3:4,1:2] <- U2
U <- U

# maltipoo never learns different scales for these guys; what's going on?
maltipoofit <- maltipoo(Y=Y, X=X, upsilon=upsilon, Theta=Theta, U=U, Xi=Xi, n_samples=100)
maltipoofit$VCScale



labraduckfit <- function(D, N, Q, coord_system, iter=NULL,  
                       alr_base=NULL, ilr_base=NULL,
                       Eta=NULL, Lambda=NULL, Sigma=NULL, Sigma_default=NULL, 
                       Y=NULL, X=NULL, upsilon=NULL, 
                       Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                       init=NULL, names_categories=NULL, names_samples=NULL, 
                       names_covariates=NULL){
  m <- new_labraduckfit(D, N, Q, coord_system, iter, alr_base, ilr_base,
                      Eta, Lambda, Sigma, Sigma_default, 
                      Y, X, upsilon, Theta, Xi,Xi_default, Gamma, 
                      init, names_categories, names_samples, 
                      names_covariates)
  verify(m)
  return(m)
}


new_labraduckfit <- function(D, N, Q, coord_system, iter=NULL, 
                           alr_base=NULL, ilr_base=NULL,
                           Eta=NULL, Lambda=NULL,Sigma=NULL, Sigma_default=NULL, 
                           Y=NULL, X=NULL, upsilon=NULL, 
                           Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                           init=NULL, names_categories=NULL, names_samples=NULL, 
                           names_covariates=NULL){
  m <- new_pibblefit(D, N, Q, coord_system, iter, alr_base, ilr_base,
                      Eta, Lambda, Sigma, Sigma_default, 
                      Y, X, upsilon, Theta, Xi,Xi_default, Gamma, 
                      init, deltainit, names_categories, names_samples, 
                      names_covariates)
  class(m) <- c("labraduckfit", "pibblefit")
}

# silent for now
verify.labraduckfit <- function(m,...){
  # verify.labraduckfit(m)
  # stopifnot(is.integer(m$P))
  # P <- m$P
  # ifnotnull(m$VCScale, check_dims(m$VCScale, P, "VCScale"))
  # ifnotnull(m$U, check_dims(m$U, c(P*m$Q, m$Q), "U"))
  # ifnotnull(m$ellinit, check_dims(m$P, P, "ellinit"))
}

# reinstate later
req.labraduckfit <- function(m, r){
  # present <- sapply(m[r], is.null)
  # if(any(present)){
  #   stop("maltipoofit object does not contain required components:", r[present])
  # }
}
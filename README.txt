State this way left in (6/19/19):

Fits on 10 best sampled individuals with periodic DLM (including local level) look really bad. The interpolation is too 
strictly periodic. There is a lot of structure in this model considering how noisy the data appears. Shifting gears to 
test a basset (GP) fit with a convolution of kernels.

rol_includes.R		all utility functions
rol_fit_labraduck.R	parses data, generates priors, etc. and calls the fit function in labraduck
rol_viz_posteriors.R	plots ordination (classical MDS) of Frobenius norm diff or Riemannian distance of Sigma 
                        posteriors between and within individuals
test_viz_posteriors.R	generates some toy covariance matrices and visualizes using both Frobenius norm of difference 
                        and Riemannian distance; an intuition-builder

State of labraduck:

Currently working well with time-invariant: F, G, W

Looks like there might be a bug with optizing either W or G while the other is fixed? This was certainly working at one 
point but last time I tried this, I was getting some weird fits (NaN popping up)? Need to re-test.

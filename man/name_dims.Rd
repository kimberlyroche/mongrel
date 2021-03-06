% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics_s3.R, R/strayfit_methods.R
\name{name_dims}
\alias{name_dims}
\alias{names_covariates}
\alias{names_samples}
\alias{names_categories}
\alias{names_coords}
\alias{names_covariates<-}
\alias{names_samples<-}
\alias{names_categories<-}
\alias{names_covariates.pibblefit}
\alias{names_samples.pibblefit}
\alias{names_categories.pibblefit}
\alias{names_coords.pibblefit}
\alias{names_covariates<-.pibblefit}
\alias{names_samples<-.pibblefit}
\alias{names_categories<-.pibblefit}
\title{Generic method for getting and setting dimension names of fit object}
\usage{
names_covariates(m)

names_samples(m)

names_categories(m)

names_coords(m)

names_covariates(m) <- value

names_samples(m) <- value

names_categories(m) <- value

\method{names_covariates}{pibblefit}(m)

\method{names_samples}{pibblefit}(m)

\method{names_categories}{pibblefit}(m)

\method{names_coords}{pibblefit}(m)

\method{names_covariates}{pibblefit}(m) <- value

\method{names_samples}{pibblefit}(m) <- value

\method{names_categories}{pibblefit}(m) <- value
}
\arguments{
\item{m}{object}

\item{value}{character vector (or NULL)}
}
\description{
Generic method for getting and setting dimension names of fit object
}
\details{
\code{names_coords} is different than \code{names_categories}. 
\code{names_categories} provides access to the basic names of each multinomial 
category. In contrast, \code{names_coords} provides access to the 
names of the coordinates in which an object is represented. These coordinate
names are based on the category names. For example, category names may be, 
(OTU1, ..., OTUD) where as coordinate names could be (log(OTU1/OTUD), etc...)
if object is in default coordinate system.
}

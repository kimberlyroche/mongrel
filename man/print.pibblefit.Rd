% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strayfit_methods.R
\name{print.pibblefit}
\alias{print.pibblefit}
\title{Print dimensions and coordinate system information for pibblefit object.}
\usage{
\method{print}{pibblefit}(x, summary = FALSE, ...)
}
\arguments{
\item{x}{an object of class pibblefit}

\item{summary}{if true also calculates and prints summary}

\item{...}{other arguments to pass to summary function}
}
\description{
Print dimensions and coordinate system information for pibblefit object.
}
\examples{
\dontrun{
fit <- pibble(Y, X)
print(fit)
}
}
\seealso{
\code{\link{summary.pibblefit}} summarizes posterior intervals
}

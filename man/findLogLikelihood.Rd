% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimisation.R
\name{findLogLikelihood}
\alias{findLogLikelihood}
\title{Find Log Likelihood}
\usage{
findLogLikelihood(data, densities, estimates)
}
\arguments{
\item{data}{Data frame with columns for Length and Age, including unknown values.}

\item{densities}{Data frame containing Gaussian probability densities of observing each observed fish length under the estimated mu and sigma.}

\item{estimates}{Data frame containing the current estimates for mu, sigma, and lamba for Age group 1, 2, and 3}
}
\value{
Log likelihood of the current estiamtes for mu, sigma, and lambda belonging to the given Fish Length data.
}
\description{
Find Log Likelihood
}
\examples{
known <- docExampleData[!is.na(docExampleData$Age), ] # known data is not NA
unknown <- docExampleData[is.na(docExampleData$Age), ] # unknown data is NA
sorted_data <- rbind(known, unknown) # combine data
init <- initialise(known, unknown, sorted_data)

exp <- expector(known, sorted_data, init)
estimates <- maximiser(sorted_data, exp$posteriors)

findLogLikelihood(docExampleData, exp$densities, estimates)
}

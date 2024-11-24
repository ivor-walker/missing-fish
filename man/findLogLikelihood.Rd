\name{findLogLikelihood}
\alias{findLogLikelihood}
\title{Find Log Likelihood}
\description{Returns the log likelihood of the current estiamtes for mu, sigma, and lambda belonging to the given Fish Length data.}
\usage{findLogLikelihood(data, densities, estimates)}
\arguments{
  \item{data}{Dataframe with columns for Length and Age, including unknown values.}
  \item{densities}{Dataframe containing Gaussian probability densities of observing each observed fish length under the estimated mu and sigma.}
  \item{estimates}{Dataframe containing the current estimates for mu, sigma, and lamba for Age group 1, 2, and 3}
}
\value{Log likelihood of the current estiamtes for mu, sigma, and lambda belonging to the given Fish Length data.}
\examples{
data <- read_csv("data.csv")

known <- data[!is.na(data$Age), ] # known data is not NA
unknown <- data[is.na(data$Age), ] # unknown data is NA
sorted_data <- rbind(known, unknown) # combine data
init <- initialise(knwon, unknown, sorted_data)

exp <- expector(known, sorted_data, init)
estimates <- maximiser(sorted_data, exp$posteriors)

findLogLikelihood(data, exp$densities, estimates)
}

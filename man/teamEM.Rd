\name{teamEM}
\alias{teamEM}
\title{Perform EM Optimisation}
\description{
Completes the EM Optimisation algorithm given FishLength data with some missing Age values. Returns estimates for the mean, sd, and lambda of the different age group classifications once optimisation is complete. Optimises to a tolerance of epsiolon, or convergence halts after maxit iterations.
}
\usage{teamEM(data, epsilon, maxit)}
\arguments{
  \item{data}{Dataframe with columns for Length and Age, including unknown values.}
  \item{epsilon}{Convergence threshold for when the change in log likelihoods is < epsilon}
  \item{maxit}{Maximum iterations allowed before algorithm halts optimisation}
}
\value{
  \item{estimates}{Dataframe of final estimates for mu, sigma, and lambda for Age groups 1, 2, and 3}
  \item{inits}{Dataframe of initialisation values for mu, sigma, and lambda for Age groups 1, 2, and 3}
  \item{converged}{TRUE/FALSE did the algorithm converge on values given the convergence threshold epsilon}
  \item{posteriors}{Dataframe of fnial posterior probabilities of each observed fish length belonging to its assigned Age group}
  \item{likelihood}{Vector of likelihood values for each iteration in the algorithm}
}
\examples{
data <- read_csv("data/docExampleData.csv")
teamEM(data)
}

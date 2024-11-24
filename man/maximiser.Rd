\name{maximiser}
\alias{maximiser}
\title{Maximisation Step for EM Algorithm}
\description{
Inputs the entire data and posterior probabilities. Returns a dataframe of updated estimate of mu, sigma, and lambda for each age group.
}
\usage{maximiser(data, posterior)}
\arguments{
  \item{sorted_data}{Dataframe with columns for Length and Age, with all observations with known Age at the top, and all unknown Age obervations at the bottom.}
  \item{posterior}{Dataframe with posterior probabilities of the length belonging to each age group}
}
\value{
  \item{mu}{Updated mu estimates}
  \item{sigma}{Updated sigma estimates}
  \item{lambda}{Updated lambda estimates}
}
\examples{
data <- read_csv("data.csv")
posteriors <- data.frame(matrix(c(1, 0, 0, 1, 1, 1, 0, 1, 0),
                        nrow = 3,
                        ncol = 3))
maximiser(data, posteriors)
}

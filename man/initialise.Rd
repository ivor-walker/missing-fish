\name{initialise}
\alias{initialise}
\title{Initialise Data for EM Algorithm}
\description{
Returns a data frame containing initial guesses for the mean, sd, and lambda for each age group. Estimates are based on the known values given in the data.
}
\usage{initialise(data)}
\arguments{
  \item{data}{Dataframe with columns for Length and Age, including unknown values.}
}
\value{
  \item{mu}{mu estimates for Age group 1, 2, 3}
  \item{sd}{sd estimates for Age group 1, 2, 3}
  \item{lambda}{lambda estimates for Age group 1, 2, 3}
}
\examples{
data <- read_csv("data.csv")
initialise(data)
}


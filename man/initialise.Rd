\name{initialise}
\alias{initialise}
\title{Initialise Data for EM Algorithm}
\description{
Returns a data frame containing initial guesses for the mean, sd, and lambda for each age group. Estimates are based on the known values given in the data.
}
\usage{initialise(data)}
\arguments{
  \item{known}{Dataframe with columns for Length and Age, only for known Age}
  \item{unknown}{Dataframe with columns for Length and Age, only for NA Age}
  \item{sorted_data}{Dataframe with columns for Length and Age, with all observations with known Age at the top, and all unknown Age obervations at the bottom.}
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


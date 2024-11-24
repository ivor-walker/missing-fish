\name{expector}
\alias{expector}
\title{Expectation Step for EM Algorithm}
\description{
Inputs the entire data and the data estimates. Returns the posterior probabilities of each length belonging to each age class, as well as the Gaussian probability densities of each length based on the current estimated parameteres.
}
\usage{expector(data, estimates)}
\arguments{
  \item{known}{Dataframe with columns for Length and Age, only for known Age}
  \item{sorted_data}{Dataframe with columns for Length and Age, with all observations with known Age at the top, and all unknown Age obervations at the bottom.}
  \item{estimates}{3x3 matrix / dataframe with row 1 containing mu estimates, row 2 containing sd   estimates, and row 3 containg lambda estimates all for age group j}
}
\value{
  \item{posteriors}{Dataframe containing posterior probabilities that a length belongs to each age   group, for each length in the data}
  \item{densities}{Gaussian probability densities of each length based on the current estimated     parameters}
}
\examples{
data <- read_csv("data.csv")
estimates <- data.frame(matrix(c(5, 10, 15, 5, 6, 7, .5, .6, .7),
                        nrow = 3,
                        ncol = 3))
colnames(estimates) <- c("mu", "sigma", "lambda")
rownames(estimates) <- c("Age1", "Age2", "Age3")

expector(data, estimates)
}


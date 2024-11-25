\name{generateDataset}
\alias{generateDataset}
\title{Generate Random Dataset for Testing}
\description{Generates a dataset similar to FishLengths with a specified proportion of NA values.}
\usage{generateDataset(seed, fish, datasize, init, distributions)}
\arguments{
  \item{seed}{Set the random seed.}
  \item{fish}{Specify the amount of fish.}
  \item{datasize}{Specify the amount of length observations. The size of the dataset.}
  \item{init}{Initialisation parameters from FishLengths dataset.}
  \item{distributions}{A length-4 vector with probabilities for the occurance of (NA, Age1, Age2, and Age3), respectively.}
}
\value{A simulated dataset with columns FishID, Length, and Age}
\examples{
data <- read_csv("data/docExampleData.csv")

known <- data[!is.na(data$Age), ] # known data is not NA
unknown <- data[is.na(data$Age), ] # unknown data is NA
sorted_data <- rbind(known, unknown) # combine data
init <- initialise(knwon, unknown, sorted_data)


generateDataset(0, 100, 1000, init, c(0.9, 0.04, 0.03, 0.03))
}

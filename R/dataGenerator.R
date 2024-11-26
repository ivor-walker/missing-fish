library(dplyr)

#' Generate Random Dataset for Testing
#'
#' @param seed Set the random seed.
#' @param fish Specify the amount of fish.
#' @param datasize Specify the amount of length observations. The size of the dataset.
#' @param inits Initialisation parameters from FishLengths dataset.
#' @param distributions A length-4 vector with probabilities for the occurance of (NA, Age1, Age2, and Age3), respectively.
#'
#' @return A simulated dataset with columns FishID, Length, and Age
#' @export
#'
#' @examples data <- readr::read_csv("data/docExampleData.csv")
#' known <- data[!is.na(data$Age), ] # known data is not NA
#' unknown <- data[is.na(data$Age), ] # unknown data is NA
#' sorted_data <- rbind(known, unknown) # combine data
#' inits <- initialise(known, unknown, sorted_data)
#'
#' generateDataset(0, 100, 1000, inits, c(0.9, 0.04, 0.03, 0.03))
generateDataset <- function(seed, fish, datasize, inits, distributions) {
  set.seed(seed)

  FishID <- 1:fish # there are `fish` many fish in the experiment

  # sample different values for Age with replacement with probability specified in `distributions`
  Age <- sample(c(NA, 1, 2, 3), size = datasize, replace = TRUE, prob = distributions)

  lengths <- numeric(datasize) # initialize lengths column

  for (i in 1:datasize) {
    if (is.na(Age[i])) { # if Age is NA for that observation
      # take a random value for length with a random mean and sd from the intialisation of FishLengths Dataset
      lengths[i] <- rnorm(1, mean = sample(inits$mu, 1), sd = sample(inits$sigma, 1))
      # else take a random value for length in accordance with the mean and sd from that given Age group in the initialisation of FishLengths Dataset
    } else {
      lengths[i] <- rnorm(1, mean = inits$mu[Age[i]], sd = inits$sigma[Age[i]])
    }
  }

  dataset <- data.frame(FishID = FishID, Length = lengths, Age = Age)
  dataset <- dataset %>% arrange(FishID)

  return(dataset)
}

#data1 <- generateDataset(0, 100, 1000, inits, c(0.9, 0.04, 0.03, 0.03))

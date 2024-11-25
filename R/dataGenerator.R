library(dplyr)

generateDataset <- function(seed, fish, datasize, init, distributions) {
  set.seed(seed)

  FishID <- 1:fish # there are `fish` many fish in the experiment

  # sample different values for Age with replacement with probability specified in `distributions`
  Age <- sample(c(NA, 1, 2, 3), size = datasize, replace = TRUE, prob = distributions)

  lengths <- numeric(datasize) # initialize lengths colummn

  for (i in 1:datasize) {
    if (is.na(Age[i])) { # if Age is NA for that observation
      # take a random value for length with a random mean and sd from the intialisation of FishLengths Dataset
      lengths[i] <- rnorm(1, mean = sample(init$mu, 1), sd = sample(init$sigma, 1))
      # else take a random value for length in accordane with the mean and sd from that given Age group in the initialisation of FishLengths Dataset
    } else {
      lengths[i] <- rnorm(1, mean = init$mu[Age[i]], sd = init$sigma[Age[i]])
    }
  }

  dataset <- data.frame(FishID = FishID, Length = lengths, Age = Age)
  dataset <- dataset %>% arrange(FishID)

  return(dataset)
}

#data1 <- generateDataset(0, 100, 1000, init, c(0.9, 0.04, 0.03, 0.03))

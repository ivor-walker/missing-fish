library(dplyr)

generateDataset <- function(seed, fish, datasize, inits, distributions) {
  set.seed(seed)
  
  FishID <- 1:fish
  
  Age <- sample(c(NA, 1, 2, 3), size = datasize, replace = TRUE, prob = distributions)
  
  lengths <- numeric(datasize)
  
  for (i in 1:datasize) {
    if (is.na(Age[i])) {
      lengths[i] <- rnorm(1, mean = sample(inits$mu, 1), sd = sample(inits$sigma, 1))
    } else {
      lengths[i] <- rnorm(1, mean = inits$mu[Age[i]], sd = inits$sigma[Age[i]])
    }
  }
  
  dataset <- data.frame(FishID = FishID, Length = lengths, Age = Age)
  dataset <- dataset %>% arrange(FishID)
  
  return(dataset)
}

data1 <- generateDataset(0, 100, 1000, inits, c(0.9, 0.04, 0.03, 0.03))
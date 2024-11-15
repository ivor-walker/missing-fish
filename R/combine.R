library(tidyverse)

#Initialisation - Narayan
initialise <- function(data) {
  initialisation <- list()
  # extract the mean lengths in each age group from the frogs with known ages
  meanLengths <- data %>%
    group_by(Age) %>%
    summarise(meanLength = mean(Length)) %>%
    select(meanLength) %>%
    unlist()
  
  # using these means, generate bounds to initialise the age groups with
  group1Cutoff <- (meanLengths[1] + meanLengths[2])/2
  group2Cutoff <- (meanLengths[2] + meanLengths[3])/2
  
  # mutate the Age column based on these initial bounds
  data <- data %>%
    mutate(Age = if_else(!is.na(Age), Age, # if the entry already has an value, don't mutate it
                         # if the entry is currently NA, change it in accordance with my bounds
                         if_else(0<=Length & Length<=group1Cutoff, 1, # if Length is within group1Cutoff, place it in age group 1
                                 if_else(group1Cutoff<Length & Length<=group2Cutoff, 2, # if length is within group1 and group2Cutoff, place it in age group 2
                                         if_else(group2Cutoff<Length, 3, Age))))) # if length is greater than group2Cutoff, place it in age group 3
  
  # summarise the initialisation
  initDataSummary <- data %>%
    group_by(Age) %>%
    summarise(meanLength = mean(Length),
              sdLength = sd(Length))
  
  
  # extract the initial parameter values from this summary
  muHat <- pull(initDataSummary[,2]) # pull the second column from initialise
  sigmaHat <- pull(initDataSummary[,3]) # pull the third row from initialise
  lambdaHat <- as.vector(table(data$Age) / nrow(data)) # compute the lambdas as proportions
  
  return(list(
    muHat = muHat,
    sigmaHat = sigmaHat,
    lambdaHat = lambdaHat
  ))
}

#Expectation - Lee
expector <- function(data, initialisation) {
  # Creating subsets of data set for groups with known/ unknown ages
  known <- data[!is.na(data$Age), ]
  unknown <- data[is.na(data$Age), ]
  
  # Combining known and unknown data (after filled) for full data set
  combined_data <- rbind(known, unknown)
  
  # Initialising information to be used throughout function
  age_groups <- sort(unique(known$Age))
  k <- length(age_groups)
  rows <- nrow(data)
  
  # Creating empty data frame to be filled with densities
  densities <- data.frame(matrix(0, nrow = nrow(data), ncol = k))
  colnames(densities) <- c("Age1", "Age2", "Age3")
  
  # Using initial guess as first estimate
  estimates <- initialisation
  
  # Calculating densities for each yi for each age group with Gaussian probability density function
  # based off estimates for mu and sigma for unknown age groups
  for (i in ((length(known$Age) + 1):rows)) {
    
    yi <- data$Length[i]
    
    densities[i, ] <- c(dnorm(yi, mean = estimates$mu[1], sd = estimates$sigma[1]),
                        dnorm(yi, mean = estimates$mu[2], sd = estimates$sigma[2]),
                        dnorm(yi, mean = estimates$mu[3], sd = estimates$sigma[3]))
  }
  
  # Creating empty data frame to store posteriors in
  posteriors <- data.frame(matrix(0, nrow = nrow(data), ncol = k))
  colnames(posteriors) <- c("Age1", "Age2", "Age3")
  
  # Calculating posteriors for each yi for each age group for unknown age groups
  for (i in ((length(known$Age) + 1):rows)) {
    
    yi <- combined_data$Length[i]
    densities_yi <- densities[i, ] # Getting densities for specific yi
    
    # Calculating total probability
    Pyi <-0
    for (j in 1:k) {
      Pyi <- Pyi + densities_yi[j]*estimates$lambda[j]
    }
    
    # Calculating and storing posteriors using Bayes Theorem
    for (l in 1:k) {
      posteriors[i, l] <- (densities_yi[l]*estimates$lambda[l])/Pyi
    }
  }
  
  for (i in 1:(length(known$Age))) {
    posteriors[i, known$Age[i]] <- 1
  }
  
  return(list(
    posteriors = posteriors,
    densities = densities
  ))
}

#Checking convergence - Ivor
converger <- function(data, expectations, epsilon, maxit) {
  converged <- FALSE
  logLikelihoods <- numeric(maxit)
  iterations <- 0
  optimised <- NULL
  
  while (!converged) {
    iterations <- iterations + 1
    optimised <- optimise(data, expectations$posteriors)
    expectations <- expector(data, optimised)
    logLikelihoods[iterations] <- findLogLikelihood(data, expectations$densities, optimised)
    converged <- checkConvergence(logLikelihoods, iterations, maxit, epsilon)
  }
  
  return(list(
    converged = converged,
    logLikelihoods = logLikelihoods,
    iterations = iterations,
    optimised = optimised
  ))
}

#Checking convergence of a specific point
checkConvergence <- function(logLikelihoods, iterations, maxit, epsilon) {
  #Must immediately break
  minIterations <- 2
  if(iterations < minIterations) {
    return(FALSE)
  }
  
  hitMaxIterations <- iterations >= maxit
  
  #If change is too small for tolerance, function has converged
  change <- abs(logLikelihoods[iterations] - logLikelihoods[iterations - 1])
  changeBelowTolerance <- change < epsilon
  
  return(hitMaxIterations || changeBelowTolerance)
}

#Compute log likelihood
findLogLikelihood <- function(data, densities, optimised) {
  logLikelihood <- 0
  n <- nrow(data)
  k <- length(optimised$mu)
  
  for(i in 1:n){
    likelihood <- sum(densities[i, ] * optimised$lambda)
    if (likelihood > 0) {
      logLikelihood <- logLikelihood + log(likelihood)
    } else {
      #To avoid log(0) errors, TODO
      logLikelihood <- logLikelihood + log(.Machine$double.eps)
    }
  }
  return(logLikelihood)
}

#Optimisation - other group member
optimise <- function(data, posterior) {
  #TODO
  
  return(list(
    mew <- mew,
    sigma <- sigma,
    lambda <- lambda
  ))  
}

#Main function
expectationMaximisation <- function(data, epsilon = 1e-08, maxit = 1000) {
  initialisation <- initialise(data)
  expectations <- expector(data, initialisation)
  convergence <- converger(data, expectations, epsilon, maxit)
  result <- arrangeResult(initialisation, convergence)
  return(result)
}

#Main method
data <- load("FishLengths.RData")
result <- expectationMaximisation(data)


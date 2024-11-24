library(tidyverse)

##### Initialisation - Narayan & Lee #####
initialise <- function(known, unknown, sorted_data) {
  age_groups <- sort(unique(known$Age))
  k <- length(age_groups) # amount of age groups
  rows <- nrow(sorted_data)
  urows <- length(unknown$Age) # amount of unknown ages

  mu_known <- tapply(known$Length, known$Age, mean) # mean length of each age in the known data

  # for each unknown row, assign the age group which has its observed mu closest to that length
  for (i in 1:urows) {
    unknown$Age[i] <- age_groups[which.min(abs(unknown$Length[i] - mu_known))]
  }

  mu <- tapply(sorted_data$Length, sorted_data$Age, mean) # mean length of each age in the combined data
  sigma <- tapply(sorted_data$Length, sorted_data$Age, sd) # sd for each age in the combined data
  lambda <- as.numeric(table(sorted_data$Age)/nrow(sorted_data)) # lambda for each age in the combined data

  inits <- data.frame(mu, sigma, lambda) # combine these initial estimates into a df
  rownames(inits) <- c("Age1", "Age2", "Age3")

  return (inits)
}

##### Expectation - Lee #####
expector <- function(known, sorted_data, estimates) {

  age_groups <- sort(unique(known$Age))
  k <- length(age_groups) # amount of age groups
  rows <- nrow(sorted_data)

  densities <- data.frame(matrix(0, nrow = rows, ncol = k)) # initialise densities object
  colnames(densities) <- c("Age1", "Age2", "Age3")

  for (i in 1:rows) {
    yi <- sorted_data$Length[i] # observed length i

    densities[i, ] <- c(dnorm(yi, mean = estimates$mu[1], sd = estimates$sigma[1]), # Gaussian pdf with mu and sd of Age 1
                        dnorm(yi, mean = estimates$mu[2], sd = estimates$sigma[2]), # Gaussian pdf with mu and sd of Age 2
                        dnorm(yi, mean = estimates$mu[3], sd = estimates$sigma[3])) # Gaussian pdf with mu and sd of Age 3
  }

  posteriors <- data.frame(matrix(0, nrow = nrow(sorted_data), ncol = k)) # initialise posteriors object

  colnames(posteriors) <- c("Age1", "Age2", "Age3")

  for (i in ((length(known$Age) + 1):rows)) { # iterating through the unknown data only
    yi <- sorted_data$Length[i] # observed length
    densities_yi <- densities[i, ] # observed Gaussian pdf

    # sum across each age group the product of the density and the lambda value
    Pyi <-0
    for (j in 1:k) {
      Pyi <- Pyi + densities_yi[j]*estimates$lambda[j]
    }

    # calculate the posterior probility of the length belonging to each age group, entering into the posteriors object accordingly
    for (l in 1:k) {
      posteriors[i, l] <- (densities_yi[l]*estimates$lambda[l])/Pyi
    }
  }

  for (i in 1:(length(known$Age))) { # iterating through the known values
    posteriors[i, known$Age[i]] <- 1 # the posterior probability of the known values is 1 since we already know it
  }

  return(list(
    posteriors = posteriors,
    densities = densities
  ))
}

##### Maximisation - Yi #####
maximiser <- function(sorted_data, posteriors) {
  k <- ncol(posteriors)
  N <- nrow(sorted_data)

  # Initialise lists to store updated parameters
  mu <- c()
  sigma <- c()
  lambda <- c()

  for (j in 1:k) {
    P_ij <- posteriors[, j] # posterior probabilities for age group j

    mu <- c(mu, sum(P_ij * sorted_data$Length) / sum(P_ij)) # compute updated mu estimate
    sigma <- c(sigma, sqrt(sum(P_ij * (sorted_data$Length - mu[j])^2) / sum(P_ij))) # compute updated sd estimate
    lambda <- c(lambda, sum(P_ij) / N) # compute updated lambda estimate
  }

  estimates <- data.frame(mu, sigma, lambda) # update estimates object with parameters above
  rownames(estimates) <- c("Age1", "Age2", "Age3")
  colnames(estimates) <- c("mu", "sigma", "lambda")
  return(estimates)
}

##### Convergence - Ivor #####
teamEM <- function(unsorted_data, epsilon = 1e-08, maxit = 1000) {
  known <- unsorted_data[!is.na(unsorted_data$Age), ] # known data is not NA
  unknown <- unsorted_data[is.na(unsorted_data$Age), ] # unknown data is NA
  sorted_data <- rbind(known, unknown) # combine data

  # initialise all variables for the loop
  converged <- FALSE
  iterations <- 0
  minIterations <- 2 # minimum 2 iterations so that the change criteria runs

  expectations <- NULL

  # initialise estimates to data initialisation
  inits <- initialise(known, unknown, sorted_data)
  estimates <- inits

  logLikelihoods <- numeric(maxit)
  change <- 0

  startTime <- Sys.time()

  while (!converged && iterations < maxit) { # while not converged and within max iterations
    iterations <- iterations + 1

    expectations <- expector(known, sorted_data, estimates) # complete expectation step
    expectationsTime <- Sys.time()
    expectationsTimeTaken <- expectationsTime - startTime

    estimates <- maximiser(sorted_data, expectations$posteriors) # complete maximisation step
    maximiserTime <- Sys.time()
    maximiserTimeTaken <- maximiserTime - expectationsTime

    logLikelihoods[iterations] <- findLogLikelihood(sorted_data, expectations$densities, estimates) # compute loglikelihood for this maximisation
    logLikelihoodTime <- Sys.time()
    likelihoodTimeTaken <- logLikelihoodTime - maximiserTime

    change <- abs(logLikelihoods[iterations] - logLikelihoods[iterations - 1])# check if function converged based on previous loglikelihood value
    converged <- change < epsilon && iterations > minIterations

    initTime <- Sys.time()
    changeTime <- initTime - startTime
    startTime <- initTime
    print(estimates)
    print(paste("iteration:", iterations, " | delta(logLikelihood):", round(change, 9), "| time for iteration to complete:", round(changeTime, 3), "s | expectations:", round(expectationsTimeTaken, 3), "s | maximiser: ", round(maximiserTimeTaken, 5),"s | logLikelihood:", round(likelihoodTimeTaken, 3), "s"))
  }

  logLikelihoods <- head(logLikelihoods, iterations)
  return(list(
    estimates = estimates,
    inits = inits,
    converged = converged,
    posteriors = expectations$posteriors,
    logLikelihoods = logLikelihoods
  ))
}

findLogLikelihood <- function(data, densities, estimates) {
  loglikelihood <- 0

  N <- nrow(data)
  K <- length(estimates$mu) # amount of age groups

  for(i in 1:N){ # iterate through the rows of the data
    likelihood <- 0

    for(k in 1:K) { # iterate through each age group
      likelihood <- likelihood + estimates$lambda[k] * densities[i, k] # add the product of the lambda estimate and density for that age group
    }

    loglikelihood <- loglikelihood + log(likelihood) # add the log of this value to the overall likelihood, since the log of a product is a sum of logs.
  }

  return(loglikelihood)
}


load("data/FishLengths.RData")
result <- teamEM(x)

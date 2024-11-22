library(tidyverse)

##### Initialisation - Narayan & Lee #####
initialise <- function(data) {
  known <- data[!is.na(data$Age), ] # known values are not NA values
  unknown <- data[is.na(data$Age), ] # unknown values are NA values

  age_groups <- sort(unique(known$Age))
  k <- length(age_groups) # amount of age groups
  rows <- nrow(data)
  urows <- length(unknown$Age) # amount of unknown ages

  mu_known <- tapply(known$Length, known$Age, mean) # mean length of each age in the known data

  # for each unknown row, assign the age group which has its observed mu closest to that length
  for (i in 1:urows) {
    unknown$Age[i] <- age_groups[which.min(abs(unknown$Length[i] - mu_known))]
  }

  combined_data <- rbind(known, unknown) # combine data

  mu <- tapply(combined_data$Length, combined_data$Age, mean) # mean length of each age in the combined data
  sigma <- tapply(combined_data$Length, combined_data$Age, sd) # sd for each age in the combined data
  lambda <- as.numeric(table(combined_data$Age)/nrow(combined_data)) # lambda for each age in the combined data

  inits <- data.frame(mu, sigma, lambda) # combine these initial estimates into a df
  rownames(inits) <- c("Age1", "Age2", "Age3")

  return (inits)
}

##### Expectation - Lee #####
expector <- function(data, estimates) {
  known <- data[!is.na(data$Age), ] # known values are not NA values
  unknown <- data[is.na(data$Age), ] # unknown values are NA values

  combined_data <- rbind(known, unknown) #combine data

  age_groups <- sort(unique(known$Age))
  k <- length(age_groups) # amount of age groups
  rows <- nrow(data)

  densities <- data.frame(matrix(0, nrow = rows, ncol = k)) # initialise densities object
  colnames(densities) <- c("Age1", "Age2", "Age3")

  for (i in 1:rows) {

    yi <- data$Length[i] # observed length i

    densities[i, ] <- c(dnorm(yi, mean = estimates$mu[1], sd = estimates$sigma[1]), # Gaussian pdf with mu and sd of Age 1
                        dnorm(yi, mean = estimates$mu[2], sd = estimates$sigma[2]), # Gaussian pdf with mu and sd of Age 2
                        dnorm(yi, mean = estimates$mu[3], sd = estimates$sigma[3])) # Gaussian pdf with mu and sd of Age 3
  }

  posteriors <- data.frame(matrix(0, nrow = nrow(data), ncol = k)) # initialise posteriors object
  colnames(posteriors) <- c("Age1", "Age2", "Age3")

  for (i in ((length(known$Age) + 1):rows)) { # iterating through the unknown data only

    yi <- combined_data$Length[i] # observed length
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
maximiser <- function(data, posteriors) {
  known <- data[!is.na(data$Age), ] # known values are not NA values
  unknown <- data[is.na(data$Age), ] # unknown values are NA values

  data <- rbind(known, unknown) # combine data

  k <- ncol(posteriors)
  N <- nrow(data)

  mu <- c()
  sigma <- c()
  lambda <- c()

  for (j in 1:k) {
    P_ij <- posteriors[, j] # posterior probabilities for age group j

    mu <- c(mu, sum(P_ij * data$Length) / sum(P_ij)) # compute updated mu estimate as defined in algorithm
    sigma <- c(sigma, sqrt(sum(P_ij * (data$Length - mu[j])^2) / sum(P_ij))) # compute sd estimate as defined in algorithm
    lambda <- c(lambda, sum(P_ij) / N) # compute lambda estimate as defined in algorithm
  }

  estimates <- data.frame(mu, sigma, lambda) # update estimates object with this maximisation
  rownames(estimates) <- c("Age1", "Age2", "Age3")

  return(estimates)
}

##### Convergence - Ivor #####
checkConvergence <- function(logLikelihoods, iterations, epsilon) {

  minIterations <- 2 # minimum 2 iterations so that the change criteria runs
  if(iterations < minIterations) {
    return(FALSE)
  }

  change <- abs(logLikelihoods[iterations] - logLikelihoods[iterations - 1]) # calculate the change in log likelihood from the last iteration
  #if (iterations %% 100 == 0) {
  #  print(change)
  #}
  changeBelowTolerance <- change < epsilon # check if the change is below tolerance

  return(changeBelowTolerance) # if returns TRUE, function has converged
}


findLogLikelihood <- function(data, densities, estimates) {
  known <- data[!is.na(data$Age), ] # known values are not NA
  unknown <- data[is.na(data$Age), ] # unknown values are NA

  data <- rbind(known, unknown) # combine data

  likelihood <- 0

  N <- nrow(data)
  K <- length(estimates$mu) # amount of age groups

  for(n in 1:N){ # iterate through the rows of the data
    a <- 0

    for(k in 1:K) { # iterate through each age group
      a <- a + estimates$lambda[k] * densities[n, k] # add the product of the lambda estimate and density for that age group
    }

    likelihood <- likelihood + log(a) # add the log of this value to the overall likelihood, since the log of a product is a sum of logs.
    }

  return(likelihood)
}

#
converger <- function(data, inits, epsilon, maxit) {
  known <- data[!is.na(data$Age), ] # known data is not NA
  unknown <- data[is.na(data$Age), ] # unknown data is NA

  data <- rbind(known, unknown) # combine data

  # initialise all variables for the loop
  iterations <- 0
  estimates <- inits # initialise estimates to data initialisation
  logLikelihoods <- numeric(maxit)
  converged <- FALSE

  while (!converged && iterations < maxit) { # while not converged and within max iterations
    iterations <- iterations + 1
    expectations <- expector(data, estimates) # complete expectation step
    estimates <- maximiser(data, expectations$posteriors) # complete maximisation step
    logLikelihoods[iterations] <- findLogLikelihood(data, expectations$densities, estimates) # compute loglikelihood for this maximisation
    converged <- checkConvergence(logLikelihoods, iterations, epsilon) # check if function converged based on previous loglikelihood value
  }

  return(list(
    iterations = iterations,
    estimates = estimates,
    expectations = expectations,
    logLikelihoods = logLikelihoods,
    converged = converged
  ))
}


##### Run Algorithm #####
teamEM <- function(data, epsilon = 1e-08, maxit = 1000) {
  inits <- initialise(data)
  convergence <- converger(data, inits, epsilon, maxit)
  return (list(convergence$estimates,
               inits,
               convergence$converged,
               convergence$expectations$posteriors,
               convergence$logLikelihoods))
}

load("FishLengths.RData")
result <- teamEM(x)

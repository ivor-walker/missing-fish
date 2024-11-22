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

#Expectation - Lee
expector <- function(data, estimates) {
  known <- data[!is.na(data$Age), ]
  unknown <- data[is.na(data$Age), ]

  combined_data <- rbind(known, unknown)

  age_groups <- sort(unique(known$Age))
  k <- length(age_groups)
  rows <- nrow(data)

  densities <- data.frame(matrix(0, nrow = rows, ncol = k))
  colnames(densities) <- c("Age1", "Age2", "Age3")

  for (i in 1:rows) {

    yi <- data$Length[i]

    densities[i, ] <- c(dnorm(yi, mean = estimates$mu[1], sd = estimates$sigma[1]),
                        dnorm(yi, mean = estimates$mu[2], sd = estimates$sigma[2]),
                        dnorm(yi, mean = estimates$mu[3], sd = estimates$sigma[3]))
  }

  posteriors <- data.frame(matrix(0, nrow = nrow(data), ncol = k))
  colnames(posteriors) <- c("Age1", "Age2", "Age3")

  for (i in ((length(known$Age) + 1):rows)) {

    yi <- combined_data$Length[i]
    densities_yi <- densities[i, ]

    Pyi <-0
    for (j in 1:k) {
      Pyi <- Pyi + densities_yi[j]*estimates$lambda[j]
    }

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

#
maximiser <- function(data, posteriors) {
  known <- data[!is.na(data$Age), ]
  unknown <- data[is.na(data$Age), ]

  data <- rbind(known, unknown)

  k <- ncol(posteriors)
  N <- nrow(data)

  mu <- c()
  sigma <- c()
  lambda <- c()

  for (j in 1:k) {
    P_ij <- posteriors[, j]

    mu <- c(mu, sum(P_ij * data$Length) / sum(P_ij))

    sigma <- c(sigma, sqrt(sum(P_ij * (data$Length - mu[j])^2) / sum(P_ij)))

    lambda <- c(lambda, sum(P_ij) / N)
  }

  estimates <- data.frame(mu, sigma, lambda)
  rownames(estimates) <- c("Age1", "Age2", "Age3")

  return(estimates)
}

#
checkConvergence <- function(logLikelihoods, iterations, epsilon) {

  minIterations <- 2
  if(iterations < minIterations) {
    return(FALSE)
  }

  change <- abs(logLikelihoods[iterations] - logLikelihoods[iterations - 1])
  if (iterations %% 100 == 0) {
    print(change)
  }
  changeBelowTolerance <- change < epsilon

  return(changeBelowTolerance)
}

#
findLogLikelihood <- function(data, densities, estimates) {
  known <- data[!is.na(data$Age), ]
  unknown <- data[is.na(data$Age), ]

  data <- rbind(known, unknown)

  likelihood <- 0

  N <- nrow(data)
  K <- length(estimates$mu)

  for(n in 1:N){
    a <- 0

    for(k in 1:K) {
      a <- a + estimates$lambda[k] * densities[n, k]
    }

    likelihood <- likelihood + log(a)
    }

  return(likelihood)
}

#
converger <- function(data, inits, epsilon, maxit) {
  known <- data[!is.na(data$Age), ]
  unknown <- data[is.na(data$Age), ]

  data <- rbind(known, unknown)

  iterations <- 0
  estimates <- inits
  logLikelihoods <- numeric(maxit)
  converged <- FALSE

  while (!converged && iterations < maxit) {
    iterations <- iterations + 1
    expectations <- expector(data, estimates)
    estimates <- maximiser(data, expectations$posteriors)
    logLikelihoods[iterations] <- findLogLikelihood(data, expectations$densities, estimates)
    converged <- checkConvergence(logLikelihoods, iterations, epsilon)
  }

  return(list(
    iterations = iterations,
    estimates = estimates,
    expectations = expectations,
    logLikelihoods = logLikelihoods,
    converged = converged
  ))
}

inits <- initialise(x)
convergence <- converger(x, inits, 1e-08, 1000)
results <- (list(convergence$estimates,
             inits,
             convergence$converged,
             convergence$expectations$posteriors,
             convergence$logLikelihoods))


#teamEM <- function(data, epsilon = 1e-08, maxit = 1000) {
#  inits <- initialise(data)
#  convergence <- converger(data, inits, epsilon, maxit)
#  return (list(convergence$estimates,
#               inits,
#               convergence$converged,
#               convergence$expectations$posteriors,
#               convergence$logLikelihoods))
#}
#
#result <- teamEM(x)

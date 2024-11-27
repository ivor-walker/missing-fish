##### Initialisation - Narayan & Lee #####
#' Initialise Data for EM Algorithm
#'
#' @param known Data frame with columns for Length and Age, only for known Age
#' @param unknown Data frame with columns for Length and Age, only for NA Age
#' @param sorted_data Data frame with columns for Length and Age, with all observations with known Age at the top, and all unknown Age obervations at the bottom.
#'
#' @return A data frame containing initial estimates for mu, sigma, and lambda for Age groups 1, 2, and 3
#' @export
#'
#' @examples
#' known <- docExampleData[!is.na(docExampleData$Age), ] # known data is not NA
#' unknown <- docExampleData[is.na(docExampleData$Age), ] # unknown data is NA
#' sorted_data <- rbind(known, unknown) # combine data
#'
#' initialise(known, unknown, sorted_data)
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

  sorted_data <- rbind(known, unknown) # combine data with newly assigned unknowns

  mu <- tapply(sorted_data$Length, sorted_data$Age, mean) # mean length of each age in the combined data
  sigma <- tapply(sorted_data$Length, sorted_data$Age, sd) # sd for each age in the combined data
  lambda <- as.numeric(table(sorted_data$Age)/nrow(sorted_data)) # lambda for each age in the combined data

  inits <- data.frame(mu, sigma, lambda) # combine these initial estimates into a df
  rownames(inits) <- c("Age1", "Age2", "Age3")

  return (inits)
}

##### Expectation - Lee #####
#' Expectation Step for EM Algorithm
#'
#' @param known Data frame with columns for Length and Age, only for known Age
#' @param sorted_data Data frame with columns for Length and Age, with all observations with known Age at the top, and all unknown Age obervations at the bottom.
#' @param estimates 3x3 matrix / dataframe with row 1 containing mu estimates, row 2 containing sd   estimates, and row 3 containg lambda estimates all for age group j
#'
#' @return A list containing posteriors and densities
#' \itemize{
#'  \item posteriors - Data frame containing posterior probabilities that a length belongs to each age group, for each length in the data
#'  \item densities - Gaussian probability densities of each length based on the current estimated parameters
#' }
#' @export
#'
#' @examples
#' estimates <- data.frame(matrix(c(5, 10, 15, 5, 6, 7, .5, .6, .7),
#'                               nrow = 3,
#'                               ncol = 3))
#' colnames(estimates) <- c("mu", "sigma", "lambda")
#' rownames(estimates) <- c("Age1", "Age2", "Age3")
#'
#' known <- docExampleData[!is.na(docExampleData$Age), ] # known data is not NA
#' unknown <- docExampleData[is.na(docExampleData$Age), ] # unknown data is NA
#' sorted_data <- rbind(known, unknown) # combine data
#'
#' expector(known, sorted_data, estimates)
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

    # calculate the posterior probability of the length belonging to each age group, entering into the posteriors object accordingly
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
#' Maximisation Step for EM Algorithm
#'
#' @param sorted_data Data frame with columns for Length and Age, with all observations with known Age at the top, and all unknown Age obervations at the bottom.
#' @param posteriors Data frame with posterior probabilities of the length belonging to each age group
#'
#' @return A data frame containing maximised estimates for mu, sigma, and lambda for Age groups 1, 2, and 3
#' @export
#'
#' @examples
#' posteriors <- data.frame(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),
#'                          nrow = 3,
#'                          ncol = 3))
#' maximiser(docExampleData, posteriors)
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
#' Perform EM Optimisation
#'
#' @param data Data frame with columns for Length and Age, including unknown values.
#' @param epsilon Convergence threshold for when the change in log likelihoods is < epsilon
#' @param maxit Maximum iterations allowed before algorithm halts optimisation
#'
#' @return A list of results from the optimisation:
#' \itemize{
#'  \item estimates - Data frame of final estimates for mu, sigma, and lambda for Age groups 1, 2, and 3
#'  \item inits - Data frame of initialisation values for mu, sigma, and lambda for Age groups 1, 2, and 3
#'  \item converged - TRUE/FALSE did the algorithm converge on values given the convergence threshold epsilon
#'  \item posteriors - Data frame of fnial posterior probabilities of each observed fish length belonging to its assigned Age group
#'  \item likelihood - Vector of likelihood values for each iteration in the algorithm
#' }
#' @export
#'
#' @examples
#' teamEM(docExampleData)
teamEM <- function(data, epsilon = 1e-08, maxit = 1000) {
  known <- data[!is.na(data$Age), ] # known data is not NA
  unknown <- data[is.na(data$Age), ] # unknown data is NA
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

  #startTime <- Sys.time()

  while (!converged && iterations < maxit) { # while not converged and within max iterations
    iterations <- iterations + 1

    expectations <- expector(known, sorted_data, estimates) # complete expectation step
    #expectationsTime <- Sys.time()
    #expectationsTimeTaken <- expectationsTime - startTime

    estimates <- maximiser(sorted_data, expectations$posteriors) # complete maximisation step
    #maximiserTime <- Sys.time()
    #maximiserTimeTaken <- maximiserTime - expectationsTime

    logLikelihoods[iterations] <- findLogLikelihood(sorted_data, expectations$densities, estimates) # compute loglikelihood for this maximisation
    #logLikelihoodTime <- Sys.time()
    #likelihoodTimeTaken <- logLikelihoodTime - maximiserTime

    change <- abs(logLikelihoods[iterations] - logLikelihoods[iterations - 1])# check if function converged based on previous loglikelihood value
    converged <- change < epsilon && iterations > minIterations

    #initTime <- Sys.time()
    #changeTime <- initTime - startTime
    #startTime <- initTime
    print(estimates)
    print(paste("iteration:", iterations, " | delta(logLikelihood):", round(change, 9)))
                #,"| time for iteration to complete:", round(changeTime, 3), "s | expectations:", round(expectationsTimeTaken, 3), "s | maximiser: ", round(maximiserTimeTaken, 5),"s | logLikelihood:", round(likelihoodTimeTaken, 3), "s"))
  }

  logLikelihoods <- head(logLikelihoods, iterations)
  return(list(
    estimates = estimates,
    inits = inits,
    converged = converged,
    posteriors = expectations$posteriors,
    likelihood = logLikelihoods
  ))
}

#' Find Log Likelihood
#'
#' @param data Data frame with columns for Length and Age, including unknown values.
#' @param densities Data frame containing Gaussian probability densities of observing each observed fish length under the estimated mu and sigma.
#' @param estimates Data frame containing the current estimates for mu, sigma, and lamba for Age group 1, 2, and 3
#'
#' @return Log likelihood of the current estiamtes for mu, sigma, and lambda belonging to the given Fish Length data.
#' @export
#'
#' @examples
#' known <- docExampleData[!is.na(docExampleData$Age), ] # known data is not NA
#' unknown <- docExampleData[is.na(docExampleData$Age), ] # unknown data is NA
#' sorted_data <- rbind(known, unknown) # combine data
#' init <- initialise(known, unknown, sorted_data)
#'
#' exp <- expector(known, sorted_data, init)
#' estimates <- maximiser(sorted_data, exp$posteriors)
#'
#' findLogLikelihood(docExampleData, exp$densities, estimates)
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

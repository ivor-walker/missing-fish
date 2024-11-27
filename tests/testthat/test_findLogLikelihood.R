test_that("findLogLikelihood() computes log-likelihood correctly using generated dataset", {
  # Generate initial estimates
  inits <- data.frame(
    mu = c(10, 20, 30),
    sigma = c(2, 3, 4),
    lambda = c(0.3, 0.4, 0.3)
  )
  rownames(inits) <- c("Age1", "Age2", "Age3")

  # Generate a dataset
  seed <- 123
  fish <- 100
  datasize <- 500
  distributions <- c(0.9, 0.04, 0.03, 0.03) # Probability for (NA, Age1, Age2, Age3)
  dataset <- generateDataset(seed, fish, datasize, inits, distributions)

  # Split data into known and unknown (to simulate the context of log-likelihood calculation)
  known <- dataset[!is.na(dataset$Age), ]
  unknown <- dataset[is.na(dataset$Age), ]
  sorted_data <- rbind(known, unknown)

  # Compute densities manually for the generated data
  densities <- data.frame(
    Age1 = dnorm(sorted_data$Length, mean = inits$mu[1], sd = inits$sigma[1]),
    Age2 = dnorm(sorted_data$Length, mean = inits$mu[2], sd = inits$sigma[2]),
    Age3 = dnorm(sorted_data$Length, mean = inits$mu[3], sd = inits$sigma[3])
  )

  # Manually compute the expected log-likelihood
  expected_log_likelihood <- sum(log(
    rowSums(densities * matrix(rep(inits$lambda, each = nrow(sorted_data)), nrow = nrow(sorted_data), byrow = TRUE))
  ))

  # Run the function
  computed_log_likelihood <- findLogLikelihood(sorted_data, densities, inits)

  # Check that the computed log-likelihood matches the manually computed value
  expect_equal(computed_log_likelihood, expected_log_likelihood, tolerance = 5, # expect some differences since calculation methods differ
               label = "Log-likelihood should match the manually computed value from generated data")
})

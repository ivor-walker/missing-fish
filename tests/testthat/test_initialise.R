test_that("initialise() correctly computes initial estimates", {
  # Sample input data
  known <- data.frame(
    Length = c(10, 20, 22, 30, 15, 25),
    Age = c(1, 2, 2, 3, 1, 3)
  )

  unknown <- data.frame(
    Length = c(12, 28),
    Age = NA
  )

  sorted_data <- rbind(known, unknown)

  # Expected values (manually computed or based on the logic of the initialise function)
  expected_mu <- c(mean(c(10, 15, 12)), mean(c(20, 22)), mean(c(30, 25, 28)))
  expected_sigma <- c(sd(c(10, 15, 12)), sd(c(20, 22)), sd(c(30, 25, 28)))
  expected_lambda <- c(3 / 8, 2 / 8, 3 / 8)

  # Run the function
  result <- initialise(known, unknown, sorted_data)

  # Check the outputs
  expect_equal(result$mu, expected_mu, tolerance = 1e-2, label = "Mu values should match")
  expect_equal(result$sigma, expected_sigma, tolerance = 1e-2, label = "Sigma values should match")
  expect_equal(result$lambda, expected_lambda, tolerance = 1e-2, label = "Lambda values should match")
})

test_that("maximiser() correctly returns an update of estimates",
          {
            # Sample input data
            sorted_data <- data.frame(
              Length = c(10, 20, 30, 12, 28, 22, 19, 16, 32, 9),
              Age = c(1, 2, 3, NA, NA, NA, NA, NA, NA, NA)
            )

            # Posterior probabilities with only 1s (deterministic assignments)
            posteriors <- data.frame(
              Age1 = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 1),
              Age2 = c(0, 1, 0, 0, 0, 1, 1, 1, 0, 0),
              Age3 = c(0, 0, 1, 0, 1, 0, 0, 0, 1, 0)
            )

            # Expected result (pre-calculated)
            expected_mu <- c(10.33333, 19.25, 30)
            expected_sigma <- c(1.247219, 2.165064, 1.632993)
            expected_lambda <- c(0.3, 0.4, 0.3)

            # Run the function
            result <- maximiser(sorted_data, posteriors)

            # Test result
            expect_equal(result$mu, expected_mu, tolerance = 1e-05, label = "mu values should match")
            expect_equal(result$sigma, expected_sigma, tolerance = 1e-05, label = "sigma values should match")
            expect_equal(result$lambda, expected_lambda, tolerance = 1e-05, label = "lambda values should match")
          })

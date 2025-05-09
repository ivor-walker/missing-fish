test_that("teamEM() is able to converge using synthetic data",
          {
            # Get estimates from the original data
            file_path <- testthat::test_path("data/x.rda")
            load(file_path)
            known <- x[!is.na(x$Age),]
            unknown <- x[is.na(x$Age),]
            sorted_data <- rbind(known, unknown)
            init <- initialise(known, unknown, sorted_data)

            # Generate datasets based on estimates above
            # Scenario 1: larger data size
            gen_data1 <- generateDataset(235, 150, 1500, init, c(0.9, 0.02, 0.04, 0.04))
            result1 <- teamEM(gen_data1)
            expect_true(result1$converged, label = "1: teamEM() returns a converged value")

            # Scenario 2: different proportion of known values
            gen_data2 <- generateDataset(263, 500, 1000, init, c(0.5, 0.1, 0.2, 0.2))
            result2 <- teamEM(gen_data2)
            expect_true(result2$converged, label = "2: teamEM() returns a converged value")

            # Scenario 3: different proportions of age groups
            gen_data3 <- generateDataset(235, 100, 1000, init, c(0.9, 0.01, 0.06, 0.03))
            result3 <- teamEM(gen_data3)
            expect_true(result3$converged, label = "3: teamEM() returns a converged value")

            # Scenario 4: a worse initial guess
            worse_init <- data.frame(matrix(c(100, 100, 100, 20, 20, 20, 0.1, 0.1, 0.8),
                                            nrow = 3,
                                            ncol = 3))
            colnames(worse_init) <- c("mu", "sigma", "lambda")
            gen_data4 <- generateDataset(276, 100, 1000, worse_init, c(0.9, 0.02, 0.04, 0.04))
            result4 <- teamEM(gen_data4)
            expect_true(result4$converged, label = "4: teamEM() returns a converged value")
          })

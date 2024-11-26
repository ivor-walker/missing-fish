test_that("teamEM() is able to converge using synthetic data",
          {
            # Get estimates from the original data
            file_path <- testthat::test_path("../../data/x.rda")
            load(file_path)
            known <- x[!is.na(x$Age),]
            unknown <- x[is.na(x$Age),]
            sorted_data <- rbind(known, unknown)
            init <- initialise(known, unknown, sorted_data)

            # Generate datasets based on estimates above
            gen_data <- generateDataset(235, 200, 2000, init, c(0.9, 0.02, 0.04, 0.04))
            result <- teamEM(gen_data)
            expect_true(result$converged, label = "teamEM() returns a converged value")
          })

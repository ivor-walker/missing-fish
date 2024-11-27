test_that("expector() correctly returns desired values",
          {
            # Get estimates from the original data
            file_path <- testthat::test_path("data/x.rda")
            load(file_path)
            known <- x[!is.na(x$Age),]
            unknown <- x[is.na(x$Age),]
            sorted_data <- rbind(known, unknown)
            init <- initialise(known, unknown, sorted_data)

            # Generate datasets based on estimates above
            gen_data <- generateDataset(12, 200, 2000, init, c(0.9, 0.02, 0.04, 0.04))
            gen_known <- gen_data[!is.na(gen_data$Age),]
            gen_unknown <- gen_data[is.na(gen_data$Age),]
            gen_sorted <- rbind(gen_known, gen_unknown)
            estimates <- initialise(gen_known, gen_unknown, gen_sorted)

            # Test expector function
            result <- expector(gen_known, gen_sorted, estimates)
            row_sums <- rowSums(result$posteriors)
            expect_equal(row_sums, rep(1, 2000), tolerance = 1e-5, label = "Row sums of posteriors should equal 1")
            expect_true(all(result$densities > 0), label = "All densities should be greater than zero")
          })


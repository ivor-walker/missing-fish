test_that("maximiser() correctly returns an update of estimates",
          {
            # Get estimates from the original data
            file_path <- testthat::test_path("../../data/x.rda")
            load(file_path)
            known <- x[!is.na(x$Age),]
            unknown <- x[is.na(x$Age),]
            sorted_data <- rbind(known, unknown)
            init <- initialise(known, unknown, sorted_data)

            # Generate datasets based on estimates above
            gen_data <- generateDataset(142, 200, 2000, init, c(0.9, 0.02, 0.04, 0.04))
            gen_known <- gen_data[!is.na(gen_data$Age),]
            gen_unknown <- gen_data[is.na(gen_data$Age),]
            gen_sorted <- rbind(gen_known, gen_unknown)
            estimates <- initialise(gen_known, gen_unknown, gen_sorted)

            # Get posterior probabilities and densities
            expect <- expector(gen_known, gen_sorted, estimates)

            # Test maximiser() function
            result1 <- maximiser(sorted_data, expect$posteriors)
          })

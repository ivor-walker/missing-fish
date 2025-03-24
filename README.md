# EMGroup1

EMGroup1 performs the Expectation-Maximisation (EM) Algorithm on Fish
Lengths Data that has age groups assigned to only a portion of the data.
The package can initialise the missing age data to estimates based on
the given data. It can also calculate the posterior probability that
each observed fish length belongs to each age group based on current
guesses for the parameters of each age group, then optimise these
parameter guesses based on those posterior probabilities. Finally, the
package can allow the parameter estimates to converge on values that
they are most likely to take. The function `teamEM` performs all steps
of the optimisation in one call, returning the optimised mu, sigma, and
labmda parameters for Age groups 1, 2, and 3.

## Function Overview

Here are examples of how to use the most important functions in the
package, namely `initialise`, `expector`, `maximise`, and `teamEM`.

#### `initialise`

Inputs: known - Data frame with columns for Length and Age, only for
known Age. unknown - Data frame with columns for Length and Age, only
for NA Age. sorted_data - Data frame with columns for Length and Age,
with all observations with known Age at the top, and all unknown Age
observations at the bottom.

Returns: A data frame containing initial estimates for mu, sigma, and
lambda for Age groups 1, 2, and 3.

``` r
library(EMGroup1)

# sort data between known and unknown age groups
known <- docExampleData[!is.na(docExampleData$Age), ] # known data is not NA
unknown <- docExampleData[is.na(docExampleData$Age), ] # unknown data is NA
sorted_data <- rbind(known, unknown) # combine data

initialise(known, unknown, sorted_data)
#>            mu    sigma    lambda
#> Age1 21.33333 3.511885 0.4285714
#> Age2 40.50000 9.192388 0.2857143
#> Age3 71.00000 5.656854 0.2857143
```

#### `expector`

Inputs: known - Data frame with columns for Length and Age, only for
known Age. sorted_data - Data frame with columns for Length and Age,
with all observations with known Age at the top, and all unknown Age
observations at the bottom. estimates - 3x3 matrix / data frame with row
1 containing mu estimates, row 2 containing sd estimates, and row 3
containing lambda estimates all for age group j.

Returns: A list containing posteriors and densities.

``` r
# create sample estimates
estimates <- data.frame(matrix(c(5, 10, 15, 5, 6, 7, .5, .6, .7),
                              nrow = 3,
                              ncol = 3))
colnames(estimates) <- c("mu", "sigma", "lambda")
rownames(estimates) <- c("Age1", "Age2", "Age3")

known <- docExampleData[!is.na(docExampleData$Age), ] # known data is not NA
unknown <- docExampleData[is.na(docExampleData$Age), ] # unknown data is NA
sorted_data <- rbind(known, unknown) # combine data

expector(known, sorted_data, estimates)
#> $posteriors
#>           Age1       Age2      Age3
#> 1 1.000000e+00 0.00000000 0.0000000
#> 2 0.000000e+00 1.00000000 0.0000000
#> 3 0.000000e+00 0.00000000 1.0000000
#> 4 0.000000e+00 0.00000000 1.0000000
#> 5 8.288755e-04 0.10856126 0.8906099
#> 6 6.753974e-03 0.21051891 0.7827271
#> 7 1.946283e-06 0.01317296 0.9868251
#> 
#> $densities
#>           Age1         Age2         Age3
#> 1 2.716594e-03 2.733501e-02 5.199096e-02
#> 2 3.802163e-17 3.673938e-10 1.651497e-06
#> 3 3.261221e-35 1.679656e-21 5.926927e-14
#> 4 2.193213e-44 2.178298e-27 6.340700e-18
#> 5 2.676605e-05 2.921383e-03 2.054255e-02
#> 6 4.768176e-04 1.238519e-02 3.947074e-02
#> 7 3.954639e-09 2.230504e-05 1.432231e-03
```

#### `maximise`

Inputs: sorted_data - Data frame with columns for Length and Age, with
all observations with known Age at the top, and all unknown Age
observations at the bottom. posteriors - Data frame with posterior
probabilities of the length belonging to each age group.

Returns: A data frame containing maximised estimates for mu, sigma, and
lambda for Age groups 1, 2, and 3.

``` r
# create sample posteriors 
posteriors <- data.frame(matrix(c(1, 0, 0, 1, 1, 1, 0, 1, 0),
                         nrow = 3,
                         ncol = 3))
maximiser(docExampleData, posteriors)
#>             mu     sigma    lambda
#> Age1 134.00000 159.25451 0.1428571
#> Age2  95.66667  89.48205 0.4285714
#> Age3  68.00000  51.47815 0.1428571
```

#### `teamEm`

Inputs: data - Data frame with columns for Length and Age, including
unknown values. epsilon - Convergence threshold for when the change in
log likelihoods is \< epsilon. Default of 1e-08. maxit - Maximum
iterations allowed before algorithm halts optimisation. Default of 1000.

Returns: A list of results from the optimisation: estimates - Data frame
of final estimates for mu, sigma, and lambda for Age groups 1, 2, and 3
inits - Data frame of initialisation values for mu, sigma, and lambda
for Age groups 1, 2, and 3 converged - TRUE/FALSE did the algorithm
converge on values given the convergence threshold epsilon posteriors -
Data frame of fnial posterior probabilities of each observed fish length
belonging to its assigned Age group likelihood - Vector of likelihood
values for each iteration in the algorithm

``` r
# run the optimisation on the example data
teamEM(docExampleData)
#>            mu    sigma    lambda
#> Age1 21.24744 2.917237 0.4121977
#> Age2 39.57833 7.387296 0.3020880
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 1  | delta(logLikelihood): "
#>            mu    sigma    lambda
#> Age1 21.22712 2.847222 0.4145558
#> Age2 39.75066 7.196402 0.2997299
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 2  | delta(logLikelihood): 0.375908167"
#>            mu    sigma    lambda
#> Age1 21.23754 2.847993 0.4161410
#> Age2 39.83460 7.121917 0.2981447
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 3  | delta(logLikelihood): 0.002586459"
#>            mu    sigma    lambda
#> Age1 21.24398 2.849448 0.4169902
#> Age2 39.87868 7.083652 0.2972955
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 4  | delta(logLikelihood): 8.2953e-05"
#>            mu    sigma    lambda
#> Age1 21.24733 2.850226 0.4174276
#> Age2 39.90143 7.063785 0.2968581
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 5  | delta(logLikelihood): 0.000167605"
#>            mu    sigma    lambda
#> Age1 21.24904 2.850625 0.4176507
#> Age2 39.91306 7.053591 0.2966350
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 6  | delta(logLikelihood): 0.000124563"
#>            mu    sigma    lambda
#> Age1 21.24991 2.850827 0.4177640
#> Age2 39.91897 7.048397 0.2965217
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 7  | delta(logLikelihood): 7.408e-05"
#>            mu   sigma    lambda
#> Age1 21.25035 2.85093 0.4178215
#> Age2 39.92196 7.04576 0.2964642
#> Age3 71.00000 4.00000 0.2857143
#> [1] "iteration: 8  | delta(logLikelihood): 4.0417e-05"
#>            mu    sigma    lambda
#> Age1 21.25057 2.850982 0.4178506
#> Age2 39.92348 7.044424 0.2964352
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 9  | delta(logLikelihood): 2.121e-05"
#>            mu    sigma    lambda
#> Age1 21.25068 2.851008 0.4178653
#> Age2 39.92425 7.043748 0.2964204
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 10  | delta(logLikelihood): 1.0925e-05"
#>            mu    sigma    lambda
#> Age1 21.25074 2.851022 0.4178727
#> Age2 39.92464 7.043406 0.2964130
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 11  | delta(logLikelihood): 5.576e-06"
#>            mu    sigma    lambda
#> Age1 21.25077 2.851028 0.4178765
#> Age2 39.92483 7.043233 0.2964092
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 12  | delta(logLikelihood): 2.833e-06"
#>            mu    sigma    lambda
#> Age1 21.25078 2.851032 0.4178784
#> Age2 39.92493 7.043146 0.2964073
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 13  | delta(logLikelihood): 1.436e-06"
#>            mu    sigma    lambda
#> Age1 21.25079 2.851034 0.4178793
#> Age2 39.92498 7.043101 0.2964064
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 14  | delta(logLikelihood): 7.27e-07"
#>            mu    sigma    lambda
#> Age1 21.25079 2.851034 0.4178798
#> Age2 39.92501 7.043079 0.2964059
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 15  | delta(logLikelihood): 3.68e-07"
#>            mu    sigma    lambda
#> Age1 21.25080 2.851035 0.4178801
#> Age2 39.92502 7.043068 0.2964056
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 16  | delta(logLikelihood): 1.86e-07"
#>            mu    sigma    lambda
#> Age1 21.25080 2.851035 0.4178802
#> Age2 39.92503 7.043062 0.2964055
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 17  | delta(logLikelihood): 9.4e-08"
#>            mu    sigma    lambda
#> Age1 21.25080 2.851035 0.4178803
#> Age2 39.92503 7.043059 0.2964054
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 18  | delta(logLikelihood): 4.8e-08"
#>            mu    sigma    lambda
#> Age1 21.25080 2.851035 0.4178803
#> Age2 39.92503 7.043058 0.2964054
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 19  | delta(logLikelihood): 2.4e-08"
#>            mu    sigma    lambda
#> Age1 21.25080 2.851035 0.4178803
#> Age2 39.92503 7.043057 0.2964054
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 20  | delta(logLikelihood): 1.2e-08"
#>            mu    sigma    lambda
#> Age1 21.25080 2.851035 0.4178803
#> Age2 39.92503 7.043056 0.2964054
#> Age3 71.00000 4.000000 0.2857143
#> [1] "iteration: 21  | delta(logLikelihood): 6e-09"
#> $estimates
#>            mu    sigma    lambda
#> Age1 21.25080 2.851035 0.4178803
#> Age2 39.92503 7.043056 0.2964054
#> Age3 71.00000 4.000000 0.2857143
#> 
#> $inits
#>            mu    sigma    lambda
#> Age1 21.33333 3.511885 0.4285714
#> Age2 40.50000 9.192388 0.2857143
#> Age3 71.00000 5.656854 0.2857143
#> 
#> $converged
#> [1] TRUE
#> 
#> $posteriors
#>           Age1        Age2         Age3
#> 1 1.0000000000 0.000000000 0.000000e+00
#> 2 0.0000000000 1.000000000 0.000000e+00
#> 3 0.0000000000 0.000000000 1.000000e+00
#> 4 0.0000000000 0.000000000 1.000000e+00
#> 5 0.9326729759 0.067327024 2.067017e-29
#> 6 0.9922637547 0.007736245 5.713125e-35
#> 7 0.0002255546 0.999774445 6.363074e-19
#> 
#> $likelihood
#>  [1] -27.47987 -27.10396 -27.10137 -27.10146 -27.10162 -27.10175 -27.10182
#>  [8] -27.10186 -27.10188 -27.10190 -27.10190 -27.10190 -27.10191 -27.10191
#> [15] -27.10191 -27.10191 -27.10191 -27.10191 -27.10191 -27.10191 -27.10191
```

## Installation

You can install the development version of EMGroup1 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MT41132024/assignment-2-em-algorithm-emgroup1")
```

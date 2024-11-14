library(dplyr)

# Gaining understanding of data set
head(x)
summary(x)
length(unique(x$FishID))
unique(x$Age)

data <- x
data <- data %>% arrange(FishID)
maxit <- 1000

#teamEM <- function(data, epsilon = 1e-08, maxit = 1000) {

# --------------------------------------- Step 1: Initialisation ---------------------------------------

# Creating subsets of data set for groups with known/ unknown ages
known <- data[!is.na(data$Age), ]
unknown <- data[is.na(data$Age), ]

# Initialising information to be used throughout function
age_groups <- sort(unique(known$Age))
k <- length(age_groups)
rows <- nrow(data)

# Calculating mean mus and sigmas for each age in known subset 
mu_known <- tapply(known$Length, known$Age, mean)
sigma_known <- tapply(known$Length, known$Age, sd) # Probably don't need this?

# Using length in unknown age data set alongside average length found for age group in known data set
# to assign age groups in unknown based on where the has smallest difference in length is observed
for (i in 1:nrow(unknown)) {
  unknown$Age[i] <- age_groups[which.min(abs(unknown$Length[i] - mu_known))]
}

# Combining known and unknown data (after filled) for full data set
combined_data <- rbind(known, unknown)
  
# Computing initial values for mu, sigma and lambda
mu <- tapply(combined_data$Length, combined_data$Age, mean)
sigma <- tapply(combined_data$Length, combined_data$Age, sd)
lambda <- as.numeric(table(combined_data$Age)/nrow(combined_data))
  
# Defining initial values data frame
inits <- data.frame(mu, sigma, lambda)
rownames(inits) <- c("Age1", "Age2", "Age3")

estimates <- inits # Using initial guess as first estimate

#for (it in 1:maxit) {
  # --------------------------------------- Step 2: Expectation ---------------------------------------
  
  # Creating empty data frame to be filled with densities
  densities <- data.frame(matrix(0, nrow = nrow(data), ncol = k))
  colnames(densities) <- c("Age1", "Age2", "Age3")
  
  # Calculating densities for each yi for each age group with Gaussian probability density function
  # based off estimates for mu and sigma for unknown age groups
  for (i in ((length(known$Age) + 1):rows)) {
    
    yi <- data$Length[i]
  
    densities[i, ] <- c(dnorm(yi, mean = estimates$mu[1], sd = estimates$sigma[1]),
                        dnorm(yi, mean = estimates$mu[2], sd = estimates$sigma[2]),
                        dnorm(yi, mean = estimates$mu[3], sd = estimates$sigma[3]))
  }
  
  # Creating empty data frame to store posteriors in
  posteriors <- data.frame(matrix(0, nrow = nrow(data), ncol = k))
  colnames(posteriors) <- c("Age1", "Age2", "Age3")
  
  # Calculating posteriors for each yi for each age group for unknown age groups
  for (i in ((length(known$Age) + 1):rows)) {
    
    yi <- combined_data$Length[i]
    densities_yi <- densities[i, ] # Getting densities for specific yi
    
    # Calculating total probability
    Pyi <-0
    for (j in 1:k) {
      Pyi <- Pyi + densities_yi[j]*estimates$lambda[j]
    }
    
    # Calculating and storing posteriors using Bayes Theorem
    for (l in 1:k) {
      posteriors[i, l] <- (densities_yi[l]*estimates$lambda[l])/Pyi
    }
  }
  
  for (i in 1:(length(known$Age))) {
    posteriors[i, known$Age[i]] <- 1
  }
  
  # --------------------------------------- Step 3: Maximisation ---------------------------------------

  for (j in 1:k) {
    
    upper_mu <- 0
    lower <- 0
    
    for (i in 1:rows) {
      upper_mu <- upper_mu + posteriors[i, j] * combined_data$Length[i]
      lower <- lower + posteriors[i, j]
    }
    estimates$mu[j] <- upper_mu/lower
    
    upper_sigma <- 0
    for (i in 1:rows) {
      upper_sigma <- upper_sigma + (posteriors[i, j] * (combined_data$Length[i]-estimates$mu[j])^2)
    }
    estimates$sigma[j] <- sqrt(upper_sigma/lower)
    
    estimates$lambda[j] <- (1/rows)*lower
  }
  
  print(estimates)
#}

# issue with the fact that we are calculating posteriors for lengths which we already know age for
# need to find a way to only go over unkown values while still looking at whole dataset
# maybe just do it over unkown dataset instead, solved issue by ordering dataset the combinining ordered
# known and unknown, then looking over only range of unkown and setting posterior for known to 1 by default

#}

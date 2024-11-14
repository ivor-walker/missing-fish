library(tidyverse)
load("FishLengths.RData")


# extract the mean lengths in each age group from the frogs with known ages
meanLengths <- x %>%
  group_by(Age) %>%
  summarise(meanLength = mean(Length)) %>%
  select(meanLength) %>%
  unlist()

# using these means, generate bounds to initialise the age groups with
group1Cutoff <- (meanLengths[1] + meanLengths[2])/2
group2Cutoff <- (meanLengths[2] + meanLengths[3])/2


# mutate the Age column based on these initial bounds
x <- x %>%
  mutate(Age = if_else(!is.na(Age), Age, # if the entry already has an value, don't mutate it
                       # if the entry is currently NA, change it in accordance with my bounds
                       if_else(0<=Length & Length<=group1Cutoff, 1, # if Length is within group1Cutoff, place it in age group 1
                       if_else(group1Cutoff<Length & Length<=group2Cutoff, 2, # if length is within group1 and group2Cutoff, place it in age group 2
                               if_else(group2Cutoff<Length, 3, Age))))) # if length is greater than group2Cutoff, place it in age group 3

# summarise the initialisation
initialise <- x %>%
  group_by(Age) %>%
  summarise(meanLength = mean(Length),
            sdLength = sd(Length))


# extract the initial parameter values from this summary
muHat <- pull(initialise[,2]) # pull the second column from initialise
sigmaHat <- pull(initialise[,3]) # pull the third row from initialise
lambdaHat <- as.vector(table(x$Age) / nrow(x)) # compute the lambdas as proportions





library(dinemites)

dataset_in <- data.frame(allele = c('A', NA, NA, NA, 'B'),
     subject = rep('A', 5),
     time = c(1, 29, 30, 39, 50))
treatments <- data.frame(subject = c('A'),
                         time = c(29))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_present_infection(dataset)
dataset <- add_persistent_column(dataset)
dataset <- add_lag_column(dataset)
dataset <- add_treatment_column(dataset, treatments)

expect_that(dataset$treatment_acute, equals(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0)))
expect_that(dataset$treatment_longitudinal, equals(c(0, 0, 0, 0, 0, 0, 1, 0, 1, 0)))

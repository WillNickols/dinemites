library(dinemites)

dataset_in <- data.frame(allele = c('A', NA, NA, NA, 'B', 'A'),
                         subject = rep('A', 6),
                         time = c(1, 29, 30, 39, 50, 70))
treatments <- data.frame(subject = c('A'),
                         time = c(29))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_present_infection(dataset)
dataset <- add_persistent_column(dataset)
dataset <- add_persistent_infection(dataset)
dataset <- add_lag_column(dataset)
dataset <- add_lag_infection(dataset)
dataset <- add_treatment_column(dataset, treatments)
dataset <- add_treatment_infection(dataset, treatments)

expect_that(dataset$treatment_acute, equals(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)))
expect_that(dataset$treatment_longitudinal, equals(c(0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0)))
expect_that(dataset$treatment_acute_infection,
            equals(c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0)))
expect_that(dataset$treatment_longitudinal_infection,
            equals(c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0)))


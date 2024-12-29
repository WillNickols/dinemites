library(dinemites)

dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
    subject = rep('A', 5),
    time = c(1, 1, 14, 28, 42))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_present_infection(dataset)
dataset <- add_lag_column(dataset)

expect_that(dataset$lag_30, equals(c(0, 0, 1, 1, 1, 1, 0, 0)))

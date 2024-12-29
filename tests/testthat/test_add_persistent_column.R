library(dinemites)

dataset_in <- data.frame(allele = c('A', 'B', NA, NA),
    subject = rep('A', 4),
    time = c(1, 2, 3, 4))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_present_infection(dataset)
dataset <- add_persistent_column(dataset)

expect_that(dataset$persistent, equals(c(0, 0, 1, 0, 1, 1, 1, 1)))

library(dinemites)

dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
    subject = rep('A', 5),
    time = c(1, 1, 2, 3, 4))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_present_infection(dataset)
dataset <- add_persistent_infection(dataset)

expect_that(dataset$persistent_infection, equals(c(0, 0, 1, 1, 1, 1, 1, 1)))

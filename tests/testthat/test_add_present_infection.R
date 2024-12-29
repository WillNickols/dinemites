library(dinemites)

dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
    subject = rep('A', 5),
    time = c(1, 1, 2, 3, 4))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_present_infection(dataset)

expect_that(dataset$present_infection, equals(c(1, 1, 0, 0, 0, 0, 0, 0)))

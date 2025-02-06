library(dinemites)

dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
    subject = rep('A', 5),
    time = c(1, 1, 3, 8, 17))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_time_gap(dataset, -14)

expect_that(dataset$time_gap, equals(c(15, 15, 2, 2, 5, 5, 9, 9)))

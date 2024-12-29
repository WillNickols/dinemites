library(dinemites)

dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
subject = rep('A', 5),
time = c(1, 1, 2, 3, 4))

qpcr_times <- data.frame(subject = rep('A', 2),
                         time = c(1, 4))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_qpcr_times(dataset, qpcr_times)

expect_that(dataset$present, equals(c(1, 1, 0, 0, 0, 0, 2, 2)))

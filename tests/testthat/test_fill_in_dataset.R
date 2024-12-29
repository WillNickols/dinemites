library(dinemites)

dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
    subject = rep('A', 5),
    time = c(1, 1, 2, 3, 4))

dataset <- fill_in_dataset(dataset_in)

expect_that(dataset$allele, equals(rep(c("A", "B"), 4)))
expect_that(dataset$subject, equals(rep("A", 8)))
expect_that(dataset$time, equals(rep(c(1,2,3,4), each = 2)))
expect_that(dataset$present, equals(c(1, 1, 0, 0, 0, 0, 0, 0)))

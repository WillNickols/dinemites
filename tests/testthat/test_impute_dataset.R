library(dinemites)

set.seed(1)
dataset_in <- data.frame(allele = c('A', 'A', 'A', NA, NA, 'B', NA, 'B'),
    subject = rep('A', 8),
    time = c(1, 2, 3, 4, 5, 6, 7, 8))

qpcr_times <- data.frame(subject = rep('A', 1), time = c(7))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_qpcr_times(dataset, qpcr_times)
imputed_mat <- impute_dataset(dataset)
dataset$present_probability <- rowMeans(imputed_mat)

expect_that(dataset$present_probability,
            equals(c(1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1)))

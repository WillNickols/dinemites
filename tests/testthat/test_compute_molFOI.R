library(dinemites)
library(dplyr)

dataset_in <- data.frame(
    locus = rep(c(1, 2), each = 5),
    allele = c('A', 'B', 'B', NA, NA, 'A', 'A', 'B', NA, NA),
    subject = rep('A', 10),
    time = rep(c(1, 1, 5, 15, 44), 2)) %>%
    mutate(allele = interaction(allele, locus))

dataset <- fill_in_dataset(dataset_in)

dataset$probability_new <-
    determine_probabilities_simple(dataset)$probability_new

expect_that(unname(unlist(compute_molFOI(dataset, method = 'sum_then_max')[,'molFOI'])),
            equals(2))

expect_that(unname(unlist(compute_molFOI(dataset, method = 'max_then_sum')[,'molFOI'])),
            equals(3))

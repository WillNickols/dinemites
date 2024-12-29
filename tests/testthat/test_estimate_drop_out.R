library(dinemites)

dataset_in <- data.frame(allele = c('A', 'B', NA, 'A', 'A'),
                         subject = rep('A', 5),
                         time = c(1, 1, 14, 28, 42))

dataset <- fill_in_dataset(dataset_in)

named_results <- c(0.5, 2, 1)
names(named_results) <- c("Drop out rate", "Evaluated gaps", "Number present in gap")
expect_that(estimate_drop_out(dataset), equals(named_results))

library(dinemites)

dataset <- data.frame(allele = rep(c('A', 'B', 'C', 'D', 'E'), 5),
    subject = rep('A', 25),
    time = rep(1:5, each = 5),
    present = c(0,0,0,0,0, 1,0,0,1,0, 1,0,0,0,0, 1,0,0,0,0, 0,0,0,0,0))

dataset$probability_new <-
    determine_probabilities_simple(dataset)$probability_new

expect_that(dataset$probability_new,
            equals(c(rep(NA, 5), 1,NA,NA,1,NA, 0,NA,NA,NA,NA, 0,NA,NA,NA,NA, rep(NA, 5))))

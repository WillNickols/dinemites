% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/determine_probabilities.R
\name{add_probability_new}
\alias{add_probability_new}
\title{add_probability_new}
\usage{
add_probability_new(dataset, probability_mat, column_name = "probability_new")
}
\arguments{
\item{dataset}{A complete longitudinal dataset with columns \code{allele},
\code{subject}, \code{time}, and \code{present}}

\item{probability_mat}{A matrix with one column per result of
\verb{determine_probabilities_*} on an imputed dataset.}

\item{column_name}{What to name the new column}
}
\value{
The dataset with a new column given by \code{column_name}
giving the probability the allele was new if present.
}
\description{
Add a column with the probability an allele is new by aggregating a
probability matrix.
}
\examples{
library(foreach)
library(doParallel)
dataset_in <- data.frame(allele = c('A', 'A', 'A', NA, NA, 'B', NA, 'B'),
    subject = rep('A', 8),
    time = c(1, 2, 3, 4, 5, 6, 7, 8))

qpcr_times <- data.frame(subject = rep('A', 1), time = c(7))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_qpcr_times(dataset, qpcr_times)

n_imputations <- 10
imputed_mat <- impute_dataset(dataset, n_imputations)

probabilities_simple_mat <-
     foreach(i = 1:n_imputations,
         .combine = cbind,
         .packages = c('dinemites', 'dplyr')) \%do\% {
         dataset_tmp <- dataset
         dataset_tmp$present <- imputed_mat[,i]
         probabilities_simple <- determine_probabilities_simple(dataset_tmp)
         probabilities_simple$probability_new
     }

dataset <- add_probability_new(dataset, probabilities_simple_mat)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_dataset.R
\name{add_probability_present}
\alias{add_probability_present}
\title{add_probability_present}
\usage{
add_probability_present(dataset, imputation_mat)
}
\arguments{
\item{dataset}{A complete longitudinal dataset with columns \code{allele},
\code{subject}, \code{time}, and \code{present}}

\item{imputation_mat}{The output of the function \code{impute_dataset}}
}
\value{
The dataset with a new column \code{probability_present}
giving the proportion of imputations in which the allele was present.
}
\description{
Add a column with the probability an allele was present based on the
imputations.
}
\examples{
dataset_in <- data.frame(allele = c('A', 'A', 'A', NA, NA, 'B', NA, 'B'),
    subject = rep('A', 8),
    time = c(1, 2, 3, 4, 5, 6, 7, 8))

qpcr_times <- data.frame(subject = rep('A', 1), time = c(7))

dataset <- fill_in_dataset(dataset_in)
dataset <- add_qpcr_times(dataset, qpcr_times)
imputed_mat <- impute_dataset(dataset)
dataset <- add_probability_present(dataset, imputed_mat)

}

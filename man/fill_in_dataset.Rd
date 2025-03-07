% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_modeling.R
\name{fill_in_dataset}
\alias{fill_in_dataset}
\title{fill_in_dataset}
\usage{
fill_in_dataset(dataset)
}
\arguments{
\item{dataset}{A longitudinal dataset with columns \code{subject}, \code{time},
\code{allele.}
One row should be present for each allele sequenced at the subject-time.
All other subject-times should be included with allele as \code{NA}. A \code{locus}
column can also be specified to give each allele's locus.}
}
\value{
A data.frame with columns \code{subject}, \code{time}, \code{allele}, and \code{present}
with \code{present} as 1 if the allele is present and 0 otherwise.
}
\description{
Fill in zeros for a dataset that contains only alleles that were present.
}
\examples{

dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
    subject = rep('A', 5),
    time = c(1, 1, 2, 3, 4))

dataset <- fill_in_dataset(dataset_in)

}

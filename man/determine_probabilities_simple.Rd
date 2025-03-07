% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/determine_probabilities.R
\name{determine_probabilities_simple}
\alias{determine_probabilities_simple}
\title{determine_probabilities_simple}
\usage{
determine_probabilities_simple(dataset, n_lags = 3, t_lag = Inf)
}
\arguments{
\item{dataset}{A complete longitudinal dataset with columns \code{allele},
\code{subject}, \code{time}, and \code{present.}}

\item{n_lags}{If the allele has been observed in both the last
\code{n_lags} samples and the last \code{t_lag} times, it is considered new.}

\item{t_lag}{If the allele has been observed in both the last
\code{n_lags} samples and the last \code{t_lag} times, it is considered new.}
}
\value{
A named list with (1) \code{probability_new} as a 0/1 vector of whether
the corresponding allele in the dataset was from a new infection and (2)
\code{fit} as \code{NULL} (for consistency with the other models).
}
\description{
Determine whether an allele is new or persistent using a simple rule.
For each subject, an allele is counted as new if it has not been observed in
the \code{n_lags} most recent samples (e.g., visits) or in the \code{t_lag} most recent
times (e.g., days). If the allele has been observed in both the last
\code{n_lags} samples and the last \code{t_lag} times, it is considered new.
}
\examples{

dataset <- data.frame(allele = rep(c('A', 'B', 'C', 'D', 'E'), 5),
    subject = rep('A', 25),
    time = rep(1:5, each = 5),
    present = c(0,0,0,0,0, 1,0,0,1,0, 1,0,0,0,0, 1,0,0,0,0, 0,0,0,0,0))

dataset$probability_new <-
    determine_probabilities_simple(dataset)$probability_new

}

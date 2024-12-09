% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/determine_probabilities.R
\name{determine_probabilities_bayesian}
\alias{determine_probabilities_bayesian}
\title{determine_probabilities_bayesian}
\usage{
determine_probabilities_bayesian(
  dataset,
  infection_persistence_covariates = NULL,
  infection_general_covariates = NULL,
  alleles_persistence_covariates = NULL,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 100,
  adapt_delta = 0.99,
  seed = 1,
  drop_out = FALSE
)
}
\arguments{
\item{dataset}{dataset}

\item{infection_persistence_covariates}{infection_persistence_covariates}

\item{infection_general_covariates}{infection_general_covariates}

\item{alleles_persistence_covariates}{alleles_persistence_covariates}

\item{chains}{chains}

\item{parallel_chains}{parallel_chains}

\item{iter_warmup}{iter_warmup}

\item{iter_sampling}{iter_sampling}

\item{refresh}{refresh}

\item{adapt_delta}{adapt_delta}

\item{seed}{seed}

\item{drop_out}{drop_out}
}
\value{
An list with probabilities
}
\description{
determine_probabilities_bayesian
}

library(testthat)
library(dinemites)
library(dplyr)
library(doParallel)
library(foreach)
library(doRNG)

set.seed(1)

##############################################################
# Run all of this to ensure there are no errors              #
# Seeds will vary across machines, so only check simple rule #
##############################################################

dataset <- read.csv(system.file(package="dinemites","extdata","dataset.tsv"), header = TRUE, sep="\t")
treatments <- read.csv(system.file(package="dinemites","extdata","treatments.tsv"), header = TRUE, sep="\t")
qPCR_only <- read.csv(system.file(package="dinemites","extdata","qPCR_only.tsv"), header = TRUE, sep="\t")

dataset <- fill_in_dataset(dataset)
dataset <- add_qpcr_times(dataset, qpcr_times = qPCR_only)

estimate_drop_out(dataset)

n_imputations <- 2
imputed_datasets <- impute_dataset(dataset, n_imputations = n_imputations)
dataset <- add_probability_present(dataset, imputed_datasets)
expect_that(rowMeans(imputed_datasets), equals(dataset$probability_present))

plot_dataset(dataset, treatments)

probabilities_simple_mat <-
    matrix(ncol = n_imputations, nrow = nrow(dataset))
probabilities_bayesian_mat <-
    matrix(ncol = n_imputations, nrow = nrow(dataset))
probabilities_clustering_mat <-
    matrix(ncol = n_imputations, nrow = nrow(dataset))

n_cores <- ifelse(detectCores() > 1, 2, 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

probabilities_simple_mat <-
    foreach(i = 1:n_imputations, .combine = cbind, .packages = c('dinemites', 'dplyr')) %dorng% {
    dataset_tmp <- dataset
    dataset_tmp$present <- imputed_datasets[,i]

    probabilities_simple <- determine_probabilities_simple(dataset_tmp)
    probabilities_simple$probability_new
}

if (instantiate::stan_cmdstan_exists()) {
    probabilities_bayesian_mat <-
        foreach(i = 1:n_imputations, .combine = cbind, .packages = c('dinemites', 'dplyr')) %dorng% {
        dataset_tmp <- dataset
        dataset_tmp$present <- imputed_datasets[,i]

        dataset_tmp <- dataset_tmp %>%
            add_present_infection() %>%
            add_persistent_column() %>%
            add_persistent_infection() %>%
            add_lag_column() %>%
            add_lag_infection() %>%
            add_treatment_column(treatments = treatments, verbose = F) %>%
            add_treatment_infection(treatments = treatments, verbose = F)

        probabilities_bayesian <- determine_probabilities_bayesian(
            dataset = dataset_tmp,
            infection_persistence_covariates =
                c("persistent_infection",
                  "lag_infection_30",
                  "treatment_acute_infection",
                  "treatment_longitudinal_infection"),
            infection_general_covariates =
                c("season"),
            alleles_persistence_covariates =
                c("persistent",
                  "lag_30",
                  "treatment_acute",
                  "treatment_longitudinal"))
        probabilities_bayesian$probability_new
    }
}

if (instantiate::stan_cmdstan_exists()) {
    probabilities_clustering_mat <-
        foreach(i = 1:n_imputations, .combine = cbind, .packages = c('dinemites', 'dplyr')) %dorng% {
        dataset_tmp <- dataset
        dataset_tmp$present <- imputed_datasets[,i]

        probabilities_clustering <- determine_probabilities_clustering(dataset = dataset_tmp)
    }
}

stopCluster(cl)

dataset <- add_probability_new(dataset, probabilities_simple_mat)
expect_that(rowMeans(probabilities_simple_mat, na.rm = T),
            equals(dataset$probability_new))

compute_total_new_COI(dataset, method = 'sum_then_max')
compute_total_new_COI(dataset, method = 'max_then_sum')

estimated_new_infections <-
    estimate_new_infections(dataset,
                            imputation_mat = imputed_datasets,
                            probability_mat = probabilities_simple_mat)

plot_dataset(dataset, treatments, estimated_new_infections)







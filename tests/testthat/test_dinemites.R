library(testthat)
library(dinemites)
library(dplyr)

dataset <- read.csv(system.file(package="dinemites","extdata","dataset.tsv"), header = TRUE, sep="\t")
treatments <- read.csv(system.file(package="dinemites","extdata","treatments.tsv"), header = TRUE, sep="\t")
qPCR_only <- read.csv(system.file(package="dinemites","extdata","qPCR_only.tsv"), header = TRUE, sep="\t")

dataset <- fill_in_dataset(dataset)
dataset <- add_qpcr_times(dataset, qpcr_times = qPCR_only)

estimate_drop_out(dataset)

dataset <- dataset %>%
    dplyr::arrange(subject, allele, time)

n_imputations <- 5
imputed_datasets <- impute_dataset(dataset, n_imputations = n_imputations)

probabilities_simple_mat <-
    matrix(ncol = n_imputations, nrow = nrow(dataset))
probabilities_bayesian_mat <-
    matrix(ncol = n_imputations, nrow = nrow(dataset))
probabilities_clustering_mat <-
    matrix(ncol = n_imputations, nrow = nrow(dataset))
for (i in 1:n_imputations) {
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

    probabilities_simple <- determine_probabilities_simple(dataset_tmp)
    probabilities_simple_mat[,i] <- probabilities_simple$probability_new

    if (instantiate::stan_cmdstan_exists()) {
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
                  "treatment_longitudinal"),
            refresh = 100)
        probabilities_bayesian_mat[,i] <- probabilities_bayesian$probability_new

        probabilities_clustering <- determine_probabilities_clustering(dataset = dataset_tmp)
        probabilities_clustering_mat[,i] <- probabilities_clustering$probability_new
    }
}

probabilities_simple_mat_check <-
    read.table(system.file(package="dinemites", "extdata", "probabilities_simple_mat.tsv"),
               sep = '\t') %>% as.matrix()
probabilities_bayesian_mat_check <-
    read.table(system.file(package="dinemites", "extdata", "probabilities_bayesian_mat.tsv"),
               sep = '\t') %>% as.matrix()
probabilities_clustering_mat_check <-
    read.table(system.file(package="dinemites", "extdata", "probabilities_clustering_mat.tsv"),
               sep = '\t') %>% as.matrix()

if (!instantiate::stan_cmdstan_exists()) {
    probabilities_bayesian_mat <-
        read.table(system.file(package="dinemites", "extdata", "probabilities_bayesian_mat.tsv"),
                   sep = '\t') %>% as.matrix()
    probabilities_clustering_mat <-
        read.table(system.file(package="dinemites", "extdata", "probabilities_clustering_mat.tsv"),
                   sep = '\t') %>% as.matrix()
}

expect_that(probabilities_simple_mat[,1],equals(probabilities_simple_mat_check[,1]))
expect_that(probabilities_bayesian_mat[,1],equals(probabilities_bayesian_mat_check[,1]))
expect_that(probabilities_clustering_mat[,1],equals(probabilities_clustering_mat_check[,1]))

expect_that(probabilities_simple_mat[,n_imputations],
            equals(probabilities_simple_mat_check[,n_imputations]))
expect_that(probabilities_bayesian_mat[,n_imputations],
            equals(probabilities_bayesian_mat_check[,n_imputations]))
expect_that(probabilities_clustering_mat[,n_imputations],
            equals(probabilities_clustering_mat_check[,n_imputations]))

dataset <- dataset %>%
    dplyr::arrange(subject, allele, time)

dataset$probability_present <- rowMeans(imputed_datasets)

dataset <- dataset %>%
    dplyr::arrange(time, subject, allele)
dataset$probability_new <- rowMeans(probabilities_bayesian_mat, na.rm = T)

plot_dataset(dataset, treatments)







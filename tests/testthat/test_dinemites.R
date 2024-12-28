# library(testthat)
# library(dinemites)
#
# dataset <- read.csv(system.file(package="dinemites","extdata","dataset.tsv"), header = TRUE, sep="\t")
# treatments <- read.csv(system.file(package="dinemites","extdata","treatments.tsv"), header = TRUE, sep="\t")
# qPCR_only <- read.csv(system.file(package="dinemites","extdata","qPCR_only.tsv"), header = TRUE, sep="\t")
#
# dataset <- fill_in_dataset(dataset)
# dataset <- add_qpcr_times(dataset, qpcr_times = qPCR_only)
#
# estimate_drop_out(dataset)
#
# dataset <- dataset %>%
#     dplyr::arrange(subject, allele, time)
#
# n_imputations <- 10
# imputed_datasets <- impute_dataset(dataset, n_imputations = n_imputations)
#
# for (i in 1:n_imputations) {
#     dataset_tmp <- dataset
#     dataset_tmp$present <- imputed_datasets[,i]
#
#     dataset_tmp <- dataset_tmp %>%
#         add_present_infection() %>%
#         add_persistent_column() %>%
#         add_persistent_infection() %>%
#         add_lag_column() %>%
#         add_lag_infection() %>%
#         add_treatment_column(treatments = treatments, verbose = F) %>%
#         add_treatment_infection(treatments = treatments, verbose = F)
#
#     probabilities_simple <- determine_probabilities_simple(dataset_tmp)
#     probabilities_simple <- determine_probabilities_bayesian(
#         dataset = dataset_tmp,
#         infection_persistence_covariates = c("persistent_infection", "lag_infection_30", "treatment_acute_infection", "treatment_longitudinal_infection"),
#         infection_general_covariates = c("season"),
#         alleles_persistence_covariates = c("persistent", "lag_30", "treatment_acute", "treatment_longitudinal"),
#         refresh = 10)
# }
#
#
#
# probabilities_simple <- determine_probabilities_simple(dataset)
#
# taxa_table <- read.table(system.file(package="maaslin3","extdata","HMP2_taxonomy.tsv"), header = TRUE, sep="\t")
# metadata <- read.table(system.file(package="maaslin3","extdata","HMP2_metadata.tsv"), header = TRUE, sep="\t")
#
# metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
# metadata$dysbiosis_state <- factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
# metadata$antibiotics <- factor(metadata$antibiotics, levels = c('No', 'Yes'))
#
# # Run MaAsLin 3
# output_tmp <- tempfile()
# set.seed(1)
# fit_out <- maaslin3(input_data = taxa_table,
#                     input_metadata = metadata,
#                     output = output_tmp,
#                     normalization = 'TSS',
#                     transform = 'LOG',
#                     formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads',
#                     save_models = FALSE,
#                     plot_summary_plot = T,
#                     plot_associations = T,
#                     max_significance = 0.1,
#                     augment = TRUE,
#                     median_comparison_abundance = TRUE,
#                     median_comparison_prevalence = FALSE,
#                     cores=1,
#                     verbosity = 'WARN')
#
# maaslin_results = read.table(file.path(output_tmp, "significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)
#
# expect_that(expected_results_run1$metadata[1:50],equals(maaslin_results$metadata[1:50]))
# expect_that(expected_results_run1$feature[1:50],equals(maaslin_results$feature[1:50]))
# expect_that(round(expected_results_run1$N[1:50],10),equals(round(maaslin_results$N[1:50],10)))
# expect_that(round(as.numeric(expected_results_run1$pval_individual[1:50]),10),
#             equals(round(as.numeric(maaslin_results$pval_individual[1:50]),10)))
# expect_that(round(as.numeric(expected_results_run1$qval_individual[1:50]),10),
#             equals(round(as.numeric(maaslin_results$qval_individual[1:50]),10)))
#
# se <- SummarizedExperiment::SummarizedExperiment(
#     assays = list(taxa_table = t(taxa_table)),
#     colData = metadata
# )
#
# fit_out <- maaslin3(input_data = se,
#                     input_metadata = metadata,
#                     output = output_tmp,
#                     normalization = 'TSS',
#                     transform = 'LOG',
#                     formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads',
#                     save_models = FALSE,
#                     plot_summary_plot = T,
#                     plot_associations = T,
#                     max_significance = 0.1,
#                     augment = TRUE,
#                     median_comparison_abundance = TRUE,
#                     median_comparison_prevalence = FALSE,
#                     cores=1,
#                     verbosity = 'WARN')
#
# tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
#     assays = list(taxa_table = t(taxa_table)),
#     colData = metadata
# )
#
# fit_out <- maaslin3(input_data = tse,
#                     input_metadata = metadata,
#                     output = output_tmp,
#                     normalization = 'TSS',
#                     transform = 'LOG',
#                     formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads',
#                     save_models = FALSE,
#                     plot_summary_plot = T,
#                     plot_associations = T,
#                     max_significance = 0.1,
#                     augment = TRUE,
#                     median_comparison_abundance = TRUE,
#                     median_comparison_prevalence = FALSE,
#                     cores=1,
#                     verbosity = 'WARN')
#
# metadata <- as(metadata, "DataFrame")
# fit_out <- maaslin3(input_data = tse,
#                     input_metadata = metadata,
#                     output = output_tmp,
#                     normalization = 'TSS',
#                     transform = 'LOG',
#                     formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads',
#                     save_models = FALSE,
#                     plot_summary_plot = T,
#                     plot_associations = T,
#                     max_significance = 0.1,
#                     augment = TRUE,
#                     median_comparison_abundance = TRUE,
#                     median_comparison_prevalence = FALSE,
#                     cores=1,
#                     verbosity = 'WARN')
#
# unlink(output_tmp, recursive = T)

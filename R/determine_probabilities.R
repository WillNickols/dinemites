inv_logit <- function(x) {
    return(exp(x) / (1 + exp(x)))
}

# dataset must have columns for allele, subject, time, present,
# any general_covariates, and any persistence_covariates already as numerics

#' determine_probabilities_simple
#'
#' @export
#' @param dataset dataset
#' @param n_lags n_lags
#' @return return
#' @import dplyr
#'
determine_probabilities_simple <- function(dataset, n_lags = 3) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    # 2 if qPCR positive only
    if (any(!dataset$present %in% c(0,1,2))) {
        stop("present column must be 0/1/2")
    }

    ordered_dataset <- dataset %>%
        dplyr::arrange(time, subject, allele)

    if (any(dataset != ordered_dataset)) {
        stop("dataset must be arranged by time, subject, allele")
    }

    if (nrow(dataset) == 0) {
        return(list(probability_new = numeric(0), fit = NULL))
    }

    dataset$allele <- factor(dataset$allele)

    dataset <- dataset %>%
        dplyr::group_by(allele, subject) %>%
        dplyr::arrange(time) %>%
        dplyr::mutate(
            # Combine results from dynamically applied lags
            recent_present = purrr::reduce(
                purrr::map(1:n_lags, ~ lag(present, .x, default = 0) == 1),
                `|`
            ),
            # Apply logic based on recent_present
            probabilities = ifelse(recent_present & present == 1, 0, 1)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-recent_present) %>%
        dplyr::arrange(time, subject, allele)

    return(list(probability_new = dataset$probabilities, fit = NULL))
}

#' determine_probabilities_bayesian
#'
#' @export
#' @param dataset dataset
#' @param infection_persistence_covariates infection_persistence_covariates
#' @param infection_general_covariates infection_general_covariates
#' @param alleles_persistence_covariates alleles_persistence_covariates
#' @param chains chains
#' @param parallel_chains parallel_chains
#' @param iter_warmup iter_warmup
#' @param iter_sampling iter_sampling
#' @param refresh refresh
#' @param adapt_delta adapt_delta
#' @param seed seed
#' @param drop_out drop_out
#' @return An list with probabilities
#' @import dplyr
#' @importFrom  stats as.formula
#' @importFrom  stats model.matrix
#'
determine_probabilities_bayesian <- function(dataset,
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
                                    drop_out = FALSE) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    # 2 if qPCR positive only
    if (any(!dataset$present %in% c(0,1))) {
        stop("present column must be 0/1")
    }

    if (!is.null(infection_persistence_covariates) &&
        !is.numeric(c(unlist(dataset[infection_persistence_covariates])))) {
        stop("infection_persistence_covariates columns must be numeric")
    }

    if (!is.null(infection_general_covariates) &&
        !is.numeric(c(unlist(dataset[infection_general_covariates])))) {
        stop("infection_general_covariates columns must be numeric")
    }

    if (!is.null(alleles_persistence_covariates) &&
        !is.numeric(c(unlist(dataset[alleles_persistence_covariates])))) {
        stop("alleles_persistence_covariates columns must be numeric")
    }

    ordered_dataset <- dataset %>%
        dplyr::arrange(time, subject, allele)

    if (any(dataset != ordered_dataset)) {
        stop("dataset must be arranged by time, subject, allele")
    }

    unique_alleles_len <- length(unique(dataset$allele))
    i <- 1

    checked <- TRUE
    while (i <= nrow(dataset)) {
        if (dataset$present[i] == 2) {
            if (!all(dataset$present[i:(i + unique_alleles_len - 1)] == 2)) {
                checked <- FALSE
                break
            }
            i <- i + unique_alleles_len
        } else {
            i <- i + 1
        }
    }
    if (!checked) {
        stop("dataset not valid: present = 2 should come in blocks of size n_alleles")
    }

    subject_time_combinations <- unique(dataset[, c("subject", "time")])

    # Iterate over the actual existing combinations
    for (i in 1:nrow(subject_time_combinations)) {
        subject_uniq <- subject_time_combinations$subject[i]
        time_uniq <- subject_time_combinations$time[i]
        if (any(apply(dataset[dataset$time == time_uniq & dataset$subject == subject_uniq,
                              c(infection_persistence_covariates, infection_general_covariates), drop=F],
                      2, function(x){length(unique(x))})) != 1) {
            stop("dataset not valid: infection_persistence_covariates and
                 infection_general_covariates should be the same at a person-time")
        }
    }

    if (nrow(dataset) == 0) {
        return(list(probability_new = numeric(0), fit = NULL))
    }

    dataset$allele <- factor(dataset$allele)

    # Check for possible speed-up from grouping observations
    columns_for_reducing <- c("allele", "present", "present_infection", infection_general_covariates,
                              infection_persistence_covariates, alleles_persistence_covariates)

    dataset_reduced <- dataset %>%
        distinct(across(all_of(columns_for_reducing))) %>%
        left_join(dataset %>% count(across(all_of(columns_for_reducing))), by = columns_for_reducing)

    tryCatch({
        mm <- model.matrix(as.formula(paste0(c("present ~ 0 + ",
                                               c(infection_persistence_covariates, infection_general_covariates)),
                                             collapse = " + ")),
                           dataset_reduced)

        if (!is.null(infection_general_covariates)) {
            mm <- model.matrix(as.formula(paste0(c("present ~ 0 + ",
                                                   c(infection_general_covariates)),
                                                 collapse = " + ")),
                               dataset_reduced)
        }

        mm <- model.matrix(as.formula(paste0(c("present ~ 0 + allele",
                                               c(alleles_persistence_covariates)),
                                             collapse = " + allele:")),
                           dataset_reduced)
    },
    error = function(e) {
        stop(e)
    })

    X_infection_general <- as.matrix(dataset_reduced[infection_general_covariates])
    X_infection_persistent <- as.matrix(dataset_reduced[infection_persistence_covariates])
    X_alleles_persistent <- as.matrix(dataset_reduced[alleles_persistence_covariates])

    N <- nrow(dataset_reduced)
    K_infection_general <- ncol(X_infection_general)
    K_infection_persistent <- ncol(X_infection_persistent)
    K_alleles_persistent <- ncol(X_alleles_persistent)
    J <- length(unique(dataset_reduced$allele))

    group <- as.numeric(dataset_reduced$allele)

    z <- dataset_reduced$present_infection
    y <- dataset_reduced$present
    w <- dataset_reduced$n

    # Prepare data for Stan
    stan_data <- list(
        N = N,
        K_infection_general = K_infection_general,
        K_infection_persistent = K_infection_persistent,
        K_alleles_persistent = K_alleles_persistent,
        J = J,
        X_infection_general = X_infection_general,
        X_infection_persistent = X_infection_persistent,
        X_alleles_persistent = X_alleles_persistent,
        group = group,
        z = z,
        y = y,
        w = w
    )

    # Compile the model
    if (!drop_out) {
        mod <- instantiate::stan_package_model(name = "model_infection_probabilities_bayesian", package = "dinemites")
    } else {
        mod <- instantiate::stan_package_model(name = "model_infection_probabilities_bayesian_drop_out", package = "dinemites")
    }


    # Fit the model
    fit <- mod$sample(
        data = stan_data,
        chains = chains,
        parallel_chains = parallel_chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        seed = seed,
        adapt_delta = adapt_delta,
        refresh = refresh
    )

    output_draws <- data.frame(posterior::as_draws_df(fit$draws()), check.names = F)

    # p_new_infection
    mm_formula <- ifelse(length(infection_general_covariates) > 0,
                         paste0(c("present ~ 1 ", c(infection_general_covariates)), collapse = " + "),
                         'present ~ 1')

    dataset_reduced <- dataset_reduced %>%
        mutate(original_row_index = row_number())

    # Join back with the original dataset to map indices
    dataset_with_indices <- dataset %>%
        left_join(dataset_reduced %>%
                      select(all_of(columns_for_reducing), original_row_index),
                  by = columns_for_reducing)

    mm <- model.matrix(as.formula(mm_formula), dataset_reduced)
    p_new_infection <- matrix(nrow = nrow(output_draws), ncol = nrow(dataset_reduced))
    for (row_num in seq(nrow(output_draws))) {
        coefs <- output_draws[row_num, "alpha_infection"]
        if (length(infection_general_covariates) > 0) {
            for (i in 1:length(infection_general_covariates)) {
                coefs <- c(coefs,
                           unlist(output_draws[row_num,
                                               grepl(paste0("^beta_infection_general\\[", i, "\\]"),
                                                     colnames(output_draws))]))
            }
        }
        p_new_infection[row_num,] <- inv_logit(mm %*% coefs)
    }

    # p_any_infection
    mm <- model.matrix(as.formula(paste0(c("present ~ 1 + ",
                                           c(infection_persistence_covariates, infection_general_covariates)),
                                         collapse = " + ")),
                       dataset_reduced)
    p_any_infection <- matrix(nrow = nrow(output_draws), ncol = nrow(dataset_reduced))
    for (row_num in seq(nrow(output_draws))) {
        coefs <- output_draws[row_num, "alpha_infection"]
        if (length(infection_persistence_covariates) > 0) {
            for (i in 1:length(infection_persistence_covariates)) {
                coefs <- c(coefs,
                           unlist(output_draws[row_num,
                                               grepl(paste0("^beta_infection_persistent\\[", i, "\\]"),
                                                     colnames(output_draws))]))
            }
        }
        if (length(infection_general_covariates) > 0) {
            for (i in 1:length(infection_general_covariates)) {
                coefs <- c(coefs,
                           unlist(output_draws[row_num,
                                               grepl(paste0("^beta_infection_general\\[", i, "\\]"),
                                                     colnames(output_draws))]))
            }
        }

        p_any_infection[row_num,] <- inv_logit(mm %*% coefs)
    }

    # p_old_infection
    p_old_infection = (p_any_infection - p_new_infection) / (1 - p_new_infection)

    # p_allele_if_new
    mm <- model.matrix(as.formula(paste0(c("present ~ 0 + allele"))), dataset_reduced)
    p_allele_if_new <- matrix(nrow = nrow(output_draws), ncol = nrow(dataset_reduced))
    for (row_num in seq(nrow(output_draws))) {
        coefs <- unlist(output_draws[row_num, grepl(paste0("^alpha_alleles_new"),
                                                    colnames(output_draws))])

        p_allele_if_new[row_num,] <- inv_logit(mm %*% coefs)
    }

    # p_allele_if_old
    mm <- model.matrix(as.formula(paste0(c("present ~ 0 ",
                                           c(alleles_persistence_covariates)),
                                         collapse = " + allele:")),
                       dataset_reduced)
    p_allele_if_old <- matrix(nrow = nrow(output_draws), ncol = nrow(dataset_reduced))

    if (drop_out) {
        mm_missed <- model.matrix(as.formula(paste0(c("present ~ 0 + allele"))), dataset_reduced)
        p_allele_if_old_missed <- matrix(nrow = nrow(output_draws), ncol = nrow(dataset_reduced))
    }
    for (row_num in seq(nrow(output_draws))) {
        coefs <- c()
        if (length(alleles_persistence_covariates) > 0) {
            for (i in 1:length(alleles_persistence_covariates)) {
                coefs <- c(coefs,
                           unlist(output_draws[row_num,
                                               grepl(paste0("^beta_alleles_old\\[.*\\,", i, "\\]"),
                                                     colnames(output_draws))]))
            }
        }
        p_allele_if_old[row_num,] <- inv_logit(mm %*% coefs)

        if (drop_out) {
            coefs_missed <- unlist(output_draws[row_num, grepl(paste0("^alpha_alleles_new"),
                                                               colnames(output_draws))])

            p_allele_if_old_missed[row_num,] <- inv_logit(mm_missed %*% coefs_missed) *
                mm_missed %*% unlist(output_draws[row_num, grepl(paste0("^drop_out"), colnames(output_draws))])

            p_allele_if_old[row_num,] <- ifelse(rowSums(abs(dataset_reduced[alleles_persistence_covariates])) > 0,
                                                p_allele_if_old[row_num,], p_allele_if_old_missed[row_num,])
        } else {
            p_allele_if_old[row_num,] <- ifelse(rowSums(abs(dataset_reduced[alleles_persistence_covariates])) > 0,
                                                p_allele_if_old[row_num,], 0)
        }
    }

    p_allele_new = p_allele_if_new * p_new_infection
    p_allele_old = p_allele_if_old * p_old_infection

    denominator <- p_allele_new + p_allele_old - p_allele_new * p_allele_old
    numerator <- p_allele_new

    probabilities_new <- colMeans(numerator/denominator)

    probabilities_new <- probabilities_new[dataset_with_indices$original_row_index]

    return(list(probability_new = probabilities_new, fit = fit))
}

#' determine_probabilities_clustering
#'
#' @export
#' @param dataset dataset
#' @param refresh refresh
#' @param seed seed
#' @return An list with probabilities
#' @import dplyr
#' @importFrom  stats as.formula
#' @importFrom  stats model.matrix
#' @importFrom utils combn
#'
determine_probabilities_clustering <- function(dataset,
                                               refresh = 100,
                                               seed = 1) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    # 2 if qPCR positive only
    if (any(!dataset$present %in% c(0,1,2))) {
        stop("present column must be 0/1/2")
    }

    if (any(dataset$present == 2)) {
        warning("2s set to 0s when using the clustering method; please impute")
        dataset$present[dataset$present == 2] <- 0
    }

    ordered_dataset <- dataset %>%
        dplyr::arrange(time, subject, allele)

    if (any(dataset != ordered_dataset)) {
        stop("dataset must be arranged by time, subject, allele")
    }

    unique_alleles_len <- length(unique(dataset$allele))
    i <- 1

    if (nrow(dataset) == 0) {
        return(list(probability_new = numeric(0), fit = NULL))
    }

    dataset$allele <- factor(dataset$allele)

    mod <- instantiate::stan_package_model(name = "model_infection_probabilities_clusters", package = "dinemites")

    dataset$probabilities <- 0
    for (subject_name in unique(dataset$subject)) {
        present_df <- dataset[dataset$present == 1 & dataset$subject == subject_name,]
        if (nrow(present_df) == 0) {
            next
        }

        allele_counts <- table(interaction(present_df$subject, present_df$time), present_df$allele)
        sample_allele_counts <- rowSums(allele_counts > 0)

        samples_with_2plus_alleles <- names(sample_allele_counts[sample_allele_counts >= 2])
        filtered_df <- present_df[interaction(present_df$subject, present_df$time) %in% samples_with_2plus_alleles, ]

        allele_pairs <- list()
        for (sample in unique(interaction(filtered_df$subject, filtered_df$time))) {
            alleles_in_sample <- filtered_df$allele[interaction(filtered_df$subject, filtered_df$time) == sample]
            if (length(alleles_in_sample) >= 2) {
                pairs <- combn(alleles_in_sample, 2, simplify = FALSE)
                allele_pairs <- c(allele_pairs, pairs)
            }
        }

        edges <- t(as.data.frame(allele_pairs))
        if (nrow(edges) == 0) {
            edges <- matrix(nrow = 0, ncol = 2)
        }
        colnames(edges) <- c("allele1", "allele2")

        counts <- as.data.frame(table(apply(edges, 1, paste, collapse = "_")))
        counts_split <- do.call(rbind, strsplit(as.character(counts$Var1), "_"))
        if (length(counts$Freq) > 0) {
            edges_with_counts <- cbind(counts_split, Freq = counts$Freq)
            colnames(edges_with_counts) <- c(colnames(edges), "Count")
            edges_with_counts <- as.data.frame(edges_with_counts)
            edges_with_counts$Count <- as.numeric(as.character(edges_with_counts$Count))
        } else {
            edges_with_counts <- matrix(nrow = 0, ncol = 3)
            colnames(edges_with_counts) <- c(colnames(edges), "Count")
            edges_with_counts <- as.data.frame(edges_with_counts)
            edges_with_counts$Count <- as.numeric(as.character(edges_with_counts$Count))
        }

        edges <- edges_with_counts

        if (nrow(edges) > 1) {
            lc <- tryCatch({
                linkcomm::getLinkCommunities(edges, plot = F, verbose = F)
            }, error = function(e) {
                if (grepl('no clusters were found in this network', e$message)) {
                    lc <- list()
                    lc$nodeclusters <- data.frame(node = (unique(unlist(edges[,c(1,2)]))))
                    lc$nodeclusters$cluster <- 1:nrow(lc$nodeclusters)
                    return(lc)
                }
            })
        } else if (nrow(edges) == 1) {
            lc <- list()
            lc$nodeclusters <- data.frame(node = unlist(edges[1,c(1,2)]),
                                          cluster = 1)
        } else {
            lc <- list()
            lc$nodeclusters <- data.frame(matrix(nrow = 0, ncol = 2))
            colnames(lc$nodeclusters) <- c("node", "cluster")
        }

        lc$nodeclusters$cluster <- as.numeric(lc$nodeclusters$cluster)
        if (any(!present_df$allele %in% lc$nodeclusters$node)) {
            df_to_append <- data.frame(
                node = unique(present_df$allele)[
                    !unique(present_df$allele) %in% lc$nodeclusters$node])
            df_to_append$cluster <- (max(lc$nodeclusters$cluster, 0) + 1):(
                max(lc$nodeclusters$cluster, 0) + nrow(df_to_append))
            lc$nodeclusters <- rbind(lc$nodeclusters, df_to_append)
        }

        df_for_fit <- dataset[dataset$allele %in% lc$nodeclusters$node & dataset$subject == subject_name,]

        N = nrow(df_for_fit)
        J = length(unique(df_for_fit$allele))
        n_times = length(unique(df_for_fit$time))
        C = length(unique(lc$nodeclusters$cluster))

        df_for_fit$allele <- factor(as.character(df_for_fit$allele))
        alleles = as.numeric(df_for_fit$allele)
        y = df_for_fit$present
        times = as.numeric(factor(df_for_fit$time))

        a_levels <- unique(lc$nodeclusters[, 1])
        b_levels <- unique(lc$nodeclusters[, 2])

        alleles_in_clusters <- matrix(0, nrow = length(a_levels), ncol = length(b_levels),
                                      dimnames = list(a_levels, b_levels))

        for (i in 1:nrow(lc$nodeclusters)) {
            a <- lc$nodeclusters[i, 1]
            b <- lc$nodeclusters[i, 2]
            alleles_in_clusters[as.character(a), b] <- 1
        }

        alleles_in_clusters <- alleles_in_clusters[match(levels(df_for_fit$allele), rownames(alleles_in_clusters)),,drop=F]

        # Prepare data for Stan
        stan_data <- list(
            N = N,
            J = J,
            n_times = n_times,
            C = C,
            alleles = alleles,
            y = y,
            times = times,
            alleles_in_clusters = alleles_in_clusters
        )

        # Fit the model
        fit <- mod$optimize(
            data = stan_data,
            seed = seed,
            refresh = refresh,
            show_messages = F
        )

        output_draws <- data.frame(posterior::as_draws_df(fit$draws()), check.names = F)

        prob_cluster_in_time = matrix(unlist(output_draws[grepl("prob_cluster_in_time", colnames(output_draws))]),
                                      nrow = n_times, ncol = C)
        prob_allele_in_cluster = matrix(unlist(output_draws[grepl("prob_allele_in_cluster", colnames(output_draws))]),
                                        nrow = J, ncol = C)

        prob_allele_in_cluster = prob_allele_in_cluster * alleles_in_clusters

        # alleles x clusters x times
        p_j_in_c_at_t_given_a_present <- array(NA, dim = c(J, C, n_times))

        for (time_point in 1:n_times) {
            for (allele_index in 1:J) {
                denom <- 1 - prod(1 - prob_allele_in_cluster[allele_index,] * prob_cluster_in_time[time_point,])
                p_j_in_c_at_t_given_a_present[allele_index, , time_point] <-
                    prob_allele_in_cluster[allele_index,] * prob_cluster_in_time[time_point,] /
                    denom
            }
        }
        p_j_in_c_at_t_given_a_present <- pmin(p_j_in_c_at_t_given_a_present, 1)

        prob_new_vec <- rep(0, N)
        for (j in 1:N) {
            time_point <- times[j]
            allele_index <- alleles[j]

            p_j_present = 1 - prod((1 - prob_allele_in_cluster[allele_index,]) *
                                       prob_cluster_in_time[time_point,] +
                                       (1 - prob_cluster_in_time[time_point,]))

            p_j_not_in_c_at_previous_times <- vector(length = C)
            for (cluster in 1:C) {
                previously_occurring_times <- times[alleles == allele_index & times < time_point & y == 1]
                if (length(previously_occurring_times) > 0) {
                    p_j_not_in_c_at_previous_times[cluster] <-
                        prod(1 - p_j_in_c_at_t_given_a_present[allele_index, cluster, previously_occurring_times])
                } else {
                    p_j_not_in_c_at_previous_times[cluster] <- 1
                }
            }

            p_j_present_and_new = 1 - prod((1 - p_j_not_in_c_at_previous_times *
                                                prob_allele_in_cluster[allele_index,]) *
                                               prob_cluster_in_time[time_point,] +
                                               (1 - prob_cluster_in_time[time_point,]))

            p_j_new_given_a_present <- p_j_present_and_new / p_j_present

            prob_new_vec[j] <- p_j_new_given_a_present
        }

        dataset$probabilities[dataset$allele %in% lc$nodeclusters$node & dataset$subject == subject_name] <-
            prob_new_vec
    }

    return(list(probability_new = dataset$probabilities, fit = NULL))
}


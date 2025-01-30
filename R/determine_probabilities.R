inv_logit <- function(x) {
    return(ifelse(x > 700, 1, exp(x) / (1 + exp(x))))
}

#' determine_probabilities_simple
#'
#' Determine whether an allele is new or persistent using a simple rule.
#' For each subject, an allele is counted as new if it has not been observed in
#' the `n_lags` most recent samples (e.g., visits) or in the `t_lag` most recent
#' times (e.g., days). If the allele has been observed in both the last
#' `n_lags` samples and the last `t_lag` times, it is considered new.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present.`
#' @param n_lags If the allele has been observed in both the last
#' `n_lags` samples and the last `t_lag` times, it is considered new.
#' @param t_lag If the allele has been observed in both the last
#' `n_lags` samples and the last `t_lag` times, it is considered new.
#' @return A named list with (1) `probability_new` as a 0/1 vector of whether
#' the corresponding allele in the dataset was from a new infection and (2)
#' `fit` as `NULL` (for consistency with the other models).
#' @examples
#'
#' dataset <- data.frame(allele = rep(c('A', 'B', 'C', 'D', 'E'), 5),
#'     subject = rep('A', 25),
#'     time = rep(1:5, each = 5),
#'     present = c(0,0,0,0,0, 1,0,0,1,0, 1,0,0,0,0, 1,0,0,0,0, 0,0,0,0,0))
#'
#' dataset$probability_new <-
#'     determine_probabilities_simple(dataset)$probability_new
#'
#' @import dplyr
#'
determine_probabilities_simple <- function(dataset,
                                           n_lags = 3,
                                           t_lag = Inf) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    if (n_lags < 1) {
        stop("n_lags must be >= 1")
    }

    if (t_lag < 0) {
        stop("t_lag must be >= 0")
    }

    if (any(!dataset$present %in% c(0,1))) {
        stop("present column must be 0/1")
    }

    dataset <- dataset %>%
        dplyr::mutate(original_row_ordering = seq(nrow(dataset))) %>%
        dplyr::arrange(.data$time, .data$subject, .data$allele)

    if (nrow(dataset) == 0) {
        return(list(probability_new = numeric(0), fit = NULL))
    }

    dataset$allele <- factor(dataset$allele)

    dataset <- dataset %>%
        dplyr::group_by(.data$allele, .data$subject) %>%
        dplyr::arrange(.data$time) %>%
        dplyr::mutate(
            # Combine results from dynamically applied lags
            recent_present = purrr::reduce(
                purrr::map(1:n_lags,
                           ~ lag(.data$present, .x, default = 0) == 1),
                `|`
            ),
            lag_present = ifelse(sapply(.data$time, function(t) {
                any(.data$present[.data$time >= (t - t_lag) &
                                      .data$time < t] == 1)
            }), 1, 0),
            probabilities =
                ifelse(.data$present == 1,
                       ifelse(.data$lag_present &
                                  .data$recent_present, 0, 1), NA)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-.data$recent_present) %>%
        dplyr::arrange(.data$original_row_ordering) %>%
        dplyr::select(-.data$original_row_ordering)

    return(list(probability_new = dataset$probabilities, fit = NULL))
}

#' determine_probabilities_bayesian
#'
#' Determine whether an allele is new or persistent using a Bayesian model.
#' A Bayesian model is parameterized with general and persistence covariates
#' to model the probability an infection (overall and allele-specific) is
#' present. The model is fit in Stan, and the parameters are used to estimate
#' the probability each allele is new.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present.`
#' @param infection_persistence_covariates A character vector of column names
#' in the dataset corresponding to covariates that increase only the
#' probability of an infection being present due to a persistent infection.
#' @param infection_general_covariates A character vector of column names
#' in the dataset corresponding to covariates that increase the
#' probability of a new infection being present
#' @param alleles_persistence_covariates A character vector of column names
#' in the dataset corresponding to covariates that increase only the
#' probability of an allele being present due to a persistent infection.
#' @param chains Number of chains for the Stan model
#' @param parallel_chains Number of parallel chains for the Stan model
#' @param iter_warmup Number of warm-up steps for the Stan model
#' @param iter_sampling Number of sampling steps for the Stan model
#' @param refresh Stan updates will be printed every `refresh` steps
#' @param adapt_delta Stan's `adapt_delta` parameter
#' @param seed Random seed for Stan
#' @param drop_out Whether to use a model that allows alleles to be called
#' persistent despite having never been previously observed due to sequencing
#' drop-out. This is recommended if more than 50% of the alleles have dropped
#' out.
#' @return A named list with (1) `probability_new` as a vector of the
#' probabilities the corresponding alleles in the dataset were from new
#' infections and (2) `fit` as the Stan model output.
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

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

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

    dataset <- dataset %>%
        dplyr::mutate(original_row_ordering = seq(nrow(dataset))) %>%
        dplyr::arrange(.data$time, .data$subject, .data$allele)

    subject_time_combinations <- unique(dataset[, c("subject", "time")])

    # Iterate over the actual existing combinations
    for (i in 1:nrow(subject_time_combinations)) {
        subject_uniq <- subject_time_combinations$subject[i]
        time_uniq <- subject_time_combinations$time[i]
        if (any(apply(dataset[
            dataset$time == time_uniq & dataset$subject == subject_uniq,
            c(infection_persistence_covariates, infection_general_covariates),
            drop=F],
                      2, function(x){length(unique(x))})) != 1) {
            stop("dataset not valid: infection_persistence_covariates and
                 infection_general_covariates should be the same
                 at a subject-time")
        }
    }

    if (nrow(dataset) == 0) {
        return(list(probability_new = numeric(0), fit = NULL))
    }

    dataset$allele <- factor(dataset$allele)

    # Check for possible speed-up from grouping observations
    columns_for_reducing <-
        c("allele", "present", "present_infection",
          infection_general_covariates,
          infection_persistence_covariates,
          alleles_persistence_covariates)

    dataset_reduced <- dataset %>%
        distinct(across(all_of(columns_for_reducing))) %>%
        left_join(dataset %>% count(across(all_of(columns_for_reducing))),
                  by = columns_for_reducing)

    tryCatch({
        mm <- model.matrix(as.formula(paste0(
            c("present ~ 0 + ",
              c(infection_persistence_covariates,
                infection_general_covariates)),
            collapse = " + ")),
            dataset_reduced)

        if (!is.null(infection_general_covariates)) {
            mm <- model.matrix(
                as.formula(paste0(c("present ~ 0 + ",
                                    c(infection_general_covariates)),
                                  collapse = " + ")),
                               dataset_reduced)
        }

        mm <- model.matrix(
            as.formula(paste0(c("present ~ 0 + allele",
                                c(alleles_persistence_covariates)),
                              collapse = " + allele:")),
                           dataset_reduced)
    },
    error = function(e) {
        stop(e)
    })

    X_infection_general <-
        as.matrix(dataset_reduced[infection_general_covariates])
    X_infection_persistent <-
        as.matrix(dataset_reduced[infection_persistence_covariates])
    X_alleles_persistent <-
        as.matrix(dataset_reduced[alleles_persistence_covariates])

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
        mod <- instantiate::stan_package_model(
            name = "model_infection_probabilities_bayesian",
            package = "dinemites")
    } else {
        mod <- instantiate::stan_package_model(
            name = "model_infection_probabilities_bayesian_drop_out",
            package = "dinemites")
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

    output_draws <- data.frame(posterior::as_draws_df(fit$draws()),
                               check.names = F)

    # Get indices for mapping back reduced dataset to the original
    dataset_reduced <- dataset_reduced %>%
        mutate(original_row_index = row_number())

    # Join back with the original dataset to map indices
    dataset_with_indices <- dataset %>%
        left_join(dataset_reduced %>%
                      select(all_of(columns_for_reducing),
                             .data$original_row_index),
                  by = columns_for_reducing)

    # p_new_infection
    mm_formula <- ifelse(length(infection_general_covariates) > 0,
                         paste0(c("present ~ 1 ",
                                  c(infection_general_covariates)),
                                collapse = " + "),
                         'present ~ 1')

    mm <- model.matrix(as.formula(mm_formula), dataset_reduced)
    p_new_infection <- matrix(nrow = nrow(output_draws),
                              ncol = nrow(dataset_reduced))
    for (row_num in seq(nrow(output_draws))) {
        coefs <- output_draws[row_num, "alpha_infection"]
        if (length(infection_general_covariates) > 0) {
            for (i in 1:length(infection_general_covariates)) {
                coefs <-
                    c(coefs,
                      unlist(output_draws[row_num,
                          grepl(paste0("^beta_infection_general\\[", i, "\\]"),
                                colnames(output_draws))]))
            }
        }
        p_new_infection[row_num,] <- inv_logit(mm %*% coefs)
    }

    # p_any_infection
    mm <- model.matrix(
        as.formula(paste0(c("present ~ 1 + ",
            c(infection_persistence_covariates, infection_general_covariates)),
            collapse = " + ")),
        dataset_reduced)
    p_any_infection <- matrix(nrow = nrow(output_draws),
                              ncol = nrow(dataset_reduced))
    for (row_num in seq(nrow(output_draws))) {
        coefs <- output_draws[row_num, "alpha_infection"]
        if (length(infection_persistence_covariates) > 0) {
            for (i in 1:length(infection_persistence_covariates)) {
                coefs <-
                    c(coefs,
                      unlist(output_draws[row_num,
                      grepl(paste0("^beta_infection_persistent\\[", i, "\\]"),
                            colnames(output_draws))]))
            }
        }
        if (length(infection_general_covariates) > 0) {
            for (i in 1:length(infection_general_covariates)) {
                coefs <-
                    c(coefs,
                      unlist(output_draws[row_num,
                      grepl(paste0("^beta_infection_general\\[", i, "\\]"),
                            colnames(output_draws))]))
            }
        }

        p_any_infection[row_num,] <- inv_logit(mm %*% coefs)
    }

    # p_old_infection
    p_old_infection = (p_any_infection - p_new_infection) /
        (1 - p_new_infection)

    # p_allele_if_new
    mm <- model.matrix(
        as.formula(paste0(c("present ~ 0 + allele"))), dataset_reduced)
    p_allele_if_new <- matrix(nrow = nrow(output_draws),
                              ncol = nrow(dataset_reduced))
    for (row_num in seq(nrow(output_draws))) {
        coefs <- unlist(output_draws[row_num,
                                     grepl(paste0("^alpha_alleles_new"),
                                           colnames(output_draws))])

        p_allele_if_new[row_num,] <- inv_logit(mm %*% coefs)
    }

    # p_allele_if_old
    mm <- model.matrix(as.formula(paste0(c("present ~ 0 ",
                                           c(alleles_persistence_covariates)),
                                         collapse = " + allele:")),
                       dataset_reduced)
    p_allele_if_old <- matrix(nrow = nrow(output_draws),
                              ncol = nrow(dataset_reduced))

    if (drop_out) {
        mm_missed <- model.matrix(
            as.formula(paste0(c("present ~ 0 + allele"))), dataset_reduced)
        p_allele_if_old_missed <- matrix(nrow = nrow(output_draws),
                                         ncol = nrow(dataset_reduced))
    }
    for (row_num in seq(nrow(output_draws))) {
        coefs <- c()
        if (length(alleles_persistence_covariates) > 0) {
            for (i in 1:length(alleles_persistence_covariates)) {
                coefs <-
                    c(coefs,
                      unlist(output_draws[row_num,
                          grepl(paste0("^beta_alleles_old\\[.*\\,", i, "\\]"),
                                colnames(output_draws))]))
            }
        }
        p_allele_if_old[row_num,] <- inv_logit(mm %*% coefs)

        if (drop_out) {
            coefs_missed <-
                unlist(output_draws[row_num, grepl(paste0("^alpha_alleles_new"),
                                                   colnames(output_draws))])

            p_allele_if_old_missed[row_num,] <-
                inv_logit(mm_missed %*% coefs_missed) *
                mm_missed %*%
                unlist(output_draws[row_num, grepl(paste0("^drop_out"),
                                                   colnames(output_draws))])

            p_allele_if_old[row_num,] <-
                ifelse(rowSums(abs(
                    dataset_reduced[alleles_persistence_covariates])) > 0,
                    p_allele_if_old[row_num,], p_allele_if_old_missed[row_num,])
        } else {
            p_allele_if_old[row_num,] <-
                ifelse(rowSums(abs(
                    dataset_reduced[alleles_persistence_covariates])) > 0,
                    p_allele_if_old[row_num,], 0)
        }
    }

    p_allele_new = p_allele_if_new * p_new_infection
    p_allele_old = p_allele_if_old * p_old_infection

    denominator <- p_allele_new + p_allele_old - p_allele_new * p_allele_old
    numerator <- p_allele_new

    probabilities_new <- colMeans(numerator/denominator)

    probabilities_new <-
        probabilities_new[dataset_with_indices$original_row_index]
    probabilities_new <-
        ifelse(dataset_with_indices$present == 1, probabilities_new, NA)
    probabilities_new <-
        probabilities_new[order(dataset$original_row_ordering)]

    return(list(probability_new = probabilities_new, fit = fit))
}

#' determine_probabilities_clustering
#'
#' Determine whether an allele is new or persistent using a clustering model.
#' Alleles are clustered according to co-occurrence, and the probability of
#' each allele occurring in each cluster and each cluster occurring at each time
#' is found through a likelihood maximization procedure. These probabilities
#' are used to estimate the probability an infection is new.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present.`
#' @param refresh Stan updates will be printed every `refresh` steps
#' @param seed Random seed for Stan
#' @return A named list with (1) `probability_new` as a vector of the
#' probabilities the corresponding alleles in the dataset were from new
#' infections and (2)
#' `fit` as the resulting parameters from the likelihood maximization
#' procedure.
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

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    if (any(!dataset$present %in% c(0,1))) {
        stop("present column must be 0/1")
    }

    dataset <- dataset %>%
        dplyr::mutate(original_row_ordering = seq(nrow(dataset))) %>%
        dplyr::arrange(.data$time, .data$subject, .data$allele)

    if (nrow(dataset) == 0) {
        return(list(probability_new = numeric(0), fit = NULL))
    }

    dataset$allele <- factor(dataset$allele)

    mod <- instantiate::stan_package_model(
        name = "model_infection_probabilities_clusters",
        package = "dinemites")

    dataset$probabilities <- 0
    for (subject_name in unique(dataset$subject)) {
        present_df <- dataset[dataset$present == 1 &
                                  dataset$subject == subject_name,]
        if (nrow(present_df) == 0) {
            next
        }

        allele_counts <-
            table(interaction(present_df$subject, present_df$time),
                  present_df$allele)
        sample_allele_counts <- rowSums(allele_counts > 0)

        samples_with_2plus_alleles <-
            names(sample_allele_counts[sample_allele_counts >= 2])
        filtered_df <- present_df[
            interaction(present_df$subject, present_df$time) %in%
                samples_with_2plus_alleles, ]

        allele_pairs <- list()
        for (sample in unique(interaction(filtered_df$subject,
                                          filtered_df$time))) {
            alleles_in_sample <- filtered_df$allele[interaction(
                filtered_df$subject, filtered_df$time) == sample]
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

        counts <- as.data.frame(table(apply(edges, 1, paste, collapse = "_jOiNiNgStRiNg_")))
        counts_split <- do.call(rbind, strsplit(as.character(counts$Var1), "_jOiNiNgStRiNg_"))
        if (length(counts$Freq) > 0) {
            edges_with_counts <- cbind(counts_split, Freq = counts$Freq)
            colnames(edges_with_counts) <- c(colnames(edges), "Count")
            edges_with_counts <- as.data.frame(edges_with_counts)
            edges_with_counts$Count <-
                as.numeric(as.character(edges_with_counts$Count))
        } else {
            edges_with_counts <- matrix(nrow = 0, ncol = 3)
            colnames(edges_with_counts) <- c(colnames(edges), "Count")
            edges_with_counts <- as.data.frame(edges_with_counts)
            edges_with_counts$Count <-
                as.numeric(as.character(edges_with_counts$Count))
        }

        edges <- edges_with_counts

        if (nrow(edges) > 1) {
            lc <- tryCatch({
                linkcomm::getLinkCommunities(edges, plot = F, verbose = F)
            }, error = function(e) {
                if (grepl('no clusters were found in this network',
                          e$message)) {
                    lc <- list()
                    lc$nodeclusters <-
                        data.frame(node = (unique(unlist(edges[,c(1,2)]))))
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

        df_for_fit <- dataset[dataset$allele %in% lc$nodeclusters$node &
                                  dataset$subject == subject_name,]

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

        alleles_in_clusters <- matrix(0, nrow = length(a_levels),
                                      ncol = length(b_levels),
                                      dimnames = list(a_levels, b_levels))

        for (i in 1:nrow(lc$nodeclusters)) {
            a <- lc$nodeclusters[i, 1]
            b <- lc$nodeclusters[i, 2]
            alleles_in_clusters[as.character(a), as.character(b)] <- 1
        }

        alleles_in_clusters <-
            alleles_in_clusters[match(levels(df_for_fit$allele),
                                      rownames(alleles_in_clusters)), , drop=F]

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

        output_draws <- data.frame(posterior::as_draws_df(fit$draws()),
                                   check.names = F)

        prob_cluster_in_time =
            matrix(unlist(output_draws[grepl("prob_cluster_in_time",
                                             colnames(output_draws))]),
                                      nrow = n_times, ncol = C)
        prob_allele_in_cluster =
            matrix(unlist(output_draws[grepl("prob_allele_in_cluster",
                                             colnames(output_draws))]),
                                        nrow = J, ncol = C)

        prob_allele_in_cluster = prob_allele_in_cluster * alleles_in_clusters

        # alleles x clusters x times
        p_a_in_c_at_t_given_a_present <- array(NA, dim = c(J, C, n_times))

        for (time_point in seq(n_times)) {
            for (allele_index in 1:J) {
                denom <- 1 - prod(1 - prob_allele_in_cluster[allele_index,] *
                                      prob_cluster_in_time[time_point,])
                p_a_in_c_at_t_given_a_present[allele_index, , time_point] <-
                    prob_allele_in_cluster[allele_index,] *
                    prob_cluster_in_time[time_point,] /
                    denom
            }
        }
        p_a_in_c_at_t_given_a_present <- pmin(p_a_in_c_at_t_given_a_present, 1)

        prob_new_vec <- rep(0, N)
        for (j in seq(N)) {
            time_point <- times[j]
            allele_index <- alleles[j]

            p_a_present = 1 - prod((1 - prob_allele_in_cluster[allele_index,]) *
                                       prob_cluster_in_time[time_point,] +
                                       (1 - prob_cluster_in_time[time_point,]))

            p_a_not_in_c_at_previous_times <- vector(length = C)
            for (cluster in 1:C) {
                previously_occurring_times <-
                    times[alleles == allele_index & times < time_point & y == 1]
                if (length(previously_occurring_times) > 0) {
                    p_a_not_in_c_at_previous_times[cluster] <-
                        prod(1 - p_a_in_c_at_t_given_a_present[
                            allele_index, cluster, previously_occurring_times])
                } else {
                    p_a_not_in_c_at_previous_times[cluster] <- 1
                }
            }

            p_a_present_and_new =
                1 - prod((1 - p_a_not_in_c_at_previous_times *
                              prob_allele_in_cluster[allele_index,]) *
                             prob_cluster_in_time[time_point,] +
                             (1 - prob_cluster_in_time[time_point,]))

            p_a_new_given_a_present <- p_a_present_and_new / p_a_present

            prob_new_vec[j] <- p_a_new_given_a_present
        }

        dataset$probabilities[dataset$allele %in% lc$nodeclusters$node &
                                  dataset$subject == subject_name] <-
            prob_new_vec
    }

    dataset$probabilities <-
        ifelse(dataset$present == 1, dataset$probabilities, NA)

    dataset <- dataset %>%
        dplyr::arrange(.data$original_row_ordering) %>%
        dplyr::select(-.data$original_row_ordering)

    return(list(probability_new = dataset$probabilities, fit = NULL))
}

#' add_probability_new
#'
#' Add a column with the probability an allele is new by aggregating a
#' probability matrix.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present`
#' @param probability_mat A matrix with one column per result of
#' `determine_probabilities_*` on an imputed dataset.
#' @param column_name What to name the new column
#' @return The dataset with a new column given by `column_name`
#' giving the probability the allele was new if present.
#'
#' @examples
#' library(foreach)
#' library(doParallel)
#' dataset_in <- data.frame(allele = c('A', 'A', 'A', NA, NA, 'B', NA, 'B'),
#'     subject = rep('A', 8),
#'     time = c(1, 2, 3, 4, 5, 6, 7, 8))
#'
#' qpcr_times <- data.frame(subject = rep('A', 1), time = c(7))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_qpcr_times(dataset, qpcr_times)
#'
#' n_imputations <- 10
#' imputed_mat <- impute_dataset(dataset, n_imputations)
#'
#' probabilities_simple_mat <-
#'      foreach(i = 1:n_imputations,
#'          .combine = cbind,
#'          .packages = c('dinemites', 'dplyr')) %do% {
#'          dataset_tmp <- dataset
#'          dataset_tmp$present <- imputed_mat[,i]
#'          probabilities_simple <- determine_probabilities_simple(dataset_tmp)
#'          probabilities_simple$probability_new
#'      }
#'
#' dataset <- add_probability_new(dataset, probabilities_simple_mat)
#'
#' @import dplyr
add_probability_new <- function(dataset,
                                probability_mat,
                                column_name = "probability_new") {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    if (any(c(column_name, 'probability_new_tmp') %in% colnames(dataset))) {
        stop(paste0("probability_new_tmp and ", column_name, " should not be ",
        "columns of the input dataset"))
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    if (nrow(dataset) != nrow(probability_mat)) {
        stop(paste0("number of rows in dataset not equal to ",
                    "number of rows in probability_mat"))
    }

    dataset <- dataset %>%
        dplyr::mutate(
            probability_new_tmp = rowMeans(probability_mat, na.rm = T))

    if (any(is.na(dataset$probability_new_tmp)[dataset$present == 1])) {
        stop("probability_new is NA at an allele that was present")
    }

    colnames(dataset)[colnames(dataset) == 'probability_new_tmp'] <-
        column_name

    return(dataset)
}




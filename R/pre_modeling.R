# Fill in zeros for a dataset that has only infections to start with
# Assuming dataset has subject, time, allele with NA allele if no infection

#' fill_in_dataset
#'
#' @export
#' @param dataset dataset
#' @return return
#' @import dplyr
fill_in_dataset <- function(dataset) {
    if (any(!c('subject', 'time', 'allele') %in% colnames(dataset))) {
        stop("dataset must contain subject, time, allele as columns")
    }

    expanded_df <- dataset %>%
        dplyr::group_by(subject) %>%
        tidyr::expand(allele = unique(dataset$allele[!is.na(dataset$allele)]),
                      time = unique(time),
                      present = 0)

    unique_dataset <- dataset %>%
        dplyr::distinct(subject, time, .keep_all = TRUE) %>%
        dplyr::select(-allele)

    expanded_df <- expanded_df %>%
        dplyr::left_join(unique_dataset, by = c("subject", "time"))

    dataset$present <- ifelse(is.na(dataset$allele), 0, 1)
    dataset <- dataset %>%
        dplyr::filter(!is.na(allele))

    dataset <- rbind(dataset, expanded_df)
    dataset <- dataset[!duplicated(paste0(dataset$allele, '-', dataset$time, '-', dataset$subject)),]
    return(dataset)
}

#' add_present_infection
#'
#' @export
#' @param dataset dataset
#' @return return
#' @import dplyr
add_present_infection <- function(dataset) {
    dataset <- dataset %>%
        dplyr::group_by(subject, time) %>%
        dplyr::mutate(present_infection = ifelse(any(present != 0), 1, 0)) %>%
        dplyr::ungroup()

    return(dataset)
}

#' add_persistent_column
#'
#' @export
#' @param dataset dataset
#' @return return
#' @import dplyr
add_persistent_column <- function(dataset) {
    dataset <- dataset %>%
        dplyr::group_by(subject, allele) %>%
        dplyr::mutate(persistent = ifelse(time > min(time[present == 1], Inf), 1, 0)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(time, subject, allele)

    return(dataset)
}

#' add_persistent_infection
#'
#' @export
#' @param dataset dataset
#' @return return
#' @import dplyr
add_persistent_infection <- function(dataset) {
    dataset <- dataset %>%
        dplyr::group_by(subject, allele) %>%
        dplyr::mutate(persistent_infection = ifelse(time > min(time[present_infection == 1], Inf), 1, 0)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(time, subject, allele)

    return(dataset)
}

#' add_lag_column
#'
#' @export
#' @param dataset dataset
#' @param lag_time lag_time
#' @return return
#' @import dplyr
add_lag_column <- function(dataset, lag_time = 30) {
    dataset <- dataset %>%
        dplyr::arrange(subject, time) %>%
        dplyr::group_by(subject, allele) %>%
        dplyr::mutate(lag = ifelse(sapply(time, function(t) {
            any(present[time >= (t - lag_time) & time < t] == 1)
        }), 1, 0)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(time, subject, allele)

    colnames(dataset)[colnames(dataset) == 'lag'] <-
        paste0('lag_', lag_time, collapse = '')

    return(dataset)
}

#' add_lag_infection
#'
#' @export
#' @param dataset dataset
#' @param lag_time lag_time
#' @return return
#' @import dplyr
add_lag_infection <- function(dataset, lag_time = 30) {
    dataset <- dataset %>%
        dplyr::arrange(subject, time) %>%
        dplyr::group_by(subject, allele) %>%
        dplyr::mutate(lag_infection = ifelse(sapply(time, function(t) {
            any(present_infection[time >= (t - lag_time) & time < t] == 1)
        }), 1, 0)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(time, subject, allele)

    colnames(dataset)[colnames(dataset) == 'lag_infection'] <-
        paste0('lag_infection_', lag_time, collapse = '')

    return(dataset)
}

# After treatment but before next new infection, zero lag/persistent and add treatment covariate

#' add_treatment_column
#'
#' @export
#' @param dataset dataset
#' @param treatments treatments
#' @param treatment_acute_start treatment_acute_start
#' @param treatment_longitudinal_start treatment_longitudinal_start
#' @param verbose verbose
#' @return return
#' @import dplyr
add_treatment_column <- function(dataset,
                                    treatments,
                                    treatment_acute_start = 1,
                                    treatment_longitudinal_start = 10,
                                    verbose = TRUE) {
    dataset$treatment_acute <- 0
    dataset$treatment_longitudinal <- 0
    for (i in seq(unique(dataset$subject))) {
        if (verbose & i %% 10 == 0) {
            message("Calculating treatment at the allele level for subject: ", i, "/", length(unique(dataset$subject)))
        }
        subject_current <- unique(dataset$subject)[i]
        if (!subject_current %in% treatments$subject) {
            next
        }
        dataset_tmp_1 <- dataset[dataset$subject == subject_current,]
        treatments_in_effect <- treatments$time[treatments$subject == subject_current]

        for (j in seq(unique(dataset_tmp_1$allele))) {
            allele_current <- unique(dataset_tmp_1$allele)[j]
            dataset_tmp_2 <- dataset_tmp_1[dataset_tmp_1$allele == allele_current,]

            treatment_acute <- logical(length(dataset_tmp_2$time))
            treatment_longitudinal <- logical(length(dataset_tmp_2$time))

            if (all(dataset_tmp_2$present == 0)) {
                dataset_tmp_2$treatment_acute <- 0
                dataset_tmp_2$treatment_longitudinal <- 0
            } else {
                for (time_tmp in seq_along(dataset_tmp_2$time)) {
                    if (any(dataset_tmp_2$present[dataset_tmp_2$time < dataset_tmp_2$time[time_tmp]] == 1)) {
                        current_time <- dataset_tmp_2$time[time_tmp]
                        recent_treatment <- max(treatments_in_effect[treatments_in_effect < current_time], -Inf)

                        # If present before within the last [treatment_acute_start, treatment_longitudinal_start)
                        treatment_acute[time_tmp] <- ifelse(current_time - recent_treatment >= treatment_acute_start &
                                                         current_time - recent_treatment < treatment_longitudinal_start, 1, 0)

                        # If present before, more than treatment_longitudinal_start ago, and no other present since
                        # the last treatment + treatment_longitudinal_start
                        recent_present <- max(dataset_tmp_2$time[dataset_tmp_2$present_infection == 1 &
                                                                     dataset_tmp_2$time < current_time], -Inf)
                        treatment_longitudinal[time_tmp] <- ifelse(current_time - recent_treatment >= treatment_longitudinal_start &
                                                                recent_present < recent_treatment + treatment_longitudinal_start, 1, 0)
                    } else {
                        treatment_acute[time_tmp] <- 0
                        treatment_longitudinal[time_tmp] <- 0
                    }
                }
                dataset_tmp_2$treatment_acute <- treatment_acute
                dataset_tmp_2$treatment_longitudinal <- treatment_longitudinal
            }

            dataset_tmp_1[dataset_tmp_1$allele == allele_current,] <- dataset_tmp_2
        }

        dataset_tmp_1[dataset_tmp_1$treatment_acute == 1 | dataset_tmp_1$treatment_longitudinal == 1,
                      grepl("^lag_", colnames(dataset_tmp_1)) | "persistent" == colnames(dataset_tmp_1)] <- 0
        dataset[dataset$subject == subject_current,] <- dataset_tmp_1
    }
    return(dataset)
}

#' add_treatment_infection
#'
#' @export
#' @param dataset dataset
#' @param treatments treatments
#' @param treatment_acute_start treatment_acute_start
#' @param treatment_longitudinal_start treatment_longitudinal_start
#' @param verbose verbose
#' @return return
#' @import dplyr
add_treatment_infection <- function(dataset,
                                    treatments,
                                    treatment_acute_start = 1,
                                    treatment_longitudinal_start = 10,
                                    verbose = TRUE) {
    dataset$treatment_acute_infection <- 0
    dataset$treatment_longitudinal_infection <- 0
    for (i in seq(unique(dataset$subject))) {
        if (verbose & i %% 10 == 0) {
            message("Calculating treatment at the infection level for subject: ", i, "/", length(unique(dataset$subject)))
        }
        subject_current <- unique(dataset$subject)[i]
        if (!subject_current %in% treatments$subject) {
            next
        }
        dataset_tmp_1 <- dataset[dataset$subject == subject_current,]
        treatment_acute <- logical(length(dataset_tmp_1$time))
        treatment_longitudinal <- logical(length(dataset_tmp_1$time))
        treatments_in_effect <- treatments$time[treatments$subject == subject_current]

        if (all(dataset_tmp_1$present_infection == 0)) {
            dataset_tmp_1$treatment_acute_infection <- 0
            dataset_tmp_1$treatment_longitudinal_infection <- 0
        } else {
            for (time_tmp in seq_along(dataset_tmp_1$time)) {
                if (any(dataset_tmp_1$present_infection[dataset_tmp_1$time < dataset_tmp_1$time[time_tmp]] == 1)) {
                    current_time <- dataset_tmp_1$time[time_tmp]

                    recent_treatment <- max(treatments_in_effect[treatments_in_effect < current_time], -Inf)
                    # If present before within the last [treatment_acute_start, treatment_longitudinal_start)
                    treatment_acute[time_tmp] <- ifelse(current_time - recent_treatment >= treatment_acute_start &
                                                     current_time - recent_treatment < treatment_longitudinal_start, 1, 0)

                    # If present before, more than treatment_longitudinal_start ago, and no other present since
                    # the last treatment + treatment_longitudinal_start
                    recent_present <- max(dataset_tmp_1$time[dataset_tmp_1$present_infection == 1 &
                                                                 dataset_tmp_1$time < current_time], -Inf)
                    treatment_longitudinal[time_tmp] <- ifelse(current_time - recent_treatment >= treatment_longitudinal_start &
                                                            recent_present < recent_treatment + treatment_longitudinal_start, 1, 0)
                } else {
                    treatment_acute[time_tmp] <- 0
                    treatment_longitudinal[time_tmp] <- 0
                }
            }
            dataset_tmp_1$treatment_acute_infection <- treatment_acute
            dataset_tmp_1$treatment_longitudinal_infection <- treatment_longitudinal
        }
        dataset_tmp_1[dataset_tmp_1$treatment_acute_infection == 1 | dataset_tmp_1$treatment_longitudinal_infection == 1,
                      grepl("^lag_infection_", colnames(dataset_tmp_1)) | grepl("^persistent_infection", colnames(dataset_tmp_1))] <- 0
        dataset[dataset$subject == subject_current,] <- dataset_tmp_1
    }
    return(dataset)
}

#' estimate_drop_out
#'
#' @export
#' @param dataset dataset
#' @param allowed_ordered_time_difference allowed_ordered_time_difference
#' @return return
#' @import dplyr
estimate_drop_out <- function(dataset, allowed_ordered_time_difference = 6) {
    gap_df <- dataset %>%
        dplyr::group_by(subject, allele) %>%
        dplyr::filter(sum(present == 1) >= 2) %>%
        dplyr::mutate(time = rank(time)) %>% # Convert to time sequence instead
        dplyr::mutate(
            min_time = min(time[present == 1], na.rm = TRUE),
            max_time = max(time[present == 1], na.rm = TRUE)
        ) %>%
        dplyr::filter(time > min_time & time < max_time,
               present %in% c(0,1),
               max_time - min_time < allowed_ordered_time_difference) %>%
        dplyr::summarise(
            proportion_present_in_range = mean(present == 1),
            n = n()
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(subject) %>%
        dplyr::summarise(proportion_present_in_range =
                      sum(proportion_present_in_range * n) / sum(n),
                  n = sum(n))

    return(1 - sum(gap_df$proportion_present_in_range * gap_df$n) /
        sum(gap_df$n))
}

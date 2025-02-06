check_alleles_unique_across_loci <- function(dataset) {
    if ('locus' %in% colnames(dataset)) {
        if (any((dataset %>%
                 dplyr::filter(!is.na(.data$allele)) %>%
                 dplyr::group_by(.data$allele) %>%
                 dplyr::summarize(n_loci = length(unique(.data$locus)),
                                  .groups = "drop"))$n_loci
                > 1)) {
            stop(paste0("alleles at different loci must not share names,",
                        " consider appending the locus name to the allele"))
        }
    }
}

check_alleles_same_across_subject_times <- function(dataset) {
    alleles_by_time <- dataset %>%
        dplyr::group_by(.data$time) %>%
        dplyr::summarise(alleles = list(sort(unique(.data$allele))),
                         .groups = "drop")

    consistent_across_time <- alleles_by_time %>%
        dplyr::summarise(all_same = length(unique(.data$alleles)) == 1) %>%
        dplyr::pull(.data$all_same)

    alleles_by_subject <- dataset %>%
        dplyr::group_by(.data$subject) %>%
        dplyr::summarise(alleles = list(sort(unique(.data$allele))), .groups = "drop")

    consistent_across_subjects <- alleles_by_subject %>%
        dplyr::summarise(all_same = length(unique(.data$alleles)) == 1) %>%
        dplyr::pull(.data$all_same)

    is_consistent <- consistent_across_time && consistent_across_subjects

    if (!is_consistent) {
        stop("All subject-times must have the same set of alleles")
    }
}

#' fill_in_dataset
#'
#' Fill in zeros for a dataset that contains only alleles that were present.
#'
#' @export
#' @param dataset A longitudinal dataset with columns `subject`, `time`,
#' `allele.`
#' One row should be present for each allele sequenced at the subject-time.
#' All other subject-times should be included with allele as `NA`. A `locus`
#' column can also be specified to give each allele's locus.
#' @return A data.frame with columns `subject`, `time`, `allele`, and `present`
#' with `present` as 1 if the allele is present and 0 otherwise.
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#'     subject = rep('A', 5),
#'     time = c(1, 1, 2, 3, 4))
#'
#' dataset <- fill_in_dataset(dataset_in)
#'
#' @import dplyr
fill_in_dataset <- function(dataset) {
    if (any(!c('subject', 'time', 'allele') %in% colnames(dataset))) {
        stop("dataset must contain subject, time, allele as columns")
    }

    if ('present' %in% colnames(dataset)) {
        stop("present should not be a column of the input dataset")
    }

    check_alleles_unique_across_loci(dataset)

    expanded_df <- dataset %>%
        dplyr::group_by(.data$subject) %>%
        tidyr::expand(allele = unique(dataset$allele[!is.na(dataset$allele)]),
                      time = unique(.data$time),
                      present = 0) %>%
        dplyr::ungroup()

    unique_dataset <- dataset %>%
        dplyr::distinct(.data$subject, .data$time, .keep_all = TRUE) %>%
        dplyr::select(-.data$allele)

    if ('locus' %in% colnames(dataset)) {
        expanded_df <- expanded_df %>%
            mutate(locus = plyr::mapvalues(expanded_df$allele,
                                           dataset$allele,
                                           dataset$locus,
                                           warn_missing = FALSE))

        unique_dataset <- unique_dataset %>%
            dplyr::select(-.data$locus)
    }

    expanded_df <- expanded_df %>%
        dplyr::left_join(unique_dataset, by = c("subject", "time"))

    dataset <- dataset %>%
        dplyr::ungroup() %>%
        dplyr::mutate(present = ifelse(is.na(.data$allele), 0, 1)) %>%
        dplyr::filter(!is.na(.data$allele))

    dataset <- rbind(dataset, expanded_df)
    dataset <- dataset[!duplicated(paste0(dataset$allele, '-',
                                          dataset$time, '-',
                                          dataset$subject)),]
    dataset <- dataset %>%
        dplyr::arrange(.data$time, .data$subject, .data$allele)

    return(dataset)
}

#' add_qpcr_times
#'
#' Add times that were only qPCR positive (without sequencing information) as
#' 2s in the `present` column. Times marked as qPCR positive that already
#' have sequencing information will be disregarded and left with their
#' sequencing values.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `subject`,
#' `time`, `allele`, and `present`. Times that are qPCR positive only should
#' already be included in the dataset with `present = 0`.
#' @param qpcr_times A data.frame of qPCR positive only samples with columns
#' `subject` and `time`
#' @return The `dataset` with `present` set to 2 for subject times in
#' `qpcr_times` if there was not already an `allele` with `present == 1`
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#' subject = rep('A', 5),
#' time = c(1, 1, 2, 3, 4))
#'
#' qpcr_times <- data.frame(subject = rep('A', 2),
#'                          time = c(1, 4))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_qpcr_times(dataset, qpcr_times)
#'
#' @import dplyr
add_qpcr_times <- function(dataset, qpcr_times) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    if (any(!c("subject", "time") %in%
            colnames(qpcr_times))) {
        stop("qpcr_times must contain the columns: subject, time")
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    qpcr_times <- qpcr_times %>%
        dplyr::group_by(.data$subject) %>%
        tidyr::expand(allele = unique(dataset$allele),
                      time = unique(.data$time),
                      present = 2) %>%
        dplyr::ungroup()

    if ('locus' %in% colnames(dataset)) {
        qpcr_times <- qpcr_times %>%
            mutate(locus = plyr::mapvalues(qpcr_times$allele,
                                           dataset$allele,
                                           dataset$locus,
                                           warn_missing = FALSE))
    }

    qpcr_times <- qpcr_times %>%
        dplyr::mutate(
            subject_time = as.character(interaction(.data$subject, .data$time)))

    dataset <- dataset %>%
        dplyr::mutate(
            subject_time = as.character(interaction(.data$subject, .data$time)))

    if (any(!qpcr_times$subject_time %in% dataset$subject_time)) {
        stop("All qpcr_times subject times must be in the dataset")
    }

    any_present_df <- dataset %>%
        dplyr::group_by(.data$subject_time) %>%
        dplyr::filter(any(.data$present == 1)) %>%
        dplyr::ungroup()

    qpcr_times <- qpcr_times %>%
        dplyr::filter(!.data$subject_time %in% any_present_df$subject_time)

    if ('locus' %in% colnames(dataset)) {
        dataset <- left_join(dataset, qpcr_times,
            by = c("time", "subject", "allele", "locus", "subject_time"))
    } else {
        dataset <- left_join(dataset, qpcr_times,
            by = c("time", "subject", "allele", "subject_time"))
    }
    dataset <- dataset %>%
        dplyr::mutate(
            present = pmax(.data$present.x, .data$present.y, na.rm = T)) %>%
        dplyr::select(
            -.data$present.x, -.data$present.y, -.data$subject_time) %>%
        dplyr::arrange(.data$time, .data$subject, .data$allele)

    return(dataset)
}

#' add_present_infection
#'
#' Add a column with an indicator of whether an infection is present at the
#' time.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present.`
#' @return The dataset with a new column `present_infection` (0/1) indicating
#' whether an infection is present at the time.
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#'     subject = rep('A', 5),
#'     time = c(1, 1, 2, 3, 4))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_present_infection(dataset)
#'
#' @import dplyr
add_present_infection <- function(dataset) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    if ('present_infection' %in% colnames(dataset)) {
        stop("present_infection should not be a column of the input dataset")
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    dataset <- dataset %>%
        dplyr::group_by(.data$subject, .data$time) %>%
        dplyr::mutate(present_infection =
                          ifelse(any(.data$present != 0), 1, 0)) %>%
        dplyr::ungroup()

    return(dataset)
}

#' add_persistent_column
#'
#' Add a column with an indicator of whether an infection with the allele
#' has occurred before.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present.`
#' @return The dataset with a new column `persistent` (0/1) indicating whether
#' the allele has been observed before.
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#'     subject = rep('A', 5),
#'     time = c(1, 1, 2, 3, 4))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_present_infection(dataset)
#' dataset <- add_persistent_column(dataset)
#'
#' @import dplyr
add_persistent_column <- function(dataset) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    if ('persistent' %in% colnames(dataset)) {
        stop("persistent should not be a column of the input dataset")
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    dataset <- dataset %>%
        dplyr::group_by(.data$subject, .data$allele) %>%
        dplyr::mutate(persistent =
            ifelse(.data$time >
                       min(.data$time[.data$present == 1], Inf), 1, 0)) %>%
        dplyr::ungroup()

    return(dataset)
}

#' add_persistent_infection
#'
#' Add a column with an indicator of whether an infection has occurred before.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present_infection`
#' @return The dataset with a new column `persistent_infection` (0/1)
#' indicating whether any infection has been observed before.
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#'     subject = rep('A', 5),
#'     time = c(1, 1, 2, 3, 4))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_present_infection(dataset)
#' dataset <- add_persistent_infection(dataset)
#'
#' @import dplyr
add_persistent_infection <- function(dataset) {
    if (any(!c("allele", "subject", "time", "present_infection") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present_infection")
    }

    if ('persistent_infection' %in% colnames(dataset)) {
        stop("persistent_infection should not be a column of the input dataset")
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    dataset <- dataset %>%
        dplyr::group_by(.data$subject, .data$allele) %>%
        dplyr::mutate(persistent_infection =
            ifelse(.data$time >
                       min(.data$time[.data$present_infection == 1], Inf),
                   1, 0)) %>%
        dplyr::ungroup()

    return(dataset)
}

#' add_lag_column
#'
#' Add a column with an indicator of whether an infection with the allele
#' has occurred in the last `lag_time` times.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present`
#' @param lag_time The time within which previous infections count towards the
#' indicator
#' @return The dataset with a new column `lag_[lag_time]` (0/1)
#' indicating whether the allele has been observed in the last `lag_time` times.
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#'     subject = rep('A', 5),
#'     time = c(1, 1, 14, 28, 42))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_present_infection(dataset)
#' dataset <- add_lag_column(dataset)
#'
#' @import dplyr
add_lag_column <- function(dataset, lag_time = 30) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    if (any(c('lag', paste0('lag_', lag_time, collapse = '')) %in%
            colnames(dataset))) {
        stop(paste0("lag and ", paste0('lag_', lag_time, collapse = ''),
            " should not be columns of the input dataset"))
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    dataset <- dataset %>%
        dplyr::mutate(original_row_ordering = seq(nrow(dataset))) %>%
        dplyr::arrange(.data$subject, .data$time) %>%
        dplyr::group_by(.data$subject, .data$allele) %>%
        dplyr::mutate(lag = ifelse(sapply(.data$time, function(t) {
            any(.data$present[.data$time >= (t - lag_time) &
                                  .data$time < t] == 1)
        }), 1, 0)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(.data$original_row_ordering) %>%
        dplyr::select(-.data$original_row_ordering)

    colnames(dataset)[colnames(dataset) == 'lag'] <-
        paste0('lag_', lag_time, collapse = '')

    return(dataset)
}

#' add_lag_infection
#'
#' Add a column with an indicator of whether any infection
#' has occurred in the last `lag_time` times.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present_infection`
#' @param lag_time The time within which previous infections count towards the
#' indicator
#' @return The dataset with a new column `lag_infection_[lag_time]` (0/1)
#' indicating whether any infection has been observed in the last `lag_time`
#' times.
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#'     subject = rep('A', 5),
#'     time = c(1, 1, 14, 28, 42))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_present_infection(dataset)
#' dataset <- add_lag_infection(dataset)
#'
#' @import dplyr
add_lag_infection <- function(dataset, lag_time = 30) {
    if (any(!c("allele", "subject", "time", "present_infection") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present_infection")
    }

    if (any(c('lag_infection', paste0('lag_infection_',
                                      lag_time, collapse = '')) %in%
            colnames(dataset))) {
        stop(paste0("lag_infection and ",
                    paste0('lag_infection_', lag_time, collapse = ''),
                    " should not be columns of the input dataset"))
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    dataset <- dataset %>%
        dplyr::mutate(original_row_ordering = seq(nrow(dataset))) %>%
        dplyr::arrange(.data$subject, .data$time) %>%
        dplyr::group_by(.data$subject, .data$allele) %>%
        dplyr::mutate(lag_infection = ifelse(sapply(.data$time, function(t) {
            any(.data$present_infection[.data$time >= (t - lag_time) &
                                      .data$time < t] == 1)
        }), 1, 0)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(.data$original_row_ordering) %>%
        dplyr::select(-.data$original_row_ordering)

    colnames(dataset)[colnames(dataset) == 'lag_infection'] <-
        paste0('lag_infection_', lag_time, collapse = '')

    return(dataset)
}

#' add_treatment_column
#'
#' Add columns with indicators of whether there is acute or longitudinal
#' (any previous) treatment for the allele.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present`
#' @param treatments A data.frame of treatments with columns `subject` and
#' `time`
#' @param treatment_acute_start Number of days after the treatment at which
#' to mark the start of the acute treatment phase, corresponding to when
#' variants are being cleared but might still be detected.
#' @param treatment_longitudinal_start Number of days after the treatment
#' at which to mark the start of the longitudinal treatment phase, corresponding
#' to when cleared variants will no longer be detected.
#' @param verbose Print an update every 10 subjects since this function can
#' take some time.
#' @return The dataset with new columns `treatment_acute` and
#' `treatment_longitudinal` (both 0/1) indicating that the infection
#' was treated between `treatment_acute_start` (inclusive) and
#' `treatment_longitudinal_start` (exclusive) times ago (acute) or at least
#' `treatment_longitudinal_start` times ago (longitudinal).
#' Columns for `persistent` and `lag_[#]` will be set to 0 for the allele
#' after treatment until the next infection.
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', NA, NA, NA, 'B'),
#'                          subject = rep('A', 5),
#'                          time = c(1, 29, 30, 39, 50))
#' treatments <- data.frame(subject = c('A'),
#'                          time = c(29))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_present_infection(dataset)
#' dataset <- add_persistent_column(dataset)
#' dataset <- add_lag_column(dataset)
#' dataset <- add_treatment_column(dataset, treatments)
#'
#' @import dplyr
add_treatment_column <- function(dataset,
                                    treatments,
                                    treatment_acute_start = 1,
                                    treatment_longitudinal_start = 10,
                                    verbose = TRUE) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    if (any(!c("subject", "time") %in%
            colnames(treatments))) {
        stop("treatments must contain the columns: subject, time")
    }

    if (!(any(grepl("^lag_", colnames(dataset)) |
          "persistent" == colnames(dataset)))) {
        warning('Columns for persistent and lag not found.
                Columns for treatment should be added after persistent and lag')
    }

    if (any(c('treatment_acute', 'treatment_longitudinal') %in%
            colnames(dataset))) {
        stop("treatment_acute and treatment_longitudinal
             should not be columns of the input dataset")
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    dataset <- dataset %>%
        dplyr::mutate(original_row_ordering = seq(nrow(dataset)))

    dataset$treatment_acute <- 0
    dataset$treatment_longitudinal <- 0
    for (i in seq(unique(dataset$subject))) {
        if (verbose & i %% 10 == 0) {
            message("Calculating treatment at the allele level for subject: ",
                    i, "/", length(unique(dataset$subject)))
        }
        subject_current <- unique(dataset$subject)[i]
        if (!subject_current %in% treatments$subject) {
            next
        }
        dataset_tmp_1 <- dataset %>%
            dplyr::filter(.data$subject == subject_current)
        treatments_in_effect <-
            treatments$time[treatments$subject == subject_current]

        for (j in seq(unique(dataset_tmp_1$allele))) {
            allele_current <- unique(dataset_tmp_1$allele)[j]
            dataset_tmp_2 <- dataset_tmp_1 %>%
                dplyr::filter(.data$allele == allele_current)

            treatment_acute <- logical(length(dataset_tmp_2$time))
            treatment_longitudinal <- logical(length(dataset_tmp_2$time))

            if (all(dataset_tmp_2$present == 0)) {
                dataset_tmp_2$treatment_acute <- 0
                dataset_tmp_2$treatment_longitudinal <- 0
            } else {
                for (time_tmp in seq_along(dataset_tmp_2$time)) {
                    if (any(dataset_tmp_2$present[dataset_tmp_2$time <
                                    dataset_tmp_2$time[time_tmp]] == 1)) {
                        current_time <- dataset_tmp_2$time[time_tmp]
                        recent_treatment <- max(treatments_in_effect[
                            treatments_in_effect < current_time], -Inf)

                        # If present before within the last
                        # [treatment_acute_start, treatment_longitudinal_start)
                        treatment_acute[time_tmp] <-
                            ifelse(current_time - recent_treatment >=
                                       treatment_acute_start &
                                current_time - recent_treatment <
                                    treatment_longitudinal_start, 1, 0)

                        # If present before, more than
                        # treatment_longitudinal_start ago, and no other
                        # present since the last
                        # treatment + treatment_longitudinal_start
                        recent_present <-
                            max(dataset_tmp_2$time[
                                dataset_tmp_2$present == 1 &
                                    dataset_tmp_2$time < current_time], -Inf)
                        treatment_longitudinal[time_tmp] <-
                            ifelse(current_time - recent_treatment >=
                                treatment_longitudinal_start &
                                recent_present <
                                recent_treatment + treatment_longitudinal_start,
                                1, 0)
                    } else {
                        treatment_acute[time_tmp] <- 0
                        treatment_longitudinal[time_tmp] <- 0
                    }
                }
                dataset_tmp_2$treatment_acute <- treatment_acute
                dataset_tmp_2$treatment_longitudinal <- treatment_longitudinal
            }

            dataset_tmp_1[dataset_tmp_1$allele == allele_current,] <-
                dataset_tmp_2
        }

        dataset_tmp_1[dataset_tmp_1$treatment_acute == 1 |
                          dataset_tmp_1$treatment_longitudinal == 1,
                      (grepl("^lag_", colnames(dataset_tmp_1)) &
                           !grepl("infection", colnames(dataset_tmp_1))) |
                          "persistent" == colnames(dataset_tmp_1)] <- 0
        dataset[dataset$subject == subject_current,] <- dataset_tmp_1
    }

    dataset <- dataset %>%
        dplyr::arrange(.data$original_row_ordering) %>%
        dplyr::select(-.data$original_row_ordering)

    return(dataset)
}

#' add_treatment_infection
#'
#' Add columns with indicators of whether there is acute or longitudinal
#' (any previous) treatment for any allele.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present_infection`
#' @param treatments A data.frame of treatments with columns `subject` and
#' `time`
#' @param treatment_acute_start Number of days after the treatment at which
#' to mark the start of the acute treatment phase, corresponding to when
#' variants are being cleared but might still be detected.
#' @param treatment_longitudinal_start Number of days after the treatment
#' at which to mark the start of the longitudinal treatment phase, corresponding
#' to when cleared variants will no longer be detected.
#' @param verbose Print an update every 10 subjects since this function can
#' take some time.
#' @return The dataset with new columns `treatment_acute_infection` and
#' `treatment_longitudinal_infection` (both 0/1) indicating that the infection
#' was treated between `treatment_acute_start` (inclusive) and
#' `treatment_longitudinal_start` (exclusive) times ago (acute) or at least
#' `treatment_longitudinal_start` times ago (longitudinal).
#' Columns for `persistent_infection` and `lag_infection_[#]` will be set to 0
#' after treatment until the next infection.
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', NA, NA, NA, 'B'),
#'                          subject = rep('A', 5),
#'                          time = c(1, 29, 30, 39, 50))
#' treatments <- data.frame(subject = c('A'),
#'                          time = c(29))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' dataset <- add_present_infection(dataset)
#' dataset <- add_persistent_column(dataset)
#' dataset <- add_lag_column(dataset)
#' dataset <- add_treatment_infection(dataset, treatments)
#'
#' @import dplyr
add_treatment_infection <- function(dataset,
                                    treatments,
                                    treatment_acute_start = 1,
                                    treatment_longitudinal_start = 10,
                                    verbose = TRUE) {
    if (any(!c("allele", "subject", "time", "present_infection") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present_infection")
    }

    if (any(!c("subject", "time") %in%
            colnames(treatments))) {
        stop("treatments must contain the columns: subject, time")
    }

    if (!any(grepl("^lag_infection_", colnames(dataset)) |
          "persistent_infection" == colnames(dataset))) {
        warning('Columns for persistent_infection and lag_infection not found.
                Columns for treatment should be added after
                persistent_infection and lag_infection')
    }

    if (any(c('treatment_acute_infection',
              'treatment_longitudinal_infection') %in% colnames(dataset))) {
        stop("treatment_acute_infection and treatment_longitudinal_infection
             should not be columns of the input dataset")
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    dataset <- dataset %>%
        dplyr::mutate(original_row_ordering = seq(nrow(dataset)))

    dataset$treatment_acute_infection <- 0
    dataset$treatment_longitudinal_infection <- 0
    for (i in seq(unique(dataset$subject))) {
        if (verbose & i %% 10 == 0) {
            message(
                "Calculating treatment at the infection level for subject: ",
                i, "/", length(unique(dataset$subject)))
        }
        subject_current <- unique(dataset$subject)[i]
        if (!subject_current %in% treatments$subject) {
            next
        }
        dataset_tmp_1 <- dataset[dataset$subject == subject_current,]
        treatment_acute <- logical(length(dataset_tmp_1$time))
        treatment_longitudinal <- logical(length(dataset_tmp_1$time))
        treatments_in_effect <-
            treatments$time[treatments$subject == subject_current]

        if (all(dataset_tmp_1$present_infection == 0)) {
            dataset_tmp_1$treatment_acute_infection <- 0
            dataset_tmp_1$treatment_longitudinal_infection <- 0
        } else {
            for (time_tmp in seq_along(dataset_tmp_1$time)) {
                if (any(dataset_tmp_1$present_infection[dataset_tmp_1$time <
                            dataset_tmp_1$time[time_tmp]] == 1)) {
                    current_time <- dataset_tmp_1$time[time_tmp]

                    recent_treatment <- max(treatments_in_effect[
                        treatments_in_effect < current_time], -Inf)
                    # If present before within the last
                    # [treatment_acute_start, treatment_longitudinal_start)
                    treatment_acute[time_tmp] <-
                        ifelse(current_time - recent_treatment >=
                                   treatment_acute_start &
                                current_time - recent_treatment <
                                   treatment_longitudinal_start, 1, 0)

                    # If present before, more than treatment_longitudinal_start
                    # ago, and no other present since
                    # the last treatment + treatment_longitudinal_start
                    recent_present <-
                        max(dataset_tmp_1$time[
                                dataset_tmp_1$present_infection == 1 &
                                dataset_tmp_1$time < current_time], -Inf)
                    treatment_longitudinal[time_tmp] <-
                        ifelse(current_time - recent_treatment >=
                                   treatment_longitudinal_start &
                               recent_treatment + treatment_longitudinal_start >
                                   recent_present, 1, 0)
                } else {
                    treatment_acute[time_tmp] <- 0
                    treatment_longitudinal[time_tmp] <- 0
                }
            }
            dataset_tmp_1$treatment_acute_infection <- treatment_acute
            dataset_tmp_1$treatment_longitudinal_infection <-
                treatment_longitudinal
        }
        dataset_tmp_1[
            dataset_tmp_1$treatment_acute_infection == 1 |
                dataset_tmp_1$treatment_longitudinal_infection == 1,
            grepl("^lag_infection_", colnames(dataset_tmp_1)) |
                grepl("^persistent_infection", colnames(dataset_tmp_1))] <- 0
        dataset[dataset$subject == subject_current,] <- dataset_tmp_1
    }

    dataset <- dataset %>%
        dplyr::arrange(.data$original_row_ordering) %>%
        dplyr::select(-.data$original_row_ordering)

    return(dataset)
}

#' estimate_drop_out
#'
#' Estimate the sequencing drop-out rate in a dataset by selecting alleles
#' with close times of first and last occurrence and checking how often the
#' allele is observed in between.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present`
#' @param gap_size Alleles with first and last occurrences separated by
#' at most `gap_size` will be included in the estimation
#' @return Numeric vector with the estimated drop-out rate, the number of
#' gap positions evaluated, and the number of alleles present in those gaps.
#'
#' @examples
#' dataset_in <- data.frame(allele = c('A', 'B', NA, 'A', 'A'),
#'     subject = rep('A', 5),
#'     time = c(1, 1, 14, 28, 42))
#'
#' dataset <- fill_in_dataset(dataset_in)
#' estimate_drop_out(dataset)
#'
#' @import dplyr
estimate_drop_out <- function(dataset, gap_size = 3) {
    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    gap_df <- dataset %>%
        dplyr::group_by(.data$subject, .data$allele) %>%
        dplyr::filter(sum(.data$present == 1) >= 2) %>%
        # Convert to time sequence instead
        dplyr::mutate(time = rank(.data$time)) %>%
        dplyr::mutate(
            min_time = min(.data$time[.data$present == 1], Inf, na.rm = TRUE),
            max_time = max(.data$time[.data$present == 1], -Inf, na.rm = TRUE)
        ) %>%
        dplyr::filter(.data$time > .data$min_time & .data$time < .data$max_time,
                      .data$present %in% c(0,1),
                      .data$max_time - .data$min_time <= gap_size) %>%
        dplyr::summarise(
            proportion_present_in_range = mean(.data$present == 1),
            n = n(), .groups = 'drop'
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$subject) %>%
        dplyr::summarise(proportion_present_in_range =
                      sum(.data$proportion_present_in_range * .data$n) /
                          sum(.data$n),
                      n = sum(.data$n),
                      .groups = 'drop')

    return_vec <- c("Drop out rate" =
                        1 - sum(gap_df$proportion_present_in_range * gap_df$n) /
                        sum(gap_df$n),
                    "Evaluated gaps" = sum(gap_df$n),
                    "Number present in gap" =
                        sum(gap_df$proportion_present_in_range * gap_df$n))

    return(return_vec)
}

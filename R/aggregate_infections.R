#' merge_probability_columns
#'
#' Take two columns of probabilities in the dataset and merge them by taking
#' the average if they are close or manually selecting when the probabilities
#' are discordant. Probabilities can be selected manually by typing "f" (first
#' probability), "s" (second probability), or a numeric value between 0 and 1.
#' This function should only be run in an interactive session or it will
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, `present`, and `probability_present` if using imputed
#' data. This dataset should also contain columns with names given by
#' `cols_to_merge`.
#' @param cols_to_merge Which columns of probabilities in `dataset` to merge.
#' @param threshold Probabilities within `threshold` will be averaged
#' rather than manually selected.
#' @return The `dataset` with the column `probability_new` given by
#' the average of the two columns when the two
#' columns are within `threshold` of each other and the manually chosen value
#' otherwise.
#'
merge_probability_columns <- function(dataset, cols_to_merge, threshold = 0.3) {
    if (length(cols_to_merge) != 2) {
        stop("cols_to_merge must have length 2")
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    prob_abs_diffs <-
        abs(dataset[,cols_to_merge[1]] - dataset[,cols_to_merge[2]])

    if (!all(is.na(dataset[,cols_to_merge[1]]) ==
             is.na(dataset[,cols_to_merge[2]]))) {
        warning(paste0("The two probability columns do not have",
                       " NAs in the same positions"))
    }

    message(paste0(sum(prob_abs_diffs > threshold & !is.na(prob_abs_diffs)),
                   " discordant probabilities to be chosen manually."))

    prob_merged <- (dataset[,cols_to_merge[1]] + dataset[,cols_to_merge[2]]) / 2

    # For plotting
    dataset$probability_new <- prob_merged
    dataset$probability_new[
        which(prob_abs_diffs > threshold & !is.na(prob_abs_diffs))] <- NA

    # Run manual selection
    for (index in which(prob_abs_diffs > threshold & !is.na(prob_abs_diffs))) {
        plot_out <-
            plot_single_subject(dataset$subject[index],
                                dataset,
                                highlight_index = index)
        print(plot_out)
        print(dataset[index,])
        read_val <- readline(prompt =
            paste0("Choose probability: ",
            cols_to_merge[1],
            " (",
            round(dataset[index,cols_to_merge[1]], 3),
            ") (f), ",
            cols_to_merge[2],
            " (",
            round(dataset[index,cols_to_merge[2]], 3),
            ") (s), ",
            "other (value 0-1): "))
        read_val <- dplyr::case_when(
            read_val == 'f' ~ dataset[index,cols_to_merge[1]],
            read_val == 's' ~ dataset[index,cols_to_merge[2]],
            TRUE ~ tryCatch(as.numeric(read_val),
                            warning = function(w) {NA_real_}))
        prob_merged[index] <- read_val
        counter <- 1
        while ((prob_merged[index] > 1 |
               prob_merged[index] < 0 |
               is.na(prob_merged[index])) &
               counter <= 3) {
            read_val <-
                readline(prompt =
                    "Choose again: first (f), second (s), other (value 0-1): ")
            read_val <- dplyr::case_when(
                read_val == 'f' ~ dataset[index,cols_to_merge[1]],
                read_val == 's' ~ dataset[index,cols_to_merge[2]],
                TRUE ~ tryCatch(as.numeric(read_val),
                                warning = function(w) {NA_real_}))
            prob_merged[index] <- read_val
            counter <- counter + 1
        }
        if (counter > 3) {
            prob_merged[index] <- NA
            warning(
                "valid probability was not entered in 3 tries, moving on...")
        }
    }

    dataset$probability_new <- prob_merged
    return(dataset)
}

#' compute_total_new_COI
#'
#' Compute the total new complexity of infection (COI, number of unique
#' genetic variants) deduplicated across loci for each subject by
#' either summing across the loci and taking the maximum or maximizing the
#' diversity on each day over the loci and then summing the maxima.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, `present`, `locus`, `probability_new`, and
#' `probability_present` if using imputed data.
#' @param method Whether to sum the per-locus COI across all time points for
#' each locus and then take the maximum (`sum_then_max`) or take the maximum
#' per-locus COI at each time point across the loci and then sum the maxima
#' (`max_then_sum`).
#' @return A dataframe with a `subject` column and a `new_COI` column.
#' @examples
#'
#' library(dplyr)
#' dataset_in <- data.frame(
#'     locus = rep(c(1, 2), each = 5),
#'     allele = c('A', 'B', 'B', NA, NA, 'A', 'A', 'B', NA, NA),
#'     subject = rep('A', 10),
#'     time = rep(c(1, 1, 5, 15, 44), 2)) %>%
#'     mutate(allele = interaction(allele, locus))
#'
#' dataset <- fill_in_dataset(dataset_in)
#'
#' dataset$probability_new <-
#'     determine_probabilities_simple(dataset)$probability_new
#'
#' compute_total_new_COI(dataset, method = 'sum_then_max')
#' compute_total_new_COI(dataset, method = 'max_then_sum')
#'
#' @import dplyr
#'
compute_total_new_COI <- function(dataset, method = 'sum_then_max') {
    if (any(!c("allele", "subject",
               "time", "present", "probability_new") %in%
            colnames(dataset))) {
        stop(paste0("dataset must contain the columns:",
             " allele, subject, time, present, probability_new"))
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    if (!'probability_present' %in% colnames(dataset)) {
        message(paste0("probability_present not in colnames(dataset), ",
                       "using present == 1"))
        dataset$probability_present <- ifelse(dataset$present == 1, 1, 0)
    }

    if (!method %in% c("sum_then_max", "max_then_sum")) {
        stop("method must be sum_then_max or max_then_sum")
    }

    if (any(is.na(dataset$probability_new) & dataset$probability_present > 0)) {
        stop("probability_new is NA at a row in which probability_present > 0")
    }

    if (!'locus' %in% colnames(dataset)) {
        new_COI_out <- dataset %>%
            dplyr::group_by(.data$subject) %>%
            dplyr::summarise(
                new_COI = sum(.data$probability_present *
                                    .data$probability_new, na.rm=T),
                .groups = "drop")
    } else if (method == 'sum_then_max') {
        new_COI_out <- dataset %>%
            dplyr::group_by(.data$subject, .data$locus, .data$time) %>%
            dplyr::summarise(
                sum_first = sum(.data$probability_present *
                                    .data$probability_new, na.rm=T),
                .groups = "drop") %>%
            dplyr::group_by(.data$subject, .data$locus) %>%
            dplyr::summarise(sum_second = sum(.data$sum_first, na.rm = T),
                             .groups = "drop") %>%
            dplyr::group_by(.data$subject) %>%
            dplyr::summarise(new_COI = max(.data$sum_second, 0, na.rm=T),
                             .groups = "drop")
    } else if (method == 'max_then_sum') {
        new_COI_out <- dataset %>%
            dplyr::group_by(.data$subject, .data$locus, .data$time) %>%
            dplyr::summarise(
                sum_first = sum(.data$probability_present *
                                    .data$probability_new, na.rm=T),
                .groups = "drop") %>%
            dplyr::group_by(.data$subject, .data$time) %>%
            dplyr::summarise(
                max_second = max(.data$sum_first, 0, na.rm = T),
                .groups = "drop") %>%
            dplyr::group_by(.data$subject) %>%
            dplyr::summarise(new_COI = sum(.data$max_second, na.rm=T),
                             .groups = "drop")
    }

    return(new_COI_out)
}

# Count new infections using the algorithm described in the manuscript
count_new_infections <- function(total_new,
                                 max_new,
                                 max_new_weighted,
                                 min_new,
                                 min_new_weighted,
                                 max_prob) {
    if (length(total_new) == 0) {
        return(0)
    }

    if (length(total_new) == 1) {
        return(max_new[1])
    }

    # Define peaks and troughs
    peak_trough_regions <- list()
    current_region <- 1

    if (length(total_new) > 2) {
        for (i in 2:(length(total_new) - 1)) {
            if (total_new[i] < total_new[i - 1] &&
                total_new[i] <= total_new[i + 1] &
                (max(total_new[current_region]) -
                 min(total_new[c(current_region, i)]) >= 1)) {
                peak_trough_regions <-
                    append(peak_trough_regions, list(current_region))
                current_region <- i
            } else {
                current_region <- c(current_region, i)
            }
        }
    }

    current_region <- c(current_region, length(total_new))
    peak_trough_regions <- append(peak_trough_regions, list(current_region))
    # End peaks and troughs

    count <- 0
    for (peak_trough_region in peak_trough_regions) {
        first_certain_index <- min(which(max_prob[peak_trough_region] == 1),
                                   length(peak_trough_region) + 1)
        if (first_certain_index <=
            length(peak_trough_region) & first_certain_index > 1) {
            max_new_used <-
                c(max_new_weighted[peak_trough_region][
                    1:(first_certain_index - 1)],
                  max_new[peak_trough_region][
                      first_certain_index:length(peak_trough_region)])
            min_new_used <-
                c(min_new_weighted[peak_trough_region][
                    1:(first_certain_index - 1)],
                  min_new[peak_trough_region][
                      first_certain_index:length(peak_trough_region)])
        } else {
            max_new_used <- max_new[peak_trough_region]
            min_new_used <- min_new[peak_trough_region]
        }

        if (length(min_new_used) != length(peak_trough_region)) {
            stop(paste0(length(min_new_used), " ", length(peak_trough_region)))
        }
        if (max(total_new[peak_trough_region]) >= 1) {
            max_new_index <- which.max(max_new_used)
            count <- count + max_new_used[max_new_index]
            count <- count + sum(min_new_used[-max_new_index])
        } else {
            count <- count + max(min_new_used)
        }
    }

    return(count)
}

#' estimate_new_infections
#'
#' Compute the total new infection events for each subject by identifying
#' peaks in the pan-locus new COI.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, `present`, `probability_new`, and
#' `probability_present` if using imputed data.
#' @return A data.frame with a `subject` column and a `new_infections` column.
#' @examples
#'
#' library(dplyr)
#' dataset_in <- data.frame(
#'     locus = rep(c(1, 2), each = 5),
#'     allele = c('A', 'B', 'B', NA, NA, 'A', 'A', 'B', NA, NA),
#'     subject = rep('A', 10),
#'     time = rep(c(1, 1, 5, 15, 44), 2)) %>%
#'     mutate(allele = interaction(allele, locus))
#'
#' dataset <- fill_in_dataset(dataset_in)
#'
#' dataset$probability_new <-
#'     determine_probabilities_simple(dataset)$probability_new
#'
#' estimate_new_infections(dataset)
#'
#' @import dplyr
#'
estimate_new_infections <- function(dataset) {
    if (any(!c("allele", "subject", "time", "present", "probability_new") %in%
            colnames(dataset))) {
        stop(paste0("dataset must contain the columns:",
             " allele, subject, time, present, probability_new"))
    }

    if (!'probability_present' %in% colnames(dataset)) {
        message(paste0("probability_present not in colnames(dataset),",
                       " using present == 1"))
        dataset$probability_present <- ifelse(dataset$present == 1, 1, 0)
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    if (any(is.na(dataset$probability_new) & dataset$probability_present > 0)) {
        stop("probability_new is NA at a row in which probability_present > 0")
    }

    new_infections <- dataset %>%
        dplyr::group_by(.data$subject, .data$time) %>%
        dplyr::summarise(
            total_new = sum((.data$probability_present * .data$probability_new)[
                .data$probability_present > 0]),
             max_new = ifelse(
                 max(.data$probability_present) > 0.1,
                 max(.data$probability_new[.data$probability_present > 0.1], 0),
                 max(.data$probability_new[.data$probability_present > 0], 0)),
             max_new_weighted =
                ifelse(max(.data$probability_present) > 0.1,
                   max(.data$probability_new[.data$probability_present > 0.1] *
                           .data$probability_present[
                               .data$probability_present > 0.1], 0),
                   max(.data$probability_new[
                       .data$probability_present > 0], 0)),
             min_new =
                ifelse(max(.data$probability_present) > 0.1,
                    max(min(.data$probability_new[
                        .data$probability_present > 0.1], 2), 0),
                    max(min(.data$probability_new[
                        .data$probability_present > 0], 2), 0)),
             min_new_weighted =
                ifelse(max(.data$probability_present) > 0.1,
                       max(min(.data$probability_new[
                           .data$probability_present > 0.1] *
                               .data$probability_present[
                                   .data$probability_present > 0.1], 2), 0),
                       max(min(.data$probability_new[
                           .data$probability_present > 0], 2), 0)),
             max_prob = max(.data$probability_present),
            .groups = "drop") %>%
        dplyr::mutate(min_new = ifelse(.data$min_new == 2, 0, .data$min_new)) %>%
        dplyr::mutate(
            min_new_weighted = ifelse(.data$min_new_weighted == 2,
                                      0, .data$min_new_weighted)) %>%
        dplyr::mutate(
            max_new = ifelse(is.na(.data$max_new), 0, .data$max_new),
            max_new_weighted = ifelse(is.na(.data$max_new_weighted),
                                      0, .data$max_new_weighted),
            min_new = ifelse(is.na(.data$min_new), 0, .data$min_new),
            min_new_weighted =
                ifelse(is.na(.data$min_new_weighted),
                       0, .data$min_new_weighted),
            max_prob = ifelse(is.na(.data$max_prob), 0, .data$max_prob)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$subject) %>%
        dplyr::arrange(.data$time) %>%
        dplyr::summarise(
            new_infections = count_new_infections(
                .data$total_new[!is.na(.data$total_new)],
                .data$max_new[!is.na(.data$total_new)],
                .data$max_new_weighted[!is.na(.data$total_new)],
                .data$min_new[!is.na(.data$total_new)],
                .data$min_new_weighted[!is.na(.data$total_new)],
                .data$max_prob[!is.na(.data$total_new)]),
            .groups = "drop")

    return(new_infections)
}

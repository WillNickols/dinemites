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

    dataset <- dataset %>%
        dplyr::mutate(original_row_ordering = seq(nrow(dataset))) %>%
        dplyr::arrange(.data$subject, .data$allele, .data$time)

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

    dataset <- dataset %>%
        dplyr::arrange(.data$original_row_ordering) %>%
        dplyr::select(-.data$original_row_ordering)

    return(dataset)
}

#' compute_molFOI
#'
#' Compute the total molecular force of infection (molFOI, number of new unique
#' genetic variants) deduplicated across loci for each subject by
#' either summing across the loci and taking the maximum or maximizing the
#' diversity on each day over the loci and then summing the maxima.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, `present`, `locus`, `probability_new`, and
#' `probability_present` if using imputed data.
#' @param method Whether to sum the per-locus molFOI across all time points for
#' each locus and then take the maximum (`sum_then_max`) or take the maximum
#' per-locus new complexity at each time point across the loci and then sum
#' the maxima
#' (`max_then_sum`).
#' @return A dataframe with a `subject` column and a `molFOI` column.
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
#' compute_molFOI(dataset, method = 'sum_then_max')
#' compute_molFOI(dataset, method = 'max_then_sum')
#'
#' @import dplyr
#'
compute_molFOI <- function(dataset, method = 'sum_then_max') {
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
        compute_molFOI_out <- dataset %>%
            dplyr::group_by(.data$subject) %>%
            dplyr::summarise(
                molFOI = sum(.data$probability_present *
                                    .data$probability_new, na.rm=T),
                .groups = "drop")
    } else if (method == 'sum_then_max') {
        compute_molFOI_out <- dataset %>%
            dplyr::group_by(.data$subject, .data$locus, .data$time) %>%
            dplyr::summarise(
                sum_first = sum(.data$probability_present *
                                    .data$probability_new, na.rm=T),
                .groups = "drop") %>%
            dplyr::group_by(.data$subject, .data$locus) %>%
            dplyr::summarise(sum_second = sum(.data$sum_first, na.rm = T),
                             .groups = "drop") %>%
            dplyr::group_by(.data$subject) %>%
            dplyr::summarise(molFOI = max(.data$sum_second, 0, na.rm=T),
                             .groups = "drop")
    } else if (method == 'max_then_sum') {
        compute_molFOI_out <- dataset %>%
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
            dplyr::summarise(molFOI = sum(.data$max_second, na.rm=T),
                             .groups = "drop")
    }

    return(compute_molFOI_out)
}

# Count new infections using the algorithm described in the manuscript
count_new_infections <- function(total_new, max_new, min_new) {
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
        max_new_used <- max_new[peak_trough_region]
        min_new_used <- min_new[peak_trough_region]

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
#' peaks in the pan-locus new alleles.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, `present` (0/1), and `probability_new`.
#' @param imputation_mat If using imputations, the output of the function
#' `impute_dataset`: a
#' matrix with an equal number of rows to `dataset` and one column per
#' imputed dataset.
#' @param probability_mat If using imputations, a matrix with the same
#' dimensions as `imputation_mat`
#' with one column per imputed dataset containing the probability the allele
#' was new if present.
#' @return A data frame with row names corresponding to the subjects and
#' one column per imputed dataset (or a single column if not using
#' imputations).
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
estimate_new_infections <- function(dataset,
                                    imputation_mat = NULL,
                                    probability_mat = NULL) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop(paste0("dataset must contain the columns:",
             " allele, subject, time, present, probability_new"))
    }

    check_alleles_unique_across_loci(dataset)
    check_alleles_same_across_subject_times(dataset)

    if (any(dataset$present == 2) &
        (is.null(imputation_mat) | is.null(probability_mat))) {
        stop(paste0("imputation_mat and probability_mat must be supplied when",
                    " dataset$present has any 2s"))
    }

    if (!any(dataset$present == 2) &
        (is.null(imputation_mat) | is.null(probability_mat))) {
        # No imputations used
        if (!'probability_new' %in% colnames(dataset)) {
            stop(paste0("probability_new must be a column in dataset if not",
            " using imputed values"))
        }

        if (any(is.na(dataset$probability_new) & dataset$present == 1)) {
            stop("probability_new is NA at a row in which present == 1")
        }

        new_infections <- dataset %>%
            dplyr::group_by(.data$subject, .data$time) %>%
            dplyr::summarise(
                total_new = sum(.data$probability_new[.data$present == 1]),
                max_new = max(.data$probability_new[.data$present == 1], 0),
                min_new = max(min(
                    .data$probability_new[.data$present == 1], 2), 0),
                .groups = "drop") %>%
            dplyr::mutate(
                min_new = ifelse(.data$min_new == 2, 0, .data$min_new)) %>%
            dplyr::mutate(
                max_new = ifelse(is.na(.data$max_new), 0, .data$max_new),
                min_new = ifelse(is.na(.data$min_new), 0, .data$min_new)) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(.data$subject) %>%
            dplyr::arrange(.data$time) %>%
            dplyr::summarise(
                new_infections = count_new_infections(
                    .data$total_new[!is.na(.data$total_new)],
                    .data$max_new[!is.na(.data$total_new)],
                    .data$min_new[!is.na(.data$total_new)]),
                .groups = "drop") %>%
            as.data.frame()

        rownames(new_infections) <- new_infections$subject
        new_infections$subject <- NULL
    } else {
        # Imputations used
        if (any(dim(imputation_mat) != dim(probability_mat))) {
            stop(paste0("imputation_mat and probability_mat must ",
            "have the same dimensions"))
        }

        if (nrow(imputation_mat) != nrow(dataset)) {
            stop(paste0("imputation_mat and dataset must have the same",
                        " number of rows"))
        }

        if (any(c('present_tmp', 'probability_new_tmp') %in% colnames(dataset))) {
            stop(paste0("present_tmp and probability_new_tmp should not be",
            " columns in dataset if using multiple imputaiton"))
        }

        n_imputations <- ncol(imputation_mat)

        new_infections_mat <- matrix(ncol = n_imputations,
                                     nrow = length(unique(dataset$subject)))
        for (i in seq(n_imputations)) {
            dataset$present_tmp <- imputation_mat[,i]
            dataset$probability_new_tmp <- probability_mat[,i]
            if (any(dataset$present == 1 & dataset$present_tmp == 0) |
                any(dataset$present == 0 & dataset$present_tmp == 1)) {
                stop(paste0("Incompatible present column of dataset and",
                " imputation matrix column. Check that dataset has rows",
                " in the same order as when the imputations were created"))
            }

            if (any(is.na(dataset$probability_new_tmp) & dataset$present_tmp == 1)) {
                stop("probability_mat is NA at a position in which imputation_mat == 1")
            }

            new_infections <- dataset %>%
                dplyr::group_by(.data$subject, .data$time) %>%
                dplyr::summarise(
                    total_new = sum(.data$probability_new_tmp[.data$present_tmp == 1]),
                    max_new = max(.data$probability_new_tmp[.data$present_tmp == 1], 0),
                    min_new = max(min(
                        .data$probability_new_tmp[.data$present_tmp == 1], 2), 0),
                    .groups = "drop") %>%
                dplyr::mutate(
                    min_new = ifelse(.data$min_new == 2, 0, .data$min_new)) %>%
                dplyr::mutate(
                    max_new = ifelse(is.na(.data$max_new), 0, .data$max_new),
                    min_new = ifelse(is.na(.data$min_new), 0, .data$min_new)) %>%
                dplyr::ungroup() %>%
                dplyr::group_by(.data$subject) %>%
                dplyr::arrange(.data$time) %>%
                dplyr::summarise(
                    new_infections = count_new_infections(
                        .data$total_new[!is.na(.data$total_new)],
                        .data$max_new[!is.na(.data$total_new)],
                        .data$min_new[!is.na(.data$total_new)]),
                    .groups = "drop")

            new_infections_mat[,i] <- new_infections$new_infections
            rownames(new_infections_mat) <- new_infections$subject
        }
        new_infections <- data.frame(new_infections_mat)
    }

    if (ncol(new_infections) == 1) {
        colnames(new_infections) <- 'new_infections'
    } else {
        colnames(new_infections) <-
            paste0('imputation_', seq(ncol(new_infections)))
    }

    return(new_infections)
}

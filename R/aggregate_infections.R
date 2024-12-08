#' compute_total_new_molFOI
#'
#' @export
#' @param dataset dataset
#' @param method method
#' @return return
#' @import dplyr
#'
compute_total_new_molFOI <- function(dataset, method = 'sum_then_max') {
    if (any(!c("allele", "locus", "subject", "time", "present", "probability_new") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, locus, subject, time, present, probability_new")
    }

    if (!'probability_present' %in% colnames(dataset)) {
        warning("probability_present not in colnames(dataset), using present == 1")
        dataset$probability_present <- ifelse(dataset$present == 1, 1, 0)
    }

    if (!method %in% c("sum_then_max", "max_then_sum")) {
        stop("method must be sum_then_max or max_then_sum")
    }

    if (method == 'sum_then_max') {
        new_molFOI_out <- dataset %>%
            dplyr::group_by(subject, locus, time) %>%
            dplyr::summarise(sum_first = sum(probability_present * probability_new, na.rm=T)) %>%
            dplyr::group_by(subject, locus) %>%
            dplyr::summarise(sum_second = sum(sum_first, na.rm = T)) %>%
            dplyr::group_by(subject) %>%
            dplyr::summarise(new_molFOI = max(sum_second, 0, na.rm=T))
    } else if (method == 'max_then_sum') {
        new_molFOI_out <- dataset %>%
            dplyr::group_by(subject, locus, time) %>%
            dplyr::summarise(sum_first = sum(probability_present * probability_new, na.rm=T)) %>%
            dplyr::group_by(subject, time) %>%
            dplyr::summarise(max_second = max(sum_first, 0, na.rm = T)) %>%
            dplyr::group_by(subject) %>%
            dplyr::summarise(new_molFOI = sum(max_second, na.rm=T))
    }

    return(new_molFOI_out)
}

count_local_maxima <- function(total_new, max_new) {
    if (length(total_new) < 2) {
        return(min(max_new[1], 1))
    }

    count <- 0
    # Check first
    if (total_new[1] > total_new[2]) {
        count <- count + min(max_new[1], 1)
    }

    # Check last
    if (total_new[length(total_new)] > total_new[length(total_new) - 1]) {
        count <- count + min(max_new[length(total_new)], 1)
    }

    for (i in 2:(length(total_new) - 1)) {
        if (total_new[i] > total_new[i-1] && total_new[i] > total_new[i+1]) {
            count <- count + min(max_new[i], 1)
        }
    }

    return(count)
}

#' estimate_new_infections
#'
#' @export
#' @param dataset dataset
#' @return return
#' @import dplyr
#'
estimate_new_infections <- function(dataset) {
    if (any(!c("allele", "subject", "time", "present", "probability_new") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present, probability_new")
    }

    if (!'probability_present' %in% colnames(dataset)) {
        warning("probability_present not in colnames(dataset), using present == 1")
        dataset$probability_present <- ifelse(dataset$present == 1, 1, 0)
    }

    new_infections <- dataset %>%
        dplyr::group_by(.data$subject, .data$time) %>%
        dplyr::summarise(total_new = sum((probability_present * probability_new)[probability_present > 0]),
                         max_new = ifelse(max(probability_present) > 0.1,
                                          max(probability_new[probability_present > 0.1], 0),
                                          max(probability_new[probability_present > 0], 0))) %>%
        dplyr::mutate(total_new = ifelse(is.na(total_new), 0, total_new),
                      max_new = ifelse(is.na(max_new), 0, max_new)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$subject) %>%
        dplyr::arrange(.data$time) %>%
        dplyr::summarise(new_infections = count_local_maxima(.data$total_new[!is.na(.data$total_new)],
                                                                  .data$max_new[!is.na(.data$total_new)]))

    return(new_infections)
}

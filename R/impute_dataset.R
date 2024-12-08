#' impute_dataset
#'
#' @export
#' @param dataset dataset
#' @param n_imputations n_imputations
#' @param k k
#' @param n_cores n_cores
#' @param seed seed
#' @return return
#' @import dplyr
#' @import doParallel
#' @import foreach
#'
impute_dataset <- function(dataset,
                           n_imputations = 10,
                           k = NULL, # Even number such that k/2 points are included on either side
                           n_cores = 1,
                           seed = 1) {
    set.seed(seed)

    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    dataset_arranged <- dataset %>%
        arrange(subject, allele, time)

    if (any(dataset_arranged != dataset, na.rm = T)) {
        stop("dataset must be arranged by subject, allele, time")
    }

    if (is.null(k)) {
        k <- min(9, floor(mean(table(dataset$subject, dataset$allele))))
    }

    if (k < 1) {
        stop("k < 1")
    }

    # 2 if qPCR positive only
    if (any(!dataset$present %in% c(0,1,2))) {
        stop("present column must be 0/1/2")
    }

    if (!any(dataset$present == 2)) {
        stop("No need to impute if no present columns are '2'")
    }

    dataset <- dataset %>%
        arrange(time)

    dataset$time <- as.numeric(as.character(dataset$time))

    # Put 2s in all positions that could be replaced
    generate_plausible_2_vectors <- function(mat) {
        mat <- mat[apply(mat, 1, function(row) !any(row == 2)), , drop = FALSE]

        col_with_one <- apply(mat, 2, function(col) any(col == 1))
        if (all(!col_with_one)) {
            return(c())
        }

        row_counts <- table(apply(mat, 1, paste, collapse = ""))

        # Do the rest of this for each row_count name...

        indices_to_replace <- which(col_with_one)
        n_replacements <- length(indices_to_replace)

        # Generate all combinations of replacements
        possible_combinations <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n_replacements)))
        possible_combinations <- possible_combinations[rowSums(possible_combinations) > 0, , drop=F]

        original <- vector(length = nrow(possible_combinations) * length(row_counts))
        result <- vector(length = nrow(possible_combinations) * length(row_counts))
        result_tallies <- vector(length = nrow(possible_combinations) * length(row_counts))
        if (nrow(possible_combinations) > 0) {
            for (j in 1:length(row_counts)) {
                for (i in 1:nrow(possible_combinations)) {
                    new_vector <- unlist(strsplit(names(row_counts)[j], ''))
                    new_vector[indices_to_replace[possible_combinations[i, ]]] <- '2'
                    original[(j - 1) * nrow(possible_combinations) + i] <- names(row_counts)[j]
                    result[(j - 1) * nrow(possible_combinations) + i] <- paste0(new_vector, collapse = '')
                    result_tallies[(j - 1) * nrow(possible_combinations) + i] <- row_counts[j]
                }
            }
        }

        return(list(original = original, result = result, result_tallies = result_tallies))
    }

    # Get all matrices with k columns that are allowed for searching
    generate_k_column_submatrices <- function(mat, k) {
        n_cols <- ncol(mat)
        if (n_cols < k) {
            stop("Some matrices have fewer entries than k")
        }
        # Get all k-column submatrices
        submatrices <- lapply(1:(n_cols - k + 1), function(i) mat[,i:(i + k - 1), drop=F])

        submatrices <- Filter(function(x) !all(x != 1), submatrices)
        submatrices <- Filter(function(x) !any(x == 2), submatrices)

        return(submatrices)
    }

    comparison_table_list <- list()
    break_now <- FALSE
    for (k_sub in 1:(k+1)) {
        message(paste0("Generating matches for k = ", k_sub))
        k_column_submatrices_list <- list()
        k_column_submatrices_list_counter <- 1
        for (subject_cur in unique(dataset$subject)) {
            present_allele_mat <- data.frame(reshape2::dcast(dataset[dataset$subject == subject_cur,], allele ~ time,
                                                             value.var = "present", fun.aggregate = sum, fill = 0))
            present_allele_mat$allele <- NULL

            if (any(present_allele_mat == 1)) {
                # Generate all look-ups
                possibleError <- tryCatch({
                    k_column_submatrices <- generate_k_column_submatrices(present_allele_mat, k_sub)
                }, error = function(e) {
                    return(e)
                })

                if(inherits(possibleError, "error")) {
                    if (grepl("Some matrices have fewer entries than k", possibleError$message)) {
                        # cat("Error: ", possibleError$message, "\nSetting k to ", k_sub - 2, "\n")
                        # k <- k - 1
                        next
                    } else {
                        stop(possibleError$message)
                    }
                }

                k_column_submatrices_list[[k_column_submatrices_list_counter]] <- k_column_submatrices
                k_column_submatrices_list_counter <- k_column_submatrices_list_counter + 1
            }
        }

        # if (break_now) {
        #     break
        # }

        k_column_submatrices_list <- do.call(c, k_column_submatrices_list)


        unknown_list <- list()
        unknown_corresponding_list <- list()
        unknown_list_tallies <- list()
        counter <- 1
        for (submatrix in k_column_submatrices_list) {
            to_add <- generate_plausible_2_vectors(submatrix)
            unknown_list[[counter]] <- to_add$result
            unknown_corresponding_list[[counter]] <- to_add$original
            unknown_list_tallies[[counter]] <- to_add$result_tallies
            counter <- counter + 1
        }

        unknown_list <- unlist(unknown_list)
        unknown_corresponding_list <- unlist(unknown_corresponding_list)
        unknown_list_tallies <- unlist(unknown_list_tallies)

        if (!is.null(unknown_list_tallies)) {
            comparison_table <- tapply(unknown_list_tallies,
                                       list(unknown_list, unknown_corresponding_list),
                                       sum, default = 0)

            comparison_table_list[[k_sub]] <- comparison_table
        } else {
            comparison_table_list[[k_sub]] <- matrix(nrow = 0, ncol = 0)
        }
    }
    # Set k since it might be below originally set k
    k <- k_sub - 2

    dataset <- dataset %>%
        arrange(subject, allele, time)

    cat(paste0("Creating imputations"))
    cl <- parallel::makeCluster(n_cores)
    registerDoParallel(cl)

    imputation_mat <- matrix(nrow = nrow(dataset), ncol = n_imputations)
    results <- foreach(i = 1:n_imputations, .packages = c('dplyr', 'reshape2')) %dopar% {
        cat(paste0("Starting imputation: ", i, "\n"))

        # Initialize temporary storage for imputed values
        imputed_vals <- c()

        # Loop over each unique subject
        for (subject_cur in unique(dataset$subject)) {
            present_allele_mat <- data.frame(reshape2::dcast(dataset[dataset$subject == subject_cur,], allele ~ time,
                                                             value.var = "present", fun.aggregate = sum, fill = 0))

            present_allele_mat$allele <- NULL
            present_allele_mat_new <- present_allele_mat

            # Loop through each row in present_allele_mat
            for (row_num in 1:nrow(present_allele_mat)) {
                current_row <- unname(unlist(present_allele_mat[row_num,]))

                position_of_interest <- 1
                length_row <- length(current_row)

                while (any(current_row == 2)) {
                    if (current_row[position_of_interest] == 2) {
                        current_k <- k
                        while(current_k > 0) {
                            # Adjust sliding window around position of interest
                            half_window <- floor(current_k / 2)
                            left_bound <- max(1, position_of_interest - half_window)
                            right_bound <- min(length_row, position_of_interest + half_window)
                            sliding_window <- left_bound:right_bound

                            search_string <- paste0(current_row[sliding_window], collapse = '')

                            if(search_string %in% rownames(comparison_table_list[[length(sliding_window)]])) {
                                replacement <- sample(names(comparison_table_list[[length(sliding_window)]][search_string,]), 1,
                                                      prob = comparison_table_list[[length(sliding_window)]][search_string,] /
                                                          sum(comparison_table_list[[length(sliding_window)]][search_string,]))
                                current_row[sliding_window] <- as.numeric(unlist(strsplit(replacement, '')))
                                break
                            } else {
                                current_k <- current_k - 1
                            }
                        }
                        if (current_k == 0) {
                            stop("current_k == 0, something went wrong in the imputation procedure")
                        }
                    }
                    position_of_interest <- position_of_interest + 1
                    if (position_of_interest > length_row) break
                }

                present_allele_mat_new[row_num, ] <- current_row
            }

            # Add imputed values to growing list
            imputed_vals <- c(imputed_vals, c(t(as.matrix(present_allele_mat_new))))
        }

        # Return the imputed values for this iteration
        imputed_vals
    }

    # Combine results back into imputation matrix
    for (i in 1:n_imputations) {
        imputation_mat[, i] <- results[[i]]
    }

    # Stop the cluster
    parallel::stopCluster(cl)

    return(imputation_mat)
}

#' plot_single_subject
#'
#' Plot an infection time course for a single subject. The plot
#' shows when each allele was present (for each locus if applicable)
#' with annotations
#' for the probability the allele is new if `probability_new` is a column. Each
#' allele is its own row and color, and points are connected if the same allele
#' occurs repeatedly.
#' If `probability_present` is a column, imputed alleles are shown with an
#' opacity alpha equal to how often they were imputed.
#' If `probability_new` is a column, the pan-locus new
#' complexity of infection at each time is shown in a strip above the plot.
#' If an `estimated_new_infections` matrix or data frame is provided,
#' the plot title will include the estimated number of new infections.
#' Red dashes indicate a time point only known to be infection positive (i.e.,
#' no sequencing).
#' Blue dashes indicate treatment at that day if the `treatment` data frame
#' is provided.
#' If `probability_new` is a column, black dots indicate infections that are
#' new with >50% probability, open dots indicate infections that are new with
#' <50% probability.
#' If the dots are not annotated, the probabilities are >80% or
#' <20% respectively, and all other dots are annotated with their probabilities
#' of being new.
#'
#' @export
#' @param subject Subject to plot
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present`. To plot imputations,
#' `probability_present` should be a column. To plot probabilities alleles
#' are new, `probability_new` should be a column. To display allele
#' prevalence, `prevalence` should be a column.
#' @param treatments A data.frame of treatments with columns `subject` and
#' `time`.
#' @param estimated_new_infections The result of `estimate_new_infections`
#' run on the dataset.
#' @param no_imputation Disregard qPCR-only times by setting times with
#' `present == 2` to `present = 0`.
#' @param output File to which the plot should be saved
#' @param height Height of the output plot
#' @param width Width of the output plot
#' @param highlight_index Index to highlight
#' @return A plot object from ggplot2 showing the infection course for the
#' subject
#'
#' @examples
#'
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#'      subject = rep('A', 5),
#'      time = c(1, 1, 5, 15, 44))
#' treatments <- data.frame(subject = c('A', 'B'),
#'                          time = c(1, 29))
#'
#' dataset <- fill_in_dataset(dataset_in)
#'
#' dataset$probability_new <-
#'     determine_probabilities_simple(dataset)$probability_new
#'
#' estimated_new_infections <- estimate_new_infections(dataset)
#'
#' plot_single_subject('A', dataset, treatments, estimated_new_infections)
#'
#' @import dplyr
#' @import ggplot2
#'
plot_single_subject <- function(subject,
                                dataset,
                                treatments = NULL,
                                estimated_new_infections = NULL,
                                no_imputation = FALSE,
                                output = NULL,
                                height = 6,
                                width = 8,
                                highlight_index = NULL) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    check_alleles_unique_across_loci(dataset)

    if (!'probability_present' %in% colnames(dataset)) {
        dataset$probability_present <- 1 * (dataset$present == 1)
    }
    if (no_imputation) {
        dataset$probability_present <-
            ifelse(dataset$present == 2, 0, dataset$probability_present)
    }
    if (!is.null(highlight_index)) {
        dataset$highlight_index <- FALSE
        dataset$highlight_index[highlight_index] <- TRUE
    }

    subject_current <- subject
    tmp_df <- dataset %>%
        filter(subject == subject_current) %>%
        mutate(allele = factor(as.character(.data$allele),
                               levels = sort(unique(as.character(
                                   .data$allele[
                                       .data$probability_present > 0]))))) %>%
        filter(.data$allele %in% .data$allele[.data$probability_present > 0])

    time_points <- unique(tmp_df$time)

    if (nrow(tmp_df) > 0) {
        if ('locus' %in% colnames(tmp_df)) {
            tmp_df <- tmp_df %>%
                mutate(locus = as.factor(.data$locus))
            tmp_df$grayscale <- ifelse(as.numeric(tmp_df$locus) %% 2 == 0,
                                       "gray90", "gray95")
        } else {
            tmp_df$grayscale <- "gray95"
        }

        if ("probability_new" %in% colnames(tmp_df)) {
            new_infections <- ifelse(
                is.null(estimated_new_infections) |
                    all(rownames(estimated_new_infections) != subject_current),
                NA,
                rowMeans(estimated_new_infections)[
                rownames(estimated_new_infections) == subject_current])
        }

        # Filled if new
        # Open if old
        # Cross if NA
        if ("probability_new" %in% colnames(tmp_df)) {
            tmp_df <- tmp_df %>%
                mutate(
                    shape = ifelse(is.na(.data$probability_new), 13,
                       ifelse(.data$probability_new > 0.5, 16, 21))) %>%
                mutate(
                    probability_label = ifelse(.data$probability_new < 0.8 &
                       .data$probability_new > 0.2 & .data$present == 1 &
                           !is.na(.data$probability_new),
                       paste0(round(100 * .data$probability_new),
                              "%   "), "")) %>%
                group_by(.data$allele) %>%
                arrange(.data$time) %>%
                ungroup()
        } else {
            tmp_df <- tmp_df %>%
                mutate(shape = 4) %>%
                mutate(probability_label = "") %>%
                group_by(.data$allele) %>%
                arrange(.data$time) %>%
                ungroup()
        }

        if ('prevalence' %in% colnames(tmp_df)) {
            tmp_df <- tmp_df %>%
                group_by(.data$allele) %>%
                mutate(prevalence =
                        ifelse(.data$probability_present > 0.5 &
                            .data$time == max(
                            .data$time[.data$probability_present > 0.5],
                            -Inf),
                            .data$prevalence, "")) %>%
                ungroup()
        }

        tmp_df_sub <- tmp_df %>%
            filter(.data$probability_present > 0) %>%
            arrange(.data$time) %>%
            group_by(.data$allele) %>%
            mutate(min_adjacent_probability = pmin(
                ifelse(.data$probability_present == 0, 1,
                       .data$probability_present),
                ifelse(lead(.data$probability_present, default = 1) == 0, 1,
                       lead(.data$probability_present, default = 1))
            )) %>%
            ungroup()

        if ("probability_new" %in% colnames(tmp_df)) {
            tmp_df_molFOI <- tmp_df %>%
                dplyr::group_by(.data$time) %>%
                dplyr::summarise(total_new =
                    sum((.data$probability_present * .data$probability_new)[
                        .data$probability_present > 0]),
                    .groups = "drop") %>%
                dplyr::filter(!is.na(.data$total_new))

            p1 <- ggplot(tmp_df_molFOI,
                         aes(x = .data$time,
                             y = .data$total_new)) +
                geom_vline(xintercept = unique(time_points),
                           color = 'darkgray',
                           linetype = 'dashed') +
                geom_vline(xintercept = tmp_df[tmp_df$present == 2,]$time,
                           color = 'red',
                           linetype = 'dashed') +
                geom_vline(xintercept = treatments[treatments$subject ==
                                                       subject_current,]$time,
                           color = 'blue',
                           linetype = '33') +
                geom_point(aes(x = .data$time,
                               y = .data$total_new,
                               fill = 'black'),
                           stroke = 1, size = 2) +
                geom_line(aes(x = .data$time,
                              y = .data$total_new),
                          linewidth = 0.5) +
                scale_x_continuous(breaks = time_points,
                                   limits = c(min(time_points),
                                              max(time_points) * 1.1)) +
                theme_bw() +
                labs(x = "Day",
                     y = "Pan-locus\nNew Alleles",
                     title = paste0("Subject: ",
                        subject,
                        ifelse(is.na(new_infections),
                               "",
                               paste0("\nEstimated new infection events: ",
                                      round(new_infections, 1)))
                        )) +
                theme(legend.position = 'none',
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      strip.background = element_blank(),
                      strip.placement = "outside",
                      strip.text.y.left = element_text(angle = 0),
                      panel.spacing = unit(0, "lines"),
                      panel.grid = element_blank()) +
                theme(panel.background =
                          element_rect(fill = "gray90", color = NA),
                      panel.grid.major = element_blank())

            if (max(tmp_df_molFOI$total_new, na.rm=T) < 1) {
                p1 <- p1 + scale_y_continuous(limits = c(0,1), breaks = c(0,1))
            } else {
                p1 <- p1 + scale_y_continuous(limits = c(0, NA))
            }
        }

        p2 <- ggplot(tmp_df_sub, aes(x = .data$time, y = .data$allele)) +
            geom_rect(data = tmp_df,
                      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                          fill = .data$grayscale)) +
            scale_fill_manual(
                values = c("gray95" = "gray95", "gray90" = "gray90")) +
            geom_vline(xintercept = unique(time_points),
                       color = 'darkgray',
                       linetype = 'dashed') +
            geom_vline(xintercept = tmp_df[tmp_df$present == 2,]$time,
                       color = 'red',
                       linetype = 'dashed') +
            geom_vline(xintercept = treatments[
                treatments$subject == subject_current,]$time,
                color = 'blue',
                linetype = '33')

        if (!is.null(highlight_index)) {
            p2 <- p2  + geom_label(
                data = tmp_df %>% dplyr::filter(.data$highlight_index),
                aes(label = '  '),
                label.padding = unit(0.2, "lines"),
                label.size = 1,
                fill = "white",
                color = "black"
            )
        }

        p2 <- p2 +
            ggnewscale::new_scale_fill() +
            geom_line(aes(color = .data$allele,
                          group = .data$allele,
                          alpha = .data$min_adjacent_probability),
                      linewidth = 1) +
            geom_point(aes(x = .data$time,
                           y = .data$allele,
                           fill = .data$allele,
                           alpha = .data$probability_present,
                           shape = .data$shape),
                       stroke = 1, size = 2) +
            scale_x_continuous(breaks = time_points,
                               limits = c(min(time_points),
                                          max(time_points) * 1.1)) +
            scale_alpha_identity() +
            scale_shape_identity()

        if ('prevalence' %in% colnames(tmp_df)) {
            p2 <- p2 +
                ggrepel::geom_text_repel(
                    box.padding = 0,
                    data = tmp_df[tmp_df$prevalence != '',],
                    aes(x = max(time_points) * 1.02, label = .data$prevalence),
                    direction = "x",
                    hjust = 1,
                    nudge_x = 1,
                    segment.size = 0,
                    size = 3)
        }

        if ("probability_new" %in% colnames(tmp_df)) {
            p2 <- p2 +
                ggrepel::geom_text_repel(
                    box.padding = 0,
                    aes(x = .data$time,
                        y = .data$allele,
                        label = .data$probability_label),
                    vjust = 2,
                    direction = 'y',
                    size = 3)
        }
        p2 <- p2 +
            theme_bw() +
            theme(legend.position = 'none',
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text.y.left = element_text(angle = 0),
                  panel.spacing = unit(0, "lines"),
                  panel.grid = element_blank())

        if ('locus' %in% colnames(tmp_df)) {
            p2 <- p2 +
                facet_wrap(~ locus, ncol = 1, scales = 'free_y',
                           strip.position = "left") +
                theme(panel.background =
                          element_rect(fill = "gray90", color = NA),
                      panel.grid.major = element_blank()) +
                labs(x = "Day", y = 'Locus')
        } else {
            p2 <- p2 +
                theme(panel.background =
                          element_rect(fill = "gray90", color = NA),
                      panel.grid.major = element_blank()) +
                labs(x = "Day", y = 'Allele')
        }

        if ("probability_new" %in% colnames(tmp_df)) {
            if ('prevalence' %in% colnames(tmp_df)) {
                plot_out <- patchwork::wrap_plots(
                    list(
                        p1,
                        p2,
                        ggplot(data.frame(l = "Prevalence", x = 1, y = 1)) +
                            geom_text(aes(.data$x, .data$y, label = .data$l),
                                      angle = 270) +
                            theme_void() +
                            coord_cartesian(clip = "off")
                    ),
                    design = c(patchwork::area(1, 1),
                               patchwork::area(2, 1),
                               patchwork::area(2, 2, 2, 2)),
                    heights = c(1, 5), widths = c(19, 1),
                    axes = 'collect_x'
                )
            } else {
                plot_out <- patchwork::wrap_plots(
                    list(p1,
                         p2),
                    nrow = 2,
                    heights = c(1, 5),
                    axes = 'collect_x'
                )
            }
        } else {
            if ('prevalence' %in% colnames(tmp_df)) {
                plot_out <- patchwork::wrap_plots(
                    list(
                        p2,
                        ggplot(data.frame(l = "Prevalence", x = 1, y = 1)) +
                            geom_text(aes(x, y, label = l), angle = 270) +
                            theme_void() +
                            coord_cartesian(clip = "off")
                        ),
                    design = c(patchwork::area(1, 1),
                               patchwork::area(1, 2, 1, 2)),
                    widths = c(19, 1),
                    axes = 'collect_x'
                )
            } else {
                plot_out <- p2
            }
        }

        if (!is.null(output)) {
            ggplot2::ggsave(
                output, plot = plot_out, height = height, width = width)
        }
        return(plot_out)
    } else {
        return(NULL)
    }
}

#' plot_dataset
#'
#' Compute allele prevalences and then call `plot_single_subject` for every
#' subject in a dataset and save the results.
#'
#' @export
#' @param dataset A complete longitudinal dataset with columns `allele`,
#' `subject`, `time`, and `present`. To plot imputations,
#' `probability_present` should be a column. To plot probabilities alleles
#' are new, `probability_new` should be a column.
#' @param treatments A data.frame of treatments with columns `subject` and
#' `time`.
#' @param estimated_new_infections The result of `estimate_new_infections`
#' run on the dataset.
#' @param no_imputation Disregard qPCR-only times by setting times with
#' `present == 2` to `present = 0`.
#' @param output Folder to which the plots should be saved
#' @param height Height of the output plot
#' @param width Width of the output plot
#' @return A list of plot object from ggplot2 showing each per-subject
#' infection course. These are also saved by subject name in the `output`
#' folder.
#'
#' @examples
#'
#' library(dplyr)
#' library(dinemites)
#' dataset_in <- data.frame(allele = c('A', 'B', NA, NA, NA),
#'      subject = rep('A', 5),
#'      time = c(1, 1, 5, 15, 44))
#' treatments <- data.frame(subject = c('A', 'B'),
#'                          time = c(1, 29))
#'
#' dataset <- fill_in_dataset(dataset_in)
#'
#' dataset$probability_new <-
#'     determine_probabilities_simple(dataset)$probability_new
#'
#' plot_dataset(dataset, treatments)
#'
#' @import dplyr
#' @import ggplot2
#'
plot_dataset <- function(dataset,
                         treatments = NULL,
                         estimated_new_infections = NULL,
                         no_imputation = FALSE,
                         output = NULL,
                         height = 6,
                         width = 8) {
    if (any(!c("allele", "subject", "time", "present") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present")
    }

    check_alleles_unique_across_loci(dataset)

    freq_df <- dataset %>%
        dplyr::group_by(.data$subject, .data$time) %>%
        dplyr::filter(any(.data$present == 1)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$allele) %>%
        dplyr::summarise(mean(.data$present == 1),
                         .groups = "drop")

    colnames(freq_df) <- c("allele", "prevalence")

    freq_table <- freq_df$prevalence
    names(freq_table) <- freq_df$allele

    dataset <- dataset %>%
        mutate(prevalence =
            paste0(round(freq_table[as.character(.data$allele)] * 100, 1),
                   "%")) %>%
        mutate(prevalence =
            ifelse(.data$prevalence == 'NA%', "", .data$prevalence)) %>%
        ungroup()

    if (!is.null(output)) {
        dir.create(output, showWarnings = F)
        message("Saving to: ", output)
    }
    plot_list <- list()
    for (subject in unique(dataset$subject)) {
        if (!is.null(output)) {
            output_file <-
                file.path(output, paste0('subject_', subject, '.png'))
        } else {
            output_file <- NULL
        }

        plot_list[[subject]] <-
            plot_single_subject(
                subject,
                dataset,
                treatments = treatments,
                estimated_new_infections = estimated_new_infections,
                no_imputation = no_imputation,
                output = output_file,
                height = height,
                width = width)
    }
    return(plot_list)
}


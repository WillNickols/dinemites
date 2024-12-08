#' plot_single_subject
#'
#' @export
#' @param subject subject
#' @param dataset dataset
#' @param output output
#' @param treatments treatments
#' @param height height
#' @param width width
#' @return return
#' @import dplyr
#' @import ggplot2
#'
plot_single_subject <- function(subject, dataset, output = NULL, treatments = NULL, height = 6, width = 8) {
    if (any(!c("allele", "subject", "time", "present", "probability_new") %in%
            colnames(dataset))) {
        stop("dataset must contain the columns:
             allele, subject, time, present, probability_new")
    }

    if (!'probability_present' %in% colnames(dataset)) {
        dataset$probability_present <- 1 * (dataset$present == 1)
    }

    subject_current <- subject
    tmp_df <- dataset %>%
        filter(subject == subject_current) %>%
        mutate(allele = factor(as.character(allele), levels = sort(unique(as.character(allele[probability_present > 0]))))) %>%
        filter(allele %in% allele[probability_present > 0])

    time_points <- unique(tmp_df$time)

    if (nrow(tmp_df) > 0) {
        if ('locus' %in% colnames(tmp_df)) {
            tmp_df <- tmp_df %>%
                mutate(locus = as.factor(locus))
            tmp_df$grayscale <- ifelse(as.numeric(tmp_df$locus) %% 2 == 0,
                                       "gray90", "gray95")
        } else {
            tmp_df$grayscale <- "gray95"
        }

        # Filled if new
        # Open if old
        # Cross if NA
        tmp_df <- tmp_df %>%
            mutate(shape = ifelse(is.na(probability_new), 13, ifelse(probability_new > 0.5, 16, 21))) %>%
            mutate(probability_label = ifelse(probability_new < 0.8 &
                                                  probability_new > 0.2 & present == 1,
                                              paste0(round(100 * probability_new), "%   "), "")) %>%
            group_by(allele) %>%
            arrange(time) %>%
            ungroup()

        if ('prevalence' %in% colnames(tmp_df)) {
            tmp_df <- tmp_df %>%
                group_by(allele) %>%
                mutate(prevalence = ifelse(probability_present > 0.5 &
                                               time == max(time[probability_present > 0.5]),
                                           prevalence, "")) %>%
                ungroup()
        }

        tmp_df_sub <- tmp_df %>%
            filter(probability_present > 0) %>%
            arrange(time) %>%
            group_by(allele) %>%
            mutate(min_adjacent_probability = pmin(
                ifelse(probability_present == 0, 1, probability_present),
                ifelse(lead(probability_present, default = 1) == 0, 1, lead(probability_present, default = 1))
            )) %>%
            ungroup()

        tmp_df_new_molFOI <- tmp_df %>%
            dplyr::group_by(.data$time) %>%
            dplyr::summarise(total_new = sum((probability_present * probability_new)[probability_present > 0])) %>%
            dplyr::ungroup()

        p1 <- ggplot(tmp_df_new_molFOI, aes(x = time, y = total_new)) +
            geom_vline(xintercept = unique(time_points),
                       color = 'darkgray',
                       linetype = 'dashed') +
            geom_vline(xintercept = tmp_df[tmp_df$present == 2,]$time,
                       color = 'red',
                       linetype = 'dashed') +
            geom_vline(xintercept = treatments[treatments$subject == subject_current,]$time,
                       color = 'blue',
                       linetype = 'dotted') +
            geom_point(aes(x = time, y = total_new, fill = 'black'),
                       stroke = 1, size = 2) +
            geom_line(aes(x = time, y = total_new), linewidth = 0.5) +
            scale_x_continuous(breaks = time_points, limits = c(min(time_points), max(time_points) * 1.1)) +
            theme_bw() +
            labs(x = "Day", y = "New molFOI") +
            theme(legend.position = 'none',
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),           # Remove facet strip background
                  strip.placement = "outside",                  # Place facet labels outside
                  strip.text.y.left = element_text(angle = 0),  # Rotate facet labels to 0 degrees
                  panel.spacing = unit(0, "lines"),             # Remove space between facets
                  panel.grid = element_blank()) +
            theme(panel.background = element_rect(fill = "gray90", color = NA), # Light gray background
                  panel.grid.major = element_blank())

        p2 <- ggplot(tmp_df_sub, aes(x = time, y = allele)) +
            geom_rect(data = tmp_df, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                                         fill = grayscale)) +
            scale_fill_manual(values = c("gray95" = "gray95", "gray90" = "gray90")) +
            geom_vline(xintercept = unique(time_points),
                       color = 'darkgray',
                       linetype = 'dashed') +
            geom_vline(xintercept = tmp_df[tmp_df$present == 2,]$time,
                       color = 'red',
                       linetype = 'dashed') +
            geom_vline(xintercept = treatments[treatments$subject == subject_current,]$time,
                       color = 'blue',
                       linetype = 'dotted') +
            ggnewscale::new_scale_fill() +
            geom_line(aes(color = allele, alpha = min_adjacent_probability), linewidth = 1) +
            geom_point(aes(x = time, y = allele, fill = allele,
                           alpha = probability_present, shape = shape),
                       stroke = 1, size = 2) +
            scale_x_continuous(breaks = time_points, limits = c(min(time_points), max(time_points) * 1.1)) +
            scale_alpha_identity() +
            scale_shape_identity()

        if ('prevalence' %in% colnames(tmp_df)) {
            p2 <- p2 +
                ggrepel::geom_text_repel(data = tmp_df[tmp_df$prevalence != '',],
                                         aes(x = max(time_points), label = prevalence),
                                         direction = "x", hjust = 1, nudge_x = 1, segment.size = 0, size = 3)
        }
        p2 <- p2 +
            ggrepel::geom_text_repel(aes(x = time, y = allele, label = probability_label),
                                     vjust = 2, direction = 'y', size = 3) +
            theme_bw() +
            theme(legend.position = 'none',
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),           # Remove facet strip background
                  strip.placement = "outside",                  # Place facet labels outside
                  strip.text.y.left = element_text(angle = 0),  # Rotate facet labels to 0 degrees
                  panel.spacing = unit(0, "lines"),             # Remove space between facets
                  panel.grid = element_blank())

        if ('locus' %in% colnames(tmp_df)) {
            p2 <- p2 +
                facet_wrap(~ locus, ncol = 1, scales = 'free_y', strip.position = "left") +
                theme(panel.background = element_rect(fill = "gray90", color = NA), # Light gray background
                      panel.grid.major = element_blank()) +
                labs(x = "Day", y = 'Locus')
        } else {
            p2 <- p2 +
                theme(panel.background = element_rect(fill = "gray90", color = NA), # Light gray background
                      panel.grid.major = element_blank()) +
                labs(x = "Day", y = 'Allele')
        }

        plot_out <- patchwork::wrap_plots(
            p1,
            p2,
            nrow = 2,
            heights = c(1, 5),
            axes = 'collect_x'
        )

        if (!is.null(output)) {
            ggplot2::ggsave(output, plot = plot_out, height = height, width = width)
        }
        return(plot_out)
    } else {
        return(NULL)
    }
}

#' plot_dataset
#'
#' @export
#' @param subject subject
#' @param output output
#' @param treatments treatments
#' @param height height
#' @param width width
#' @return return
#' @import dplyr
#' @import data.table
#'
plot_dataset <- function(dataset, output, treatments = NULL, height = 6, width = 8) {
    freq_df <- dataset %>%
        dplyr::group_by(subject, time) %>%
        dplyr::filter(any(present == 1)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(allele) %>%
        dplyr::summarise(mean(present == 1))

    colnames(freq_df) <- c("allele", "prevalence")

    freq_table <- freq_df$prevalence
    names(freq_table) <- freq_df$allele

    dataset <- dataset %>%
        mutate(prevalence = paste0(round(freq_table[as.character(allele)] * 100, 1), "%")) %>%
        mutate(prevalence = ifelse(prevalence == 'NA%', "", prevalence)) %>%
        ungroup()

    dir.create(output, showWarnings = F)
    message("Saving to: ", output)
    for (subject in unique(dataset$subject)) {
        output_file <- file.path(output, paste0('subject_', subject, '.png'))
        plot_single_subject(subject, dataset, output = output_file, treatments = treatments, height = height, width = width)
    }
}


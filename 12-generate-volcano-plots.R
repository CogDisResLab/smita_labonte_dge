# Generate volcano plots for the datasets
#

library(tidyverse)

make_volcano <-
  function(filename,
           lfc_threshold = c(1, -1),
           pval_threshold = 0.05) {
    filepath <- file.path("Galaxy Output", "AI", filename)

    brain_region <-
      str_remove(filename, "_(Female|Male|Overall|Medication)_DGE.csv")

    comparison <- str_extract(filename, "(Female|Male|Overall|Medication)")

    figure_file <- file.path(
      "figures",
      "HISAT2",
      str_glue("{brain_region}-{str_to_title(comparison)}-Volcano.png")
    )

    if (length(lfc_threshold) == 1) {
      up_threshold  <-  lfc_threshold
      down_threshold <-  -lfc_threshold
    } else {
      up_threshold <- lfc_threshold[1]
      down_threshold <- lfc_threshold[2]
    }

    dataset <- read_csv(filepath, show_col_types = FALSE) |>
      select(SYMBOL, logFC, PValue) |>
      mutate(
        log_Pvalue = -log10(PValue),
        Significant = case_when(
          logFC >= up_threshold & PValue < pval_threshold ~ "SigUp",
          logFC <= up_threshold &
            PValue < pval_threshold ~ "SigDown",
          logFC >= up_threshold ~ "Up",
          logFC <= down_threshold ~ "Down",
          TRUE ~ "NonSig"
        )
      )

    set_x_limits <- function(vec) {
      min_val <- floor(min(vec))
      max_val <- ceiling(max(vec))

      val <- max(abs(min_val), abs(max_val))

      c(-val, val)
    }

    lfc_limits <-
      set_x_limits(dataset$logFC)
    pval_limit <- c(0, ceiling(dataset$log_Pvalue))

    g <-
      ggplot(dataset, aes(x = logFC, y = log_Pvalue, color = Significant))

    p <- g + geom_point(size = 3) +
      theme_minimal() +
      geom_vline(xintercept = lfc_threshold[1], color = "grey60") +
      geom_vline(xintercept = lfc_threshold[2], color = "grey60") +
      geom_hline(yintercept = -log10(pval_threshold), color = "grey60") +
      geom_vline(xintercept = 0, color = "black") +
      scale_x_continuous(breaks = seq(lfc_limits[1], lfc_limits[2], 1),
                         limits = lfc_limits) +
      scale_y_continuous(breaks = seq(pval_limit[1], pval_limit[2], 1),
                         limits = pval_limit) +
      scale_color_manual(
        breaks = c("SigUp", "SigDown", "Up", "Down", "NonSig"),
        values = c(
          "red",
          "blue",
          scales::muted("red", l = 20, c = 60),
          scales::muted("blue", l = 20, c = 60),
          "grey50"
        )
      ) +
      guides(color = "none") +
      xlab(expression(Log[2]~Fold~Change)) +
      ylab(expression(-Log[10]~italic(p-value))) +
      ggtitle(str_glue("{str_replace(brain_region, '_', ' ')} Volcano Plot"), subtitle = str_glue("{str_to_title(comparison)} Comparison")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, face = "bold"))

    ggsave(figure_file, plot = p, width = 10, height = 8, units = "in", dpi = 300, bg = "white")

    p

  }


files <- list.files(file.path("results", "HISAT2"), "dge")

plots <- files |>
  map(make_volcano)

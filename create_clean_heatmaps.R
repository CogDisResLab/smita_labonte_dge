# Generate heatmaps for each comparison for all brain regions.

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggalign)
  library(RColorBrewer)
})

# Write the function to generate the heatmap
human_result_path <- "results/4. Top 5% Heat Maps Post SVA/"

human_files <- list.files(human_result_path, "csv", full.names = TRUE, recursive = TRUE) |>
  set_names(~
    case_when(
      is.na(str_extract(.x, "L10005%_(\\w+)_Matrix.csv", 1)) ~ "Human_Overall",
      TRUE ~ paste0("Human_", str_extract(.x, "L10005%_(\\w+)_Matrix.csv", 1))
    ))

mouse_result_path <- "results/4. Mouse Top 5% Heat Maps No SVA"

mouse_files <- list.files(mouse_result_path, "csv", full.names = TRUE, recursive = TRUE) |>
  set_names(~
    case_when(
      is.na(str_extract(.x, "L10005%_(\\w+)_Matrix.csv", 1)) ~ "Mouse_Overall",
      TRUE ~ paste0("Mouse_", str_extract(.x, "L10005%_(\\w+)_Matrix.csv", 1))
    ))

all_files <- c(human_files, mouse_files)

clean_data <- function(filepath) {
  read_csv(filepath) |>
    select(-ends_with("Pval")) |>
    column_to_rownames("Name_GeneSymbol") |>
    as.matrix()
}

generate_heatmap <- function(filepath) {
  dataset <- clean_data(filepath)

  plot <- pheatmap(dataset,
    scale = "row", legend = TRUE, show_rownames = FALSE
  )

  plot
}

save_heatmap <- function(plot, id) {
  png_path <- file.path("figures", str_glue("{id}_heatmap_clean.png"))

  svg_path <- file.path("figures", str_glue("{id}_heatmap_clean.svg"))

  ggsave(png_path, plot)
  ggsave(svg_path, plot)
}

plots <- all_files |>
  map(generate_heatmap) |>
  map(~ pluck(.x, "gtable")) |>
  imap(~ save_heatmap(.x, .y))

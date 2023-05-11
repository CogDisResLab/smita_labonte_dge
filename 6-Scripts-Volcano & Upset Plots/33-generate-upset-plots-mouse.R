# Overlap of Differentially Expressed Genes
#

library(tidyverse)
library(UpSetR)

files <- list.files(file.path("results/2. Mouse DGEs (Galaxy) No SVA"), "csv", recursive = TRUE)

filenames <- str_remove(files, "_DGE.csv") %>%
  str_remove("\\w+/")

filepaths <- file.path("results/2. Mouse DGEs (Galaxy) No SVA", files)

results <- filepaths |>
  set_names(filenames) |>
  map( ~ read_csv(.x, show_col_types = FALSE)) |>
  imap( ~ mutate(.x, Dataset = .y)) |>
  map( ~ filter(.x, logFC >= 1, PValue < 0.05)) |>
  map( ~ pull(.x, GeneID))

Mouse_Nucleus_Accumbens <-
  list(
    Mouse_Nucleus_Accumbens_female = results$`Mouse Mouse_NA_Females`,
    Mouse_Nucleus_Accumbens_male = results$`Mouse Mouse_NA_Males`,
    Mouse_Nucleus_Accumbens_overall = results$`Mouse Mouse_NA_Overall`
  )

Mouse_Prefrontal_Cortex <-
  list(
    Mouse_Prefrontal_Cortex_female = results$`Mouse Mouse_PFC_Females`,
    Mouse_Prefrontal_Cortex_male = results$`Mouse Mouse_PFC_Males`,
    Mouse_Prefrontal_Cortex_overall = results$`Mouse Mouse_PFC_Overall`
  )

upset_MN <- fromList(Mouse_Nucleus_Accumbens)
upset_MP <- fromList(Mouse_Prefrontal_Cortex)

upset_plot_MN <- upset(upset_MN, nsets = 3)
upset_plot_MP <- upset(upset_MP, nsets = 3)

png(filename = "results/8. Mouse Upset Plots/Mouse_Nucleus_Accumbens_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_MN
dev.off()

png(filename = "results/8. Mouse Upset Plots/Mouse_Prefrontal_Cortex_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_MP
dev.off()

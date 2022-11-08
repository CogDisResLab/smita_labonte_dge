# Overlap of Differentially Expressed Genes
#

library(tidyverse)
library(UpSetR)

files <- list.files(file.path("results", "kallisto"), "dge")

filenames <- str_remove(files, "_dge.csv")

filepaths <- file.path("results", "kallisto", files)

results <- filepaths |>
  set_names(filenames) |>
  map( ~ read_csv(.x, show_col_types = FALSE)) |>
  imap( ~ mutate(.x, Dataset = .y)) |>
  map( ~ filter(.x, logFC >= 1, PValue < 0.05)) |>
  map( ~ pull(.x, SYMBOL))


Anterior_Insula <-
  list(
    Anterior_Insula_female = results$Anterior_Insula_female,
    Anterior_Insula_male = results$Anterior_Insula_male,
    Anterior_Insula_overall = results$Anterior_Insula_overall
  )

Cingulate_Gyrus <-
  list(
    Cingulate_Gyrus_female = results$Cingulate_Gyrus_female,
    Cingulate_Gyrus_male = results$Cingulate_Gyrus_male,
    Cingulate_Gyrus_overall = results$Cingulate_Gyrus_overall
  )

DLPFC <-
  list(DLPFC_female = results$DLPFC_female,
       DLPFC_male = results$DLPFC_male,
       DLPFC_overall = results$DLPFC_overall)

Nucleus_Accumbens <-
  list(
    Nucleus_Accumbens_female = results$Nucleus_Accumbens_female,
    Nucleus_Accumbens_male = results$Nucleus_Accumbens_male,
    Nucleus_Accumbens_overall = results$Nucleus_Accumbens_overall
  )

Orbitofrontal_Cortex <-
  list(
    Orbitofrontal_Cortex_female = results$Orbitofrontal_Cortex_female,
    Orbitofrontal_Cortex_male = results$Orbitofrontal_Cortex_male,
    Orbitofrontal_Cortex_overall = results$Orbitofrontal_Cortex_overall
  )

upset_A <- fromList(Anterior_Insula)
upset_C <- fromList(Cingulate_Gyrus)
upset_D <- fromList(DLPFC)
upset_N <- fromList(Nucleus_Accumbens)
upset_O <- fromList(Orbitofrontal_Cortex)


upset_plot_A <- upset(upset_A, nsets = 3)
upset_plot_C <- upset(upset_C, nsets = 3)
upset_plot_D <- upset(upset_D, nsets = 3)
upset_plot_N <- upset(upset_N, nsets = 3)
upset_plot_O <- upset(upset_O, nsets = 3)


png(filename = "figures/kallisto/Anterior_Insula_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_A
dev.off()

png(filename = "figures/kallisto/Cingulate_Gyrus_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_C
dev.off()

png(filename = "figures/kallisto/DLPFC_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_D
dev.off()

png(filename = "figures/kallisto/Nucleus_Accumbens_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_N
dev.off()

png(filename = "figures/kallisto/Orbitofrontal_Cortex_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_O
dev.off()

# Overlap of Differentially Expressed Genes
#

library(tidyverse)
library(UpSetR)

files <- list.files(file.path("Galaxy Output"), "csv", recursive = TRUE)

filenames <- str_remove(files, "_DGE.csv") %>%
  str_remove("\\w+/")

filepaths <- file.path("Galaxy Output", files)

results <- filepaths |>
  set_names(filenames) |>
  map( ~ read_csv(.x, show_col_types = FALSE)) |>
  imap( ~ mutate(.x, Dataset = .y)) |>
  map( ~ filter(.x, logFC >= 1, PValue < 0.05)) |>
  map( ~ pull(.x, GeneID))


Anterior_Insula <-
  list(
    Anterior_Insula_female = results$AI_Female,
    Anterior_Insula_male = results$AI_Male,
    Anterior_Insula_overall = results$AI_Overall,
    Anterior_Insula_medication = results$AI_Medication
  )

Cingulate_Gyrus <-
  list(
    Cingulate_Gyrus_female = results$CG_Female,
    Cingulate_Gyrus_male = results$CG_Male,
    Cingulate_Gyrus_overall = results$CG_Overall,
    Cingulate_Gyrus_medication = results$CG_Medication
  )

DLPFC <-
  list(DLPFC_female = results$DLPFC_Female,
       DLPFC_male = results$DLPFC_Male,
       DLPFC_overall = results$DLPFC_Overall,
       DLPFC_medication = results$DLPFC_Medication
       )

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

Nucleus_Accumbens <-
  list(
    Nucleus_Accumbens_female = results$NA_Female,
    Nucleus_Accumbens_male = results$NA_Male,
    Nucleus_Accumbens_overall = results$NA_Overall,
    Nucleus_Accumbens_medication = results$NA_Medication
  )

Orbitofrontal_Cortex <-
  list(
    Orbitofrontal_Cortex_female = results$OFC_Female,
    Orbitofrontal_Cortex_male = results$OFC_Male,
    Orbitofrontal_Cortex_overall = results$OFC_Overall,
    Orbitofrontal_Cortex_medication = results$OFC_Medication
  )

Subbiculum <-
  list(
    Subbiculum_female = results$Sub_Female,
    Subbiculum_male = results$Sub_Male,
    Subbiculum_overall = results$Sub_Overall,
    Subbiculum_medication = results$Sub_Medication
  )

upset_A <- fromList(Anterior_Insula)
upset_C <- fromList(Cingulate_Gyrus)
upset_D <- fromList(DLPFC)
upset_MN <- fromList(Mouse_Nucleus_Accumbens)
upset_MP <- fromList(Mouse_Prefrontal_Cortex)
upset_N <- fromList(Nucleus_Accumbens)
upset_O <- fromList(Orbitofrontal_Cortex)
upset_S <- fromList(Subbiculum)


upset_plot_A <- upset(upset_A, nsets = 4)
upset_plot_C <- upset(upset_C, nsets = 4)
upset_plot_D <- upset(upset_D, nsets = 4)
upset_plot_MN <- upset(upset_MN, nsets = 3)
upset_plot_MP <- upset(upset_MP, nsets = 3)
upset_plot_N <- upset(upset_N, nsets = 4)
upset_plot_O <- upset(upset_O, nsets = 4)
upset_plot_S <- upset(upset_S, nsets = 4)


png(filename = "figures/Upset_Plots/Anterior_Insula_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_A
dev.off()

png(filename = "figures/Upset_Plots/Cingulate_Gyrus_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_C
dev.off()

png(filename = "figures/Upset_Plots/DLPFC_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_D
dev.off()

png(filename = "figures/Upset_Plots/Mouse_Nucleus_Accumbens_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_MN
dev.off()

png(filename = "figures/Upset_Plots/Mouse_Prefrontal_Cortex_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_MP
dev.off()

png(filename = "figures/Upset_Plots/Nucleus_Accumbens_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_N
dev.off()

png(filename = "figures/Upset_Plots/Orbitofrontal_Cortex_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_O
dev.off()

png(filename = "figures/Upset_Plots/Subbiculum_Upset.png", width = 10, height = 8, units = "in", res = 300)
upset_plot_S
dev.off()

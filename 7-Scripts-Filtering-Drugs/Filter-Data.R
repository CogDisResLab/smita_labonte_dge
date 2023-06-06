# Filter and match drugs

library(tidyverse)
library(readxl)

data_files <- list.files("results/6. DrugFindR Output/", "xlsx", recursive = TRUE, full.names = TRUE)

read_drugfindr <- function(name) {
  if (str_detect(name, "L1000_\\D")) {
    df <- read_excel(name, sheet = "Discordant") |>
      mutate(gene_direction = "none")
  } else {
    up <- read_excel(name, sheet = "UpDiscordant") |>
      mutate(gene_direction = "up")

    dn <- read_excel(name, sheet = "DownDiscordant") |>
      mutate(gene_direction = "dn")

    df <- bind_rows(up, dn)
  }

  df
}

process_single_file <- function(name, region) {

  dataset <- read_drugfindr(name)

  drugs <- dataset |>
    filter(str_detect(compound, "CHEMBL", negate = T),
           str_detect(compound, "SCHEMBL", negate = T),
           str_detect(compound, "^\\d+", negate = T),
           str_detect(compound, "^[A-Z]\\d*\\w*\\-?\\s?\\d+", negate = T),
           str_detect(compound, "[Ii]nhibitor", negate = T),
           str_detect(compound, "^Broad", negate = T),
           str_detect(compound, "^BRD*", negate = T),
           str_detect(compound, "^UNII", negate = T),
           str_detect(compound, "omer", negate = T),
           str_detect(compound, "^Tyrphostin", negate = T),
           str_detect(compound, "N-[\\{\\(\\[]", negate = T),
           str_detect(compound, "^[\\{\\(\\[]+\\d", negate = T),
           str_detect(compound, "^GNF\\-", negate = T),
           str_detect(compound, "CVF\\-", negate = T),
           str_detect(compound, "NVP\\-", negate = T),
           str_detect(compound, "^N,N", negate = T),
           str_detect(compound, "^EMF\\-", negate = T),
           str_detect(compound, "^GSK\\s?", negate = T),
    ) |>
    mutate(region = region)
}

process_comparison <- function(files, name) {
  out <- files |>
    set_names(~ str_extract(.x, "[A-Z]{2,5}")) |>
    imap(~ process_single_file(.x, .y)) |>
    map(~ mutate(.x, comparison = name))

  out
}

quintiliz <- function(vector, value) {

}

combinations <- expand_grid(a = c("L1000", "L1000_5%", "L1000_10%"), b = c("Overall", "Males", "Females", "Males_Meds"), c = "drugfindr") |>
  mutate(filter = str_glue("{a}_{b}_{c}"),
         name_a = case_when(
           a == "L1000" ~ "all",
           str_detect(a, "5%") ~ "05p",
           str_detect(a, "10%") ~ "10p"
         ),
         name_b = str_to_lower(b),
         name = str_glue("{name_a}_{name_b}")) |>
  select(name, filter) |>
  deframe()

data_files <- combinations |>
  map(~ list.files("results/6. DrugFindR Output/", pattern = .x, recursive = TRUE, full.names = TRUE, ignore.case = TRUE))

output <- data_files |>
  imap(~ process_comparison(.x, .y)) |>
  map(~ imap_dfr(.x, ~ count(.x, cellline, time, concentration), .id = "region")) |>
  imap(~ write_csv(.x, str_glue("results/cell-line-exploration/{.y}_comparison_drugs.csv")))



x <- output[["all_overall"]] |>
  group_by(region) |>
  mutate(quintile_rank = ntile(desc(n), 5),
         triple = str_glue("{cellline}-{time}-{concentration}")) |>
  select(-cellline, -time, -concentration, -n) |>
  distinct() |>
  pivot_wider(names_from = region, values_from = quintile_rank, values_fill = 100) |>
  pivot_longer(cols =  c(AI, CG, DLPFC, NAC, OFC, SUB), names_to = "region", values_to = "quintile_rank") |>
  ungroup() |>
  group_by(triple) |>
  mutate(mean_rank = mean(quintile_rank, na.rm = TRUE)) |>
  pivot_wider(names_from = region, values_from = quintile_rank, values_fill = 100) |>
  arrange(mean_rank)


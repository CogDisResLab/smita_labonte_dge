# Filter and match drugs

library(tidyverse)
library(readxl)
library(writexl)
library(httr)
library(jsonlite)

read_drugfindr <- function(name) {
  if (str_detect(name, "L1000_\\D")) {
    df <- read_excel(name, sheet = "Discordant") |>
      mutate(gene_direction = "none")
  } else {
    up_dis <- read_excel(name, sheet = "UpDiscordant") |>
      mutate(gene_direction = "up")

    dn_dis <- read_excel(name, sheet = "DownDiscordant") |>
      mutate(gene_direction = "dn")

    up_con <- read_excel(name, sheet = "UpConcordant") |>
      mutate(gene_direction = "up")

    dn_con <- read_excel(name, sheet = "DownConcordant") |>
      mutate(gene_direction = "dn")

    df <- bind_rows(up_dis, dn_dis, up_con, dn_con)
  }

  df
}

process_single_file <- function(name, region) {
  dataset <- read_drugfindr(name)

  drugs <- dataset |>
    filter(
      str_detect(compound, fixed("CHEMBL"), negate = TRUE),
      str_detect(compound, fixed("SCHEMBL"), negate = TRUE),
      str_detect(compound, "^\\d+", negate = TRUE), # nolint: nonportable_path_linter.
      str_detect(compound, "^[A-Z]\\d*\\w*\\-?\\s?\\d+", negate = TRUE),
      str_detect(compound, "[Ii]nhibitor", negate = TRUE),
      str_starts(compound, fixed("Broad"), negate = TRUE),
      str_starts(compound, fixed("BRD"), negate = TRUE),
      str_starts(compound, fixed("UNII"), negate = TRUE),
      str_detect(compound, fixed("omer"), negate = TRUE),
      str_starts(compound, fixed("Tyrphostin"), negate = TRUE),
      str_detect(compound, "N-[\\{\\(\\[]", negate = TRUE), # nolint: nonportable_path_linter.
      str_detect(compound, "^[\\{\\(\\[]+\\d", negate = TRUE), # nolint: nonportable_path_linter.
      str_starts(compound, fixed("GNF-"), negate = TRUE),
      str_detect(compound, fixed("CVF-"), negate = TRUE),
      str_detect(compound, fixed("NVP-"), negate = TRUE),
      str_starts(compound, fixed("N,N"), negate = TRUE),
      str_starts(compound, fixed("EMF-"), negate = TRUE),
      str_detect(compound, "^GSK\\s?", negate = TRUE)
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

quintilize <- function(dataset) {
  dataset |>
    group_by(region) |>
    mutate(
      quintile_rank = ntile(desc(n), 5L),
      triple = str_glue("{cellline}-{time}-{concentration}")
    ) |>
    select(-cellline, -time, -concentration, -n) |>
    distinct() |>
    pivot_wider(names_from = region, values_from = quintile_rank, values_fill = 100L) |>
    pivot_longer(cols = c(AI, CG, DLPFC, NAC, OFC, SUB), names_to = "region", values_to = "quintile_rank") |>
    ungroup() |>
    group_by(triple) |>
    mutate(mean_rank = mean(quintile_rank, na.rm = TRUE)) |>
    pivot_wider(names_from = region, values_from = quintile_rank, values_fill = 100L) |>
    arrange(mean_rank) |>
    ungroup()
}

filterize <- function(dataset) {
  dataset |>
    mutate(mean_rank_rank = ntile(mean_rank, 5L)) |>
    filter(mean_rank_rank == 1L) |>
    select(-mean_rank_rank)
}

intersect_triples <- function(data_list) {
  data_list |>
    map(~ pull(.x, "triple")) |>
    reduce(intersect)
}

get_ilincs_metadata <- function(X) {
  url <- "http://www.ilincs.org/api/SignatureMeta/findMany"
  body <- list(signatures = toJSON(X))

  metadata <- POST(url, body = body, encode = "json") %>%
    content(as = "text") %>%
    fromJSON() %>%
    pluck("data") %>%
    select(TargetSignature = signatureid, tissue, integratedMoas, GeneTargets) # where(~ any(!is.na(.x)))
}

combinations <- expand_grid(
  a = c("L1000", "L1000_5%", "L1000_10%"),
  b = c("Overall", "Males", "Females", "Males_Meds"), c = "drugfindr"
) |>
  mutate(
    filter = str_glue("{a}_{b}_{c}"),
    name_a = case_when(
      a == "L1000" ~ "entire",
      str_detect(a, fixed("5%")) ~ "05p",
      str_detect(a, fixed("10%")) ~ "10p"
    ),
    name_b = str_to_lower(b),
    name = str_glue("{name_a}_{name_b}")
  ) |>
  select(name, filter) |>
  deframe()

data_files <- combinations |>
  map(~ list.files("results/6. DrugFindR Output/",
    pattern = .x, recursive = TRUE,
    full.names = TRUE, ignore.case = TRUE
  ))

# output <- data_files |>
#   imap(~ process_comparison(.x, .y)) |>
#   map(~ imap_dfr(.x, ~ count(.x, cellline, time, concentration), .id = "region")) |>
#   imap(~ write_csv(
#     .x,
#     str_glue("results/cell-line-exploration/triple_comparison_count/{.y}_comparison_drugs.csv")
#   )) |>
#   map(~ quintilize(.x)) |>
#   imap(~ write_csv(
#     .x,
#     str_glue("results/cell-line-exploration/triple_comparison_rank/{.y}_comparison_triple_rank.csv")
#   )) |>
#   map(~ filterize(.x)) |>
#   assign("all_output_data", value = _) |>
#   imap(~ write_csv(
#     .x,
#     str_glue("results/cell-line-exploration/triple_comparison_top_quintile/{.y}_comparison_triple_top_quintile.csv")
#   ))

# all_output_data |>
#   write_xlsx("results/cell-line-exploration/triple_ranked_comparisons.xlsx")

# cell_lines_entire <- output[str_detect(names(output), "entire_")] |>
#   intersect_triples() |>
#   str_c(collapse = "\n") |>
#   write_file("data/common_cell_lines_overall.txt") |>
#   str_split("\\n", simplify = TRUE)


# cell_lines_05p <- output[str_detect(names(output), fixed("05p_"))] |>
#   intersect_triples() |>
#   str_c(collapse = "\n") |>
#   write_file("data/common_cell_lines_05p.txt") |>
#   str_split("\\n", simplify = TRUE)

# cell_lines_10p <- output[str_detect(names(output), fixed("10p_"))] |>
#   intersect_triples() |>
#   str_c(collapse = "\n") |>
#   write_file("data/common_cell_lines_10p.txt") |>
#   str_split("\\n", simplify = TRUE)


# filtered_entire <- data_files[str_detect(names(data_files), fixed("entire"))] |>
#   imap(~ process_comparison(.x, .y)) |>
#   map(~ imap_dfr(.x, ~ mutate(.x, triple = str_glue("{cellline}-{time}-{concentration}")), .id = "region")) |>
#   map_dfr(~ filter(.x, triple %in% cell_lines_entire), .id = "source") |>
#   write_csv("results/cell-line-exploration/filtered_entire_drug_list.csv")

# filtered_05p <- data_files[str_detect(names(data_files), fixed("05p"))] |>
#   imap(~ process_comparison(.x, .y)) |>
#   map(~ imap_dfr(.x, ~ mutate(.x, triple = str_glue("{cellline}-{time}-{concentration}")), .id = "region")) |>
#   map_dfr(~ filter(.x, triple %in% cell_lines_05p), .id = "source") |>
#   write_csv("results/cell-line-exploration/filtered_05p_drug_list.csv")

# filtered_10p <- data_files[str_detect(names(data_files), fixed("10p"))] |>
#   imap(~ process_comparison(.x, .y)) |>
#   map(~ imap_dfr(.x, ~ mutate(.x, triple = str_glue("{cellline}-{time}-{concentration}")), .id = "region")) |>
#   map_dfr(~ filter(.x, triple %in% cell_lines_10p), .id = "source") |>
#   write_csv("results/cell-line-exploration/filtered_10p_drug_list.csv")

# Next steps will follow what we did in the previous paper/3-pod but with adjustments for your specific study.  
#
# Using the top 5% signature list, for each region and each comparison (e.g. AI Fm CTL v MDD): Include only signatures < -0.50.

filtered_05p_sinead <- data_files[str_detect(names(data_files), fixed("05p"))] |>
  imap(~ process_comparison(.x, .y)) |>
  map(~ bind_rows(.x, .id = "region")) |>
  bind_rows(.id = "comparison") |>
  select(signatureid, compound, cellline, similarity, region, comparison) |>
  group_by(region, comparison, compound, cellline) |>
  filter(abs(similarity) == max(abs(similarity))) |>
  ungroup() |>
  write_csv("results/cell-line-exploration/filtered_05p_sinead_drug_list_full.csv")

metadata <- get_ilincs_metadata(unique(filtered_05p_sinead$signatureid)) |>
  select(signatureid = TargetSignature, mechanism = integratedMoas)

final_filtered_05p_sinead <- filtered_05p_sinead |>
  left_join(metadata, by = "signatureid") |>
  select(-signatureid) |>
  filter(abs(similarity) > 0.321) |> # This is the only change from the previous script
  nest(.by = c(region, comparison, compound)) |>
  mutate(
    num_celllines = map_int(data, ~ length(unique(.x[["cellline"]]))),
    num_concordant = map_int(data, ~ sum(.x[["similarity"]] > 0L)),
    num_discordant = map_int(data, ~ sum(.x[["similarity"]] < 0L)),
    ratio = round(num_concordant / num_discordant, 8)
  ) |>
  filter(num_celllines > 1L) |>
  mutate(pivoted = map(data, ~ pivot_wider(.x, names_from = cellline, values_from = similarity))) |>
  select(-data, -num_celllines) |>
  unnest(pivoted) |>
  write_csv("results/cell-line-exploration/filtered_05p_sinead_drug_list_by_celllines.csv") |>
  nest(.by = c(region, comparison)) |>
  mutate(
    filename = str_glue("results/cell-line-exploration/filtered_05p_sinead_drug_list_{region}_{comparison}.csv")
  ) |>
  mutate(data = map2(data, filename, ~ write_csv(.x, .y)))

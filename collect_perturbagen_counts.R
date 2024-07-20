# Collect the perturbagen counts for each brain region and comparison and
# summarise them in a single table.

library(tidyverse)
library(readxl)

count_rows <- function(file) {
  sheets <- excel_sheets(file) |>
    set_names()
  
  counted <- sheets |>
    map(~ read_excel(file, sheet = .x) |>
          nrow()) |>
    enframe(name = "sheet", value = "count") |>
    unnest(cols = c(count)) |>
    mutate(file = basename(file))
}

files <- list.files("results", "drugfindr.xlsx", recursive = TRUE, full.names = TRUE) |> set_names(basename)

counted_files <- files |>
  map(count_rows) |>
  bind_rows()

counted_augmented <- counted_files |>
  mutate(region = str_extract(file, "(\\w+)_L1000", 1),
        percentage = str_extract(file, "L1000_(\\d+%)?"),
        percentage = if_else(percentage == "L1000_", "100%", percentage),
        comparison = str_extract(file, "L1000_(\\d+%)?_(\\w+)_drugfindr.xlsx", 2),
        comparison = case_when(
          comparison == "males" ~ "Males",
          comparison == "males_meds" ~ "Males_Meds",
          .default = comparison
        ),
        species = if_else(str_detect(file, "Mouse"), "Mouse", "Human")) |>
        select(-file)

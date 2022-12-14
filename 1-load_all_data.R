# Load all data and create input files

library(tidyverse)

files_kallisto <-
  list.files(file.path("raw", "kallisto"), , pattern = "Count")

process_count_data <- function(filename, aligner = "kallisto") {
  brain_region <- str_remove(filename, "_CLEAN_Count_Matrix.csv")

  filepath <- file.path("raw", aligner, filename)

  outfile_name <- file.path("data",
                            aligner,
                            str_glue("{brain_region}_clean_count.csv"))

  count <- read_csv(filepath, show_col_types = FALSE) |>
    group_by(SYMBOL) |>
    summarise(across(where(is.numeric), sum)) |>
    write_csv(outfile_name)
}

files_kallisto |>
  walk( ~ process_count_data(.x))

metadata <-
  read_csv("raw/GSE102556-SRA-Metadata.txt", show_col_types = FALSE) |>
  select(
    Run,
    tissue,
    phenotype,
    gender,
    drugs,
    Cause_of_death,
    medication,
    medication_type,
    ph,
    pmi,
    rin
  ) |>
  write_csv("data/GSE102556-metadata.csv")

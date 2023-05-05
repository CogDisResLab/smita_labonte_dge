# Extract mouse metadata features we want

library(tidyverse)

metadata <-
  read_csv("raw/GSE102556-SRA-Mouse-Metadata.txt", show_col_types = FALSE) |>
  select(
    Run,
    TISSUE,
    Phenotype,
    gender,
    RIN
  ) |>
  write_csv("data/GSE102556-SRA-Mouse-Metadata.csv")

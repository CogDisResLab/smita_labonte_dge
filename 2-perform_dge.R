# Differential Expression Analysis

library(tidyverse)
library(edgeR)

perform_dge <- function(filename, aligner = "kallisto") {
  filepath <- file.path("data", aligner, filename)

  brain_region <- filename |>
    str_remove("_clean_count.csv")

  count_data <- read_csv(filepath, show_col_types = FALSE)

  metadata <-
    read_csv("data/GSE102556-metadata.csv", show_col_types = FALSE) |>
    select(Run, phenotype, gender, medication) |>
    mutate(
      medication = case_when(
        medication == "yes" ~ "on",
        medication == "no" ~ "off",
        is.na(medication) ~ "na",
        TRUE ~ NA_character_
      )
    ) |>
    filter(Run %in% names(count_data)) |>
    arrange(phenotype, gender) |>
    mutate(group = str_c(gender, phenotype, medication, sep = "_"),
           rowname = Run) |>
    column_to_rownames()

  count_subset <- count_data |>
    select(SYMBOL, all_of(metadata$Run))


  dge <-
    DGEList(
      counts = count_subset[, 2:ncol(count_subset)],
      samples = metadata,
      genes = count_subset[1],
      group = metadata$group
    )

  keep <- filterByExpr(dge)

  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

  design_frame <- model.frame(~ 0 + group, data = metadata)

  design <- model.matrix(~ 0 + group, data = design_frame)

  dge_filtered <- calcNormFactors(dge_filtered)

  dge_filtered <- estimateDisp(dge_filtered, design = design)

  fit <- glmQLFit(dge_filtered, design = design)
}

get_specified_groups <-
  function(design,
           phenotype,
           sex = NULL,
           medication = NULL) {
    group_names <- colnames(design)

    p <- str_detect(group_names, phenotype)

    if (!is.null(sex)) {
      s <- str_detect(group_names, str_c("group", sex))
    } else {
      s <- rep(TRUE, length(p))
    }

    if (!is.null(medication)) {
      m <- str_detect(group_names, medication)
    } else {
      m <- rep(TRUE, length(p))
    }

    selected <- group_names[p & s & m]

    str_c(selected, collapse = " + ")

  }

perform_contrasts <- function(fit_obj) {
  design <- fit_obj$design

  contrast <- makeContrasts(
    contrasts = c(
      med = str_glue(
        '(({get_specified_groups(design, "MDD", medication = "on")}) / 2) - ({get_specified_groups(design, "MDD", medication = "off")})'
      ),
      male = str_glue(
        '({get_specified_groups(design, "MDD", "male")}) - ({get_specified_groups(design, "CTRL", "male")})'
      ),
      female = str_glue(
        '({get_specified_groups(design, "MDD", "female")}) - ({get_specified_groups(design, "CTRL", "female")})'
      ),
      overall = str_glue(
        '({get_specified_groups(design, "MDD")}) - ({get_specified_groups(design, "CTRL")})'
      )
    ),
    levels = design
  )

  medication_dge <- glmQLFTest(fit_obj, contrast = contrast[, 1])
  male_dge <- glmQLFTest(fit_obj, contrast = contrast[, 2])
  female_dge <- glmQLFTest(fit_obj, contrast = contrast[, 3])
  overall_dge <- glmQLFTest(fit_obj, contrast = contrast[, 4])

  out <- list(
    male_dge = topTags(male_dge, n = Inf)$table,
    female_dge = topTags(female_dge, n = Inf)$table,
    overall_dge = topTags(overall_dge, n = Inf)$table,
    medication_dge = topTags(medication_dge, n = Inf)$table
  )

  out
}


files_kallisto <- list.files(file.path("data", "kallisto"))

brain_regions <- files_kallisto |>
  str_remove("_clean_count.csv")

files_kallisto <- files_kallisto |> set_names(brain_regions)

result <- files_kallisto |>
  map(~ perform_dge(.x)) |>
  map(~ perform_contrasts(.x)) |>
  unlist(recursive = FALSE) |>
  imap(~ write_csv(.x,
                   file.path(
                     "results",
                     "kallisto",
                     str_glue('{str_replace(.y, "\\\\.", "_")}.csv')
                   )))

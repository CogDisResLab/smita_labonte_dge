# Differential Expression Analysis

library(tidyverse)
library(edgeR)

perform_dge <- function(filename, aligner = "kallisto") {

  filepath <- file.path("data", aligner, filename)

  brain_region <- filename |>
    str_remove("_clean_count.csv")

  count_data <- read_csv(filepath, show_col_types = FALSE)

  metadata <- read_csv("data/GSE102556-metadata.csv", show_col_types = FALSE) |>
    select(Run, phenotype, gender, medication) |>
    mutate(medication = if_else(is.na(medication), "Unknown", str_to_title(medication))) |>
    filter(Run %in% names(count_data)) |>
    arrange(phenotype, gender) |>
    mutate(group = str_c(gender, phenotype, medication, sep = "_"),
           rowname = Run) |>
    column_to_rownames()

  count_subset <- count_data |>
    select(SYMBOL, all_of(metadata$Run))


  dge <- DGEList(count = count_subset[,2:ncol(count_subset)], samples = metadata, genes = count_subset[1], group = metadata$group)

  keep <- filterByExpr(dge)

  dge_filtered <- dge[keep, ,keep.lib.sizes=FALSE]

  design_frame <- model.frame(~ 0 + group, data = metadata)

  design <- model.matrix(~ 0 + group, data = design_frame)

  dge_filtered <- calcNormFactors(dge_filtered)

  dge_filtered <- estimateDisp(dge_filtered, design = design)

  fit <- glmQLFit(dge_filtered, design = design)
}

perform_contrasts <- function(fit_obj) {

  design <- fit_obj$design

  contrast <- makeContrasts(
    male = (groupmale_MDD_Yes + groupmale_MDD_No + groupmale_MDD_Unknown) - (groupmale_CTRL_No + groupmale_CTRL_Unknown),
    female = (groupfemale_MDD_Yes + groupfemale_MDD_No + groupfemale_MDD_Unknown) - (groupfemale_CTRL_Yes + groupfemale_CTRL_No + groupfemale_CTRL_Unknown),
    overall = (groupmale_MDD_Yes + groupmale_MDD_No + groupmale_MDD_Unknown + groupfemale_MDD_Yes + groupfemale_MDD_No + groupfemale_MDD_Unknown) - (groupmale_CTRL_No + groupmale_CTRL_Unknown + groupfemale_CTRL_Yes + groupfemale_CTRL_No + groupfemale_CTRL_Unknown),
    mdd_medication = (groupfemale_MDD_Yes + groupmale_MDD_Yes) - (groupmale_MDD_No + groupfemale_MDD_No),
    levels = design)

  male_dge <- glmQLFTest(fit_obj, contrast = contrast[, "male"])
  female_dge <- glmQLFTest(fit_obj, contrast = contrast[, "female"])
  overall_dge <- glmQLFTest(fit_obj, contrast = contrast[, "overall"])
  medication_dge <- glmQLFTest(fit_obj, contrast = contrast[, "mdd_medication"])

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
  map(~ safely(perform_contrasts(.x))) |>
  unlist(recursive = FALSE) |>
  imap(~ write_csv(
    .x,
    file.path(
      "results",
      "kallisto",
      str_glue('{str_replace(.y, "\\\\.", "_")}.csv')
    )
  ))

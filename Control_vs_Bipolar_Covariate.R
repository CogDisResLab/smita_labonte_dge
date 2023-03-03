library(DEFormats)
library(tidyverse)
library(edgeR)
library(variancePartition)
library(BiocParallel)
library(scales)
library(org.Hs.eg.db)

df <-
  read_csv("Galaxy Count Matrix Mapped/NA/NA_Overall_Mapped_Final.csv") |>
  mutate(ENTREZID = as.character(Geneid))

symbols <-
  select(
    org.Hs.eg.db,
    df$ENTREZID,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "ENTREZID"
  )

symboled_df <- df |>
  inner_join(symbols) |>
  dplyr::select(-Geneid, -ENTREZID, SYMBOL, starts_with("SRR")) |>
  dplyr::relocate(SYMBOL, .before = starts_with("SRR")) |>
  group_by(SYMBOL) |>
  summarise(across(starts_with("SRR"), ~ sum(.x))) |>
  ungroup() |>
  filter(!is.na(SYMBOL)) |>
  mutate(rowname = SYMBOL) |>
  column_to_rownames() |>
  rename_with(.cols = SYMBOL, ~ str_to_title(.x))


metadata <- read_csv("data/GSE102556-NA-metadata.csv") |>
  mutate(
    tissue = as_factor(str_extract(tissue, "\\((\\w+)\\)", group = 1)),
    gender = as_factor(gender),
    diagnosis = as_factor(phenotype),
    cod = as_factor(str_to_lower(Cause_of_death)),
    medication = if_else(is.na(medication), "na", medication),
    drug_use = if_else(is.na(drugs), "na", drugs),
    medication_type = if_else(is.na(medication_type), "na", medication_type),
    medication_type = as_factor(str_to_lower(medication_type)),
    medication = as_factor(medication),
    drug_use = as_factor(drug_use),
    rowname = Run
  ) |>
  column_to_rownames() |>
  dplyr::select(run = Run,
                tissue,
                gender,
                diagnosis,
                medication,
                medication_type,
                cod,
                drug_use,
                ph,
                pmi,
                rin)

group <- factor(metadata$diagnosis)


dge <- DGEList(
  counts = symboled_df[, 2:51],
  samples = metadata,
  genes = symboled_df[,1],
  group = group
)

keep <- filterByExpr(dge, group = group)

dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

dge_filtered <- calcNormFactors(dge_filtered)

dge_normalized <- cpm(dge_filtered, normalized.lib.sizes = FALSE)

# First find canonical correlation between the variables highly correlated will be removed

formula <- ~ gender + diagnosis + medication + medication_type + cod + drug_use + ph + pmi + rin

correlation_pairs <- canCorPairs(formula, as.data.frame(metadata))
plotCorrMatrix(correlation_pairs)

metadata_rescaled <- metadata |>
  mutate(across(where(is.numeric), ~ rescale(.x, to = c(0, 100))))


variance_formula <-
  ~ (1|gender) + (1|diagnosis) + (1|medication) + (1|cod) + (1|drug_use) + ph + pmi + rin

varPart <- fitExtractVarPartModel(dge_normalized, variance_formula, metadata)

apply(varPart, 2, \(x) round(mean(x), 4) * 100)

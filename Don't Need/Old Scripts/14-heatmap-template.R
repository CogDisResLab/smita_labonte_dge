# Generate Heatmaps
#

library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)

ofc_dge <- read_csv("Galaxy Output/OFC/OFC_Overall_dge.csv")

ofc_count <- read_csv("Galaxy Count Matrix Output/OFC/OFC_Overall_Counts.csv")

ofc_top_dge <- ofc_dge |>
  arrange(FDR) |>
  slice_head(n = 10)

ofc_bot_dge <- ofc_dge |>
  arrange(FDR) |>
  slice_tail(n = 10)

filtered_count <- count |>
  filter(SYMBOL %in% c(ofc_top_dge$SYMBOL, ofc_bot_dge$SYMBOL))

count_matrix <- filtered_count |>
  column_to_rownames("SYMBOL") |>
  as.matrix()

ctmx <- count_matrix + 1

lctmtx <- log(ctmx)

# heatmap(lctmtx)

pheatmap(lctmtx)

# ggdata <- ctmx |>
#   as.data.frame() |>
#   rownames_to_column("SYMBOL") |>
#   pivot_longer(cols = starts_with("SRR"), names_to = "Sample", values_to = "Count") |>
#   mutate(LogCount = log(Count))
#
# g <- ggplot(data = ggdata, mapping = aes(Sample, SYMBOL, fill = LogCount))
#
# g + geom_tile() +
#   scale_fill_gradient(low = "white", high = "red")
#
# kaleidoscope_metadata <- read_csv("~/Downloads/kaleidoscope_lookup_updated_meta.csv")
#
# kaleidoscope_lookup_data <- read_csv("~/Downloads/SCZ_Lookup.csv")
#
# x <- kaleidoscope_lookup_data |>
#   select(HGNC_Symbol, DataSet, Log2FC) |>
#   pivot_wider(names_from = "DataSet", values_from = "Log2FC", values_fill = 0) |>
#   column_to_rownames("HGNC_Symbol") |>
#   as.matrix()
#
# y <- read_csv("~/Downloads/combined_annotated_source_pathway.csv") |>
#   pivot_longer(cols = -Gene, names_to = "Pathway", values_to = "Database")
#
# heatmap(x)

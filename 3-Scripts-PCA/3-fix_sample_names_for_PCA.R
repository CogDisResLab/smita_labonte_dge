# Remapping the names of the samples so data numbers in count matrix are mapped to SRR IDs

library(tidyverse)

# AI Overall ------------------------------------------------------------------
# Read the Mapping File for SRR to Galaxy Data numbers
mapping_AI <- read_csv("Galaxy Count Matrix Mapped/AI/AI_Mappings.csv")

# Read the Overall Count Data
counts_AI <- read_csv("Galaxy Count Matrix Output/AI/AI_Overall_Counts.csv")

# Change the Column number to the SRR ID for each sample
clean_counts_AI <- counts_AI %>%
  pivot_longer(cols = starts_with("Column"), names_to = "Dataset") %>%
  mutate(Dataset_Name = {str_extract(Dataset, " \\d{3}:") %>% str_replace(" ", "") %>% str_remove(":") %>% as.numeric()}) %>%
  inner_join(mapping_AI, by = c(Dataset_Name = "HISAT")) %>%
  select(Geneid, SRR, value) %>%
  pivot_wider(names_from = SRR)

# Create the data matrix that will be used in PCA plot generation code
clean_counts_AI %>% write_csv("Galaxy Count Matrix Mapped/AI/AI_Overall_Mapped_Final.csv")

# To filter the data matrix to only include males
# clean_counts %>% filter(gender == "male")


# CG Overall ------------------------------------------------------------------
# Read the Mapping File for SRR to Galaxy Data numbers
mapping_CG <- read_csv("Galaxy Count Matrix Mapped/CG/CG_Mappings.csv")

# Read the Overall Count Data
counts_CG <- read_csv("Galaxy Count Matrix Output/CG/CG_Overall_Counts.csv")

# Change the Column number to the SRR ID for each sample
clean_counts_CG <- counts_CG %>%
  pivot_longer(cols = starts_with("Column"), names_to = "Dataset") %>%
  mutate(Dataset_Name = {str_extract(Dataset, " \\d{3}:") %>% str_replace(" ", "") %>% str_remove(":") %>% as.numeric()}) %>%
  inner_join(mapping_CG, by = c(Dataset_Name = "HISAT")) %>%
  select(Geneid, SRR, value) %>%
  pivot_wider(names_from = SRR)

# Create the data matrix that will be used in PCA plot generation code
clean_counts_CG %>% write_csv("Galaxy Count Matrix Mapped/CG/CG_Overall_Mapped_Final.csv")


# DLPFC Overall ------------------------------------------------------------------
# Read the Mapping File for SRR to Galaxy Data numbers
mapping_DLPFC <- read_csv("Galaxy Count Matrix Mapped/DLPFC/DLPFC_Mappings.csv")

# Read the Overall Count Data
counts_DLPFC <- read_csv("Galaxy Count Matrix Output/DLPFC/DLPFC_Overall_Counts.csv")

# Change the Column number to the SRR ID for each sample
clean_counts_DLPFC <- counts_DLPFC %>%
  pivot_longer(cols = starts_with("Column"), names_to = "Dataset") %>%
  mutate(Dataset_Name = {str_extract(Dataset, " \\d{3}:") %>% str_replace(" ", "") %>% str_remove(":") %>% as.numeric()}) %>%
  inner_join(mapping_DLPFC, by = c(Dataset_Name = "HISAT")) %>%
  select(Geneid, SRR, value) %>%
  pivot_wider(names_from = SRR)

# Create the data matrix that will be used in PCA plot generation code
clean_counts_DLPFC %>% write_csv("Galaxy Count Matrix Mapped/DLPFC/DLPFC_Overall_Mapped_Final.csv")


# NA Overall ------------------------------------------------------------------
# Read the Mapping File for SRR to Galaxy Data numbers
mapping_NA <- read_csv("Galaxy Count Matrix Mapped/NA/NA_Mappings.csv")

# Read the Overall Count Data
counts_NA <- read_csv("Galaxy Count Matrix Output/NA/NA_Overall_Counts.csv")

# Change the Column number to the SRR ID for each sample
clean_counts_NA <- counts_NA %>%
  pivot_longer(cols = starts_with("Column"), names_to = "Dataset") %>%
  mutate(Dataset_Name = {str_extract(Dataset, " \\d{3}:") %>% str_replace(" ", "") %>% str_remove(":") %>% as.numeric()}) %>%
  inner_join(mapping_NA, by = c(Dataset_Name = "HISAT")) %>%
  select(Geneid, SRR, value) %>%
  pivot_wider(names_from = SRR)

# Create the data matrix that will be used in PCA plot generation code
clean_counts_NA %>% write_csv("Galaxy Count Matrix Mapped/NA/NA_Overall_Mapped_Final.csv")


# OFC Overall ------------------------------------------------------------------
# Read the Mapping File for SRR to Galaxy Data numbers
mapping_OFC <- read_csv("Galaxy Count Matrix Mapped/OFC/OFC_Mappings.csv")

# Read the Overall Count Data
counts_OFC <- read_csv("Galaxy Count Matrix Output/OFC/OFC_Overall_Counts.csv")

# Change the Column number to the SRR ID for each sample
clean_counts_OFC <- counts_OFC %>%
  pivot_longer(cols = starts_with("Column"), names_to = "Dataset") %>%
  mutate(Dataset_Name = {str_extract(Dataset, " \\d{3}:") %>% str_replace(" ", "") %>% str_remove(":") %>% as.numeric()}) %>%
  inner_join(mapping_OFC, by = c(Dataset_Name = "HISAT")) %>%
  select(Geneid, SRR, value) %>%
  pivot_wider(names_from = SRR)

# Create the data matrix that will be used in PCA plot generation code
clean_counts_OFC %>% write_csv("Galaxy Count Matrix Mapped/OFC/OFC_Overall_Mapped_Final.csv")


# SUB Overall ------------------------------------------------------------------
# Read the Mapping File for SRR to Galaxy Data numbers
mapping_SUB <- read_csv("Galaxy Count Matrix Mapped/Sub/Sub_Mappings.csv")

# Read the Overall Count Data
counts_SUB <- read_csv("Galaxy Count Matrix Output/Sub/Sub_Overall_Counts.csv")

# Change the Column number to the SRR ID for each sample
clean_counts_SUB <- counts_SUB %>%
  pivot_longer(cols = starts_with("Column"), names_to = "Dataset") %>%
  mutate(Dataset_Name = {str_extract(Dataset, " \\d{3}:") %>% str_replace(" ", "") %>% str_remove(":") %>% as.numeric()}) %>%
  inner_join(mapping_SUB, by = c(Dataset_Name = "HISAT")) %>%
  select(Geneid, SRR, value) %>%
  pivot_wider(names_from = SRR)

# Create the data matrix that will be used in PCA plot generation code
clean_counts_SUB %>% write_csv("Galaxy Count Matrix Mapped/Sub/Sub_Overall_Mapped_Final.csv")


# Mouse NA Overall ------------------------------------------------------------------
# Read the Mapping File for SRR to Galaxy Data numbers
mapping_Mouse_NA <- read_csv("Galaxy Count Matrix Mapped/Mouse NA/NA_Mouse_Mappings.csv")

# Read the Overall Count Data
counts_Mouse_NA <- read_csv("Galaxy Count Matrix Output/Mouse NA/Mouse_NA_Overall_Counts.csv")

# Change the Column number to the SRR ID for each sample
clean_counts_Mouse_NA <- counts_Mouse_NA %>%
  pivot_longer(cols = starts_with("Column"), names_to = "Dataset") %>%
  mutate(Dataset_Name = {str_extract(Dataset, " \\d{3}:") %>% str_replace(" ", "") %>% str_remove(":") %>% as.numeric()}) %>%
  inner_join(mapping_Mouse_NA, by = c(Dataset_Name = "HISAT")) %>%
  select(Geneid, SRR, value) %>%
  pivot_wider(names_from = SRR)

# Create the data matrix that will be used in PCA plot generation code
clean_counts_Mouse_NA %>% write_csv("Galaxy Count Matrix Mapped/Mouse NA/NA_Mouse_Overall_Mapped_Final.csv")


# Mouse PFC Overall ------------------------------------------------------------------
# Read the Mapping File for SRR to Galaxy Data numbers
mapping_Mouse_PFC <- read_csv("Galaxy Count Matrix Mapped/Mouse PFC/PFC_Mouse_Mappings.csv")

# Read the Overall Count Data
counts_Mouse_PFC <- read_csv("Galaxy Count Matrix Output/Mouse PFC/Mouse_PFC_Overall_Counts.csv")

# Change the Column number to the SRR ID for each sample
clean_counts_Mouse_PFC <- counts_Mouse_PFC %>%
  pivot_longer(cols = starts_with("Column"), names_to = "Dataset") %>%
  mutate(Dataset_Name = {str_extract(Dataset, " \\d{3}:") %>% str_replace(" ", "") %>% str_remove(":") %>% as.numeric()}) %>%
  inner_join(mapping_Mouse_PFC, by = c(Dataset_Name = "HISAT")) %>%
  select(Geneid, SRR, value) %>%
  pivot_wider(names_from = SRR)

# Create the data matrix that will be used in PCA plot generation code
clean_counts_Mouse_PFC %>% write_csv("Galaxy Count Matrix Mapped/Mouse PFC/PFC_Mouse_Overall_Mapped_Final.csv")

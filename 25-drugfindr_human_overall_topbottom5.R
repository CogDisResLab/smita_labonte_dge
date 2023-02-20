#DrugFindR on TopBottom5

library(tidyverse)
library(drugfindR)
library(org.Hs.eg.db)

dge <- read.csv("results/Human_DGE-L1000_Matrix_Overall.csv")
# NOTE: the generated Overall_drugfindR_topbottom5.csv file in results should be
# ignored because it used the Human_DGE-L1000_Matrix_Overall_topbottom5 file as input,
# which was already filtered and the filter_prop below filtered it more -
# this is not what we want (remember smaller n (genes) = greater # CPs)

# NOTE: This script does not accurately generate top/bottom CPs for top and bottom 5%
# of genes. This script needs to be fixed and DrugFindR needs to be run again

# Run DrugFindR on Anterior Insula -------------------------------------------

#generate drugs for top 5%
drugfind_results_AI <- investigate_signature(dge, "CP",
                                          similarity_threshold = 0.2,
                                          filter_prop = 0.95,
                                          gene_column = "Name_GeneSymbol",
                                          logfc_column = "AI",
                                          pval_column = NA,
                                          source_name = "AI_Overall_Up")
#generate drugs for bottom 5%
drugfind_results_AI2 <- investigate_signature(dge, "CP",
                                            similarity_threshold = 0.2,
                                            filter_prop = 0.95,
                                            discordant = TRUE,
                                            gene_column = "Name_GeneSymbol",
                                            logfc_column = "AI",
                                            pval_column = NA,
                                            source_name = "AI_Overall_Dn")

#combine top and bottom drugs
drugfindAI <- bind_rows(drugfind_results_AI, drugfind_results_AI2)


# Run DrugFindR on Cingulate Gyrus -------------------------------------------

#generate drugs for top 5%
drugfind_results_CG <- investigate_signature(dge, "CP",
                                             similarity_threshold = 0.2,
                                             filter_prop = 0.95,
                                             gene_column = "Name_GeneSymbol",
                                             logfc_column = "CG",
                                             pval_column = NA,
                                             source_name = "CG_Overall_Up")
#generate drugs for bottom 5%
drugfind_results_CG2 <- investigate_signature(dge, "CP",
                                              similarity_threshold = 0.2,
                                              filter_prop = 0.95,
                                              discordant = TRUE,
                                              gene_column = "Name_GeneSymbol",
                                              logfc_column = "CG",
                                              pval_column = NA,
                                              source_name = "CG_Overall_Dn")

#combine top and bottom drugs
drugfindCG <- bind_rows(drugfind_results_CG, drugfind_results_CG2)


# Run DrugFindR on DLPFC -----------------------------------------------------

#generate drugs for top 5%
drugfind_results_DLPFC <- investigate_signature(dge, "CP",
                                             similarity_threshold = 0.2,
                                             filter_prop = 0.95,
                                             gene_column = "Name_GeneSymbol",
                                             logfc_column = "DLPFC",
                                             pval_column = NA,
                                             source_name = "DLPFC_Overall_Up")
#generate drugs for bottom 5%
drugfind_results_DLPFC2 <- investigate_signature(dge, "CP",
                                              similarity_threshold = 0.2,
                                              filter_prop = 0.95,
                                              discordant = TRUE,
                                              gene_column = "Name_GeneSymbol",
                                              logfc_column = "DLPFC",
                                              pval_column = NA,
                                              source_name = "DLPFC_Overall_Dn")

#combine top and bottom drugs
drugfindDLPFC <- bind_rows(drugfind_results_DLPFC, drugfind_results_DLPFC2)


# Run DrugFindR on Nucleus Accumbens -----------------------------------------

#generate drugs for top 5%
drugfind_results_NAA <- investigate_signature(dge, "CP",
                                                similarity_threshold = 0.2,
                                                filter_prop = 0.95,
                                                gene_column = "Name_GeneSymbol",
                                                logfc_column = "NAA",
                                                pval_column = NA,
                                                source_name = "NAA_Overall_Up")
#generate drugs for bottom 5%
drugfind_results_NAA2 <- investigate_signature(dge, "CP",
                                                 similarity_threshold = 0.2,
                                                 filter_prop = 0.95,
                                                 discordant = TRUE,
                                                 gene_column = "Name_GeneSymbol",
                                                 logfc_column = "NAA",
                                                 pval_column = NA,
                                                 source_name = "NAA_Overall_Dn")

#combine top and bottom drugs
drugfindNAA <- bind_rows(drugfind_results_NAA, drugfind_results_NAA2)


# Run DrugFindR on Orbitofrontal Cortex --------------------------------------

#generate drugs for top 5%
drugfind_results_OFC <- investigate_signature(dge, "CP",
                                              similarity_threshold = 0.2,
                                              filter_prop = 0.95,
                                              gene_column = "Name_GeneSymbol",
                                              logfc_column = "OFC",
                                              pval_column = NA,
                                              source_name = "OFC_Overall_Up")
#generate drugs for bottom 5%
drugfind_results_OFC2 <- investigate_signature(dge, "CP",
                                               similarity_threshold = 0.2,
                                               filter_prop = 0.95,
                                               discordant = TRUE,
                                               gene_column = "Name_GeneSymbol",
                                               logfc_column = "OFC",
                                               pval_column = NA,
                                               source_name = "OFC_Overall_Dn")

#combine top and bottom drugs
drugfindOFC <- bind_rows(drugfind_results_OFC, drugfind_results_OFC2)


# Run DrugFindR on Subbiculum ------------------------------------------------

#generate drugs for top 5%
drugfind_results_SUB <- investigate_signature(dge, "CP",
                                              similarity_threshold = 0.2,
                                              filter_prop = 0.95,
                                              gene_column = "Name_GeneSymbol",
                                              logfc_column = "SUB",
                                              pval_column = NA,
                                              source_name = "SUB_Overall_Up")
#generate drugs for bottom 5%
drugfind_results_SUB2 <- investigate_signature(dge, "CP",
                                               similarity_threshold = 0.2,
                                               filter_prop = 0.95,
                                               discordant = TRUE,
                                               gene_column = "Name_GeneSymbol",
                                               logfc_column = "SUB",
                                               pval_column = NA,
                                               source_name = "SUB_Overall_Dn")

#combine top and bottom drugs
drugfindSUB <- bind_rows(drugfind_results_SUB, drugfind_results_SUB2)


# Combine Dataframes ----------------------------------------------------------

#combine all top/bottom drug dataframes from each brain region
drugfindAll <- drugfindAI |>
  full_join(drugfindCG) |>
  full_join(drugfindDLPFC) |>
  full_join(drugfindNAA) |>
  full_join(drugfindOFC) |>
  full_join(drugfindSUB) |>

  write_csv("figures/DrugFindR_Output/Human/Overall_drugfindR_nofilterprop.csv")


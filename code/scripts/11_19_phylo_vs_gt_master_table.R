# ============================================================================
# 11.19 — PHYLO VS GT MASTER COMPARISON TABLE
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================================
# 1. SOURCE HELPER FUNCTIONS
# ============================================================================

source("11_19_helper_functions.R")

# ============================================================================
# 2. SETUP DIRECTORIES
# ============================================================================

OUTPUT_DIR <- "../11_19_output/tables/absolute_relative_differences/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("📂 Output directory: ", OUTPUT_DIR)

# ============================================================================
# 3. LOAD DATA
# ============================================================================

message("📂 Loading data...")

gt_full <- read.csv("../output/11_7_complete_50sims/MASTER_ground_truth_bridging_results.csv")
phylo_full <- read.csv("../output/11_7_complete_50sims/MASTER_phylogenetic_bridging_30_results.csv")

message("✅ Data loaded")

# ============================================================================
# 4. PREPARE DATASETS
# ============================================================================

message("🔧 Preparing datasets...")

# Ground truth
gt_phylosims <- gt_full %>%
  pivot_wider(names_from = method, values_from = value) %>%
  mutate(MSM.MSMW = `MSM+MSMW`, MSMW.MSW = `MSMW+MSW`) %>%
  select(-`MSM+MSMW`, -`MSMW+MSW`) %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(product = MSM.MSMW * MSMW.MSW) %>%
  select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, product)

# Phylogenetic
phylo_full_btn_ntwk <- phylo_full %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(
    MSM.MSMW = phylo_30_MSM.MSMW,
    MSMW.MSW = phylo_30_MSMW.MSW,
    product = MSM.MSMW * MSMW.MSW
  ) %>%
  select(sim_id, p_msmw_w, dataset_id, population, dataset_type, 
         downsample_rate, replicate, MSM.MSMW, MSMW.MSW, product)

message("✅ Datasets prepared")

# ============================================================================
# 5. CREATE MASTER TABLE
# ============================================================================

message("📋 Creating master comparison table...")

master_table <- create_master_summary_table(
  phylo_data = phylo_full_btn_ntwk,
  gt_data = gt_phylosims,
  output_dir = OUTPUT_DIR,
  write_csv = TRUE,
  pop_order = c("all", "detected_only", "symptomatic_only", "symptomatic_men",
                "high_activity_symptomatic_men", "random_subsample_high_activity"),
  pop_labels = c(
    "all" = "All Sequences",
    "detected_only" = "Detected Only",
    "symptomatic_only" = "Symptomatic Only",
    "symptomatic_men" = "Symptomatic Men",
    "high_activity_symptomatic_men" = "High Activity Symptomatic Men",
    "random_subsample_high_activity" = "Random Subsample High Activity"
  ),
  downsample_order = c(100, 50, 25, 10, 5, 4, 3, 2, 1),
  rho_order = c(0.05, 0.25, 0.5, 0.75, 0.95)
)

message("✅ Master table created")

# ============================================================================
# 6. SUMMARY
# ============================================================================

message("\n", strrep("=", 85))
message("✨ PHYLO VS GT MASTER TABLE COMPLETE!")
message(strrep("=", 85))

message("\n📊 TABLE DIMENSIONS:")
message("   Rows: ", nrow(master_table))
message("   Populations: ", n_distinct(master_table$Population))
message("   Downsample rates: ", n_distinct(master_table$Downsample_Percent))
message("   Rho values: ", n_distinct(master_table$Rho))

message("\n📋 PREVIEW (first 15 rows):")
print(head(master_table, 15))

message("\n📁 Saved to: ../11_19_output/tables/phylo_vs_gt_master_table.csv")
message(strrep("=", 85))
message("✅ Script complete!")
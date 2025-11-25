# ============================================================================
# 11.19 — STATISTICAL TESTS: KW & DUNN
# ============================================================================
# Purpose: Run KW & Dunn tests on ground truth & phylo data
# Output: All test results saved to ../11_19_output/tables/
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FSA)
  library(purrr)
  library(broom)
})

# ============================================================================
# 1. SOURCE HELPER FUNCTIONS
# ============================================================================

source("11_19_helper_functions.R")

# ============================================================================
# 2. SETUP DIRECTORIES
# ============================================================================

OUTPUT_DIR <- "../11_19_output/tables/kw_dunn_rho_tests/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("📂 Output directory: ", OUTPUT_DIR)

# ============================================================================
# 3. LOAD DATA
# ============================================================================

message("📂 Loading data...")

# Parameter sweep (all 4 batches)
param_sweep <- bind_rows(
  read.csv("../output/11_4_param_sweep/batch_1/bridging_values.csv") %>% mutate(batch = 1),
  read.csv("../output/11_4_param_sweep/batch_2/bridging_values.csv") %>% mutate(batch = 2),
  read.csv("../output/11_4_param_sweep/batch_3/bridging_values.csv") %>% mutate(batch = 3),
  read.csv("../output/11_4_param_sweep/batch_4/bridging_values.csv") %>% mutate(batch = 4)
)

# Phylogenetic full data
phylo_full <- read.csv("../output/11_7_complete_50sims/MASTER_phylogenetic_bridging_30_results.csv")

message("✅ Data loaded")

# ============================================================================
# 4. PREPARE DATASETS
# ============================================================================

message("🔧 Preparing datasets...")

# Parameter sweep for GT KW tests
param_sweep_btn_ntwk <- param_sweep %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(product = MSM.MSMW * MSMW.MSW) %>%
  select(p_msmw_w, replicate, batch, MSM.MSMW, MSMW.MSW, product) %>%
  mutate(p_msmw_w = factor(p_msmw_w))

# Phylogenetic between-network bridging
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
# 5. GROUND TRUTH: KW & DUNN TESTS (Parameter Sweep)
# ============================================================================

message("📊 Running Ground Truth KW & Dunn tests...")

ps_kw_results <- analyze_kw_dunn(
  df = param_sweep_btn_ntwk,
  measures = c("MSM.MSMW", "MSMW.MSW", "product"),
  group_col = "p_msmw_w",
  adjust = "bh",
  write_csv = TRUE,
  out_prefix = file.path(OUTPUT_DIR, "gt_param_sweep_kw_dunn"),
  eps_display = 1e-300
)

message("✅ Ground Truth tests saved")

# ============================================================================
# 6. PHYLO: KW & DUNN TESTS BY STRATA (RHO VALUES)
# ============================================================================

message("📊 Running Phylo KW & Dunn tests by rho strata...")

phylo_kw <- analyze_kw_dunn_by_strata(
  df = phylo_full_btn_ntwk,
  group_col = "p_msmw_w",
  pop_col = "population",
  ds_col = "downsample_rate",
  measures_full = c("MSM.MSMW", "MSMW.MSW", "product"),
  measures_other = "product",
  adjust_kw = "BH",
  adjust_dunn = "bh",
  write_csv = TRUE,
  out_prefix = file.path(OUTPUT_DIR, "phylo_kw_dunn"),
  eps_display = 1e-300
)

message("✅ Phylo KW tests saved")

# ============================================================================
# 7. SUMMARY
# ============================================================================

message("\n", strrep("=", 70))
message("✨ ALL TESTS COMPLETE!")
message(strrep("=", 70))
message("\n📊 Output files saved:")
message("   • gt_param_sweep_kw_dunn_kw_summary.csv")
message("   • gt_param_sweep_kw_dunn_dunn_summary.csv")
message("   • phylo_kw_dunn_KW_summary.csv")
message("   • phylo_kw_dunn_Dunn_summary.csv")
message("\n📁 Location: ", OUTPUT_DIR)
message(strrep("=", 70))
message("\n✅ Script complete!")
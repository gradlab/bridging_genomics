# ============================================================================
# 11.19 — GROUND TRUTH BRIDGING ANALYSIS & PLOTS
# ============================================================================
# Purpose: Load parameter sweep, create ground truth plots, save summary tables
# Output: Figures & tables to ../11_19_output/
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# ============================================================================
# 1. SOURCE HELPER FUNCTIONS
# ============================================================================

source("11_19_helper_functions.R")

# ============================================================================
# 2. SETUP DIRECTORIES
# ============================================================================

FIGURES_DIR <- "../11_19_output/figures/ground_truth_bridging"
TABLES_DIR <- "../11_19_output/tables/ground_truth_bridging"

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

message("📂 Output directories created:")
message("   Figures: ", FIGURES_DIR)
message("   Tables:  ", TABLES_DIR)

# ============================================================================
# 3. LOAD PARAMETER SWEEP DATA
# ============================================================================

message("\n📂 Loading parameter sweep data...")

param_sweep <- bind_rows(
  read.csv("../output/11_4_param_sweep/batch_1/bridging_values.csv") %>% mutate(batch = 1),
  read.csv("../output/11_4_param_sweep/batch_2/bridging_values.csv") %>% mutate(batch = 2),
  read.csv("../output/11_4_param_sweep/batch_3/bridging_values.csv") %>% mutate(batch = 3),
  read.csv("../output/11_4_param_sweep/batch_4/bridging_values.csv") %>% mutate(batch = 4)
)

message("✅ Parameter sweep loaded: ", nrow(param_sweep), " rows")

# ============================================================================
# 4. PREPARE DATA FOR ANALYSIS & PLOTTING
# ============================================================================

message("\n🔧 Preparing data...")

# Filter to between-network measures only
param_sweep_btn_ntwk <- param_sweep %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(product = MSM.MSMW * MSMW.MSW) %>%
  select(p_msmw_w, replicate, batch, MSM.MSMW, MSMW.MSW, product)

message("✅ Filtered to between network measures: ", nrow(param_sweep_btn_ntwk), " rows")

# ============================================================================
# 5. SAVE DATA TABLES
# ============================================================================

message("\n📋 Creating summary tables...")

# Table 1: Replication count by rho
rep_count <- param_sweep_btn_ntwk %>%
  group_by(p_msmw_w) %>%
  summarize(
    count = n(),
    unique_replicates = n_distinct(replicate),
    unique_batches = n_distinct(batch),
    .groups = "drop"
  )

write_csv(rep_count, file.path(TABLES_DIR, "gt_param_sweep_rep_count.csv"))
message("   ✅ gt_param_sweep_rep_count.csv")

# Table 2: Summary statistics by rho
summary_stats <- param_sweep_btn_ntwk %>%
  mutate(product = MSM.MSMW * MSMW.MSW) %>%
  pivot_longer(cols = c(MSM.MSMW, MSMW.MSW, product), 
               names_to = "measure", values_to = "value") %>%
  group_by(p_msmw_w, measure) %>%
  summarize(
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    measure_label = case_when(
      measure == "MSM.MSMW" ~ "MSMW with MSM",
      measure == "MSMW.MSW" ~ "MSMW with MSW",
      measure == "product" ~ "Product (Bridge Index)",
      TRUE ~ measure
    )
  ) %>%
  select(p_msmw_w, measure_label, everything(), -measure) %>%
  arrange(measure_label, p_msmw_w)

write_csv(summary_stats, file.path(TABLES_DIR, "gt_param_sweep_summary_stats.csv"))
message("   ✅ gt_param_sweep_summary_stats.csv")

# ============================================================================
# 6. CREATE GROUND TRUTH PLOTS
# ============================================================================

message("\n🎨 Creating ground truth plots...")

# Reformat data for plotting (add measure column back for compatibility)
param_sweep_for_plots <- param_sweep %>%
  filter(measure == "total_between_network_proportion")

# Create all three plots
plots <- create_bridging_plots_single(
  data = param_sweep_for_plots,
  output_dir = FIGURES_DIR,
  palette_name = "bold_scientific_rainbow",
  # file naming
  scatter_filename = "01_gt_bridging_scatter.tiff",
  summary_filename = "02_gt_bridging_summary.tiff",
  combined_filename = "03_gt_bridging_combined.tiff",
  # dimensions
  width1 = 7, height1 = 5,
  width2 = 12, height2 = 5,
  width_grid = 12, height_grid = 5,
  dpi = 300,
  # labels & titles
  main_title = "Ground Truth: Between Network Bridging Analysis",
  scatter_xlab = "MSMW with MSM Bridging Proportion",
  scatter_ylab = "MSMW with MSW Bridging Proportion",
  summary_title = "Between Network Bridging Metrics by ρ",
  summary_subtitle = "Ground Truth | Median ± Interquartile Range",
  summary_xlab = bquote(rho ~ " (probability MSMW seek MSW partner)"),
  summary_ylab = "Bridging Value"
)

message("✅ Plots created:")
message("   • Scatter plot")
message("   • Summary plot (faceted)")
message("   • Combined plot")

# ============================================================================
# 7. SAVE INDIVIDUAL PLOT COMPONENTS (optional)
# ============================================================================

message("\n💾 Saving individual components...")

# Save scatter plot at different size
ggsave(file.path(FIGURES_DIR, "01_gt_bridging_scatter_large.tiff"),
       plot = plots$scatter, width = 10, height = 7, dpi = 300, bg = "white")

# Save summary plot at different size
ggsave(file.path(FIGURES_DIR, "02_gt_bridging_summary_large.tiff"),
       plot = plots$summary, width = 14, height = 6, dpi = 300, bg = "white")

message("✅ Additional sizes saved")

# ============================================================================
# 8. CREATE DETAILED BREAKDOWN TABLES
# ============================================================================

message("\n📊 Creating detailed breakdown tables...")

# Table 3: By batch
by_batch <- param_sweep_btn_ntwk %>%
  pivot_longer(cols = c(MSM.MSMW, MSMW.MSW, product),
               names_to = "measure", values_to = "value") %>%
  group_by(batch, p_msmw_w, measure) %>%
  summarize(
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(by_batch, file.path(TABLES_DIR, "gt_param_sweep_by_batch.csv"))
message("   ✅ gt_param_sweep_by_batch.csv")

# Table 4: All raw values (for external analysis)
raw_values <- param_sweep_btn_ntwk %>%
  pivot_longer(cols = c(MSM.MSMW, MSMW.MSW, product),
               names_to = "measure", values_to = "value") %>%
  arrange(p_msmw_w, measure, replicate)

write_csv(raw_values, file.path(TABLES_DIR, "gt_param_sweep_raw_values.csv"))
message("   ✅ gt_param_sweep_raw_values.csv")

# ============================================================================
# 9. SUMMARY REPORT
# ============================================================================

message("\n" %+% strrep("=", 75))
message("✨ GROUND TRUTH PLOTS & ANALYSIS COMPLETE!")
message(strrep("=", 75))

message("\n📊 SUMMARY STATISTICS:")
print(summary_stats)

message("\n📈 REPLICATION COUNT:")
print(rep_count)

message("\n📁 OUTPUT LOCATIONS:")
message(strrep("-", 75))
message("Figures saved to:")
message("  ", FIGURES_DIR)
message("\nTables saved to:")
message("  ", TABLES_DIR)
message("\nFiles created:")
message("  ✅ 01_gt_bridging_scatter.tiff")
message("  ✅ 01_gt_bridging_scatter_large.tiff")
message("  ✅ 02_gt_bridging_summary.tiff")
message("  ✅ 02_gt_bridging_summary_large.tiff")
message("  ✅ 03_gt_bridging_combined.tiff")
message("  ✅ gt_param_sweep_rep_count.csv")
message("  ✅ gt_param_sweep_summary_stats.csv")
message("  ✅ gt_param_sweep_by_batch.csv")
message("  ✅ gt_param_sweep_raw_values.csv")

message(strrep("=", 75))
message("✅ Script complete! All outputs ready for analysis.")
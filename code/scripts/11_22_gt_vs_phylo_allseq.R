# ============================================================================
# 11.22 — GROUND TRUTH VS PHYLO (ALL SEQUENCES) COMPARISON
# ============================================================================
# Purpose: Compare GT to phylo using full sequences, create clean plot
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# ============================================================================
# 1. SETUP DIRECTORIES
# ============================================================================

FIGURES_DIR <- "../11_22_output/figures/gt_vs_phylo_allseq"
TABLES_DIR <- "../11_22_output/tables/gt_vs_phylo_allseq"

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

message("📂 Output directories:")
message("   Figures: ", FIGURES_DIR)
message("   Tables:  ", TABLES_DIR)

# ============================================================================
# 2. LOAD DATA
# ============================================================================

message("\n📂 Loading data...")

gt_full <- read.csv("../output/11_7_complete_50sims/MASTER_ground_truth_bridging_results.csv")
phylo_full <- read.csv("../output/11_7_complete_50sims/MASTER_phylogenetic_bridging_30_results.csv")

message("✅ Data loaded")

# ============================================================================
# 3. PREPARE DATASETS
# ============================================================================

message("\n🔧 Preparing datasets...")

# Ground truth between-network bridging
gt_phylosims <- gt_full %>%
  pivot_wider(names_from = method, values_from = value) %>%
  mutate(MSM.MSMW = `MSM+MSMW`, MSMW.MSW = `MSMW+MSW`) %>%
  select(-`MSM+MSMW`, -`MSMW+MSW`) %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(product = MSM.MSMW * MSMW.MSW) %>%
  select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, product)

# Phylo data for full sequences
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
# 4. PREPARE PLOTTING DATA
# ============================================================================

message("\n📋 Preparing plotting data...")

# Phylo data: full sequences only, complete dataset (100%)
phylo_fullseq <- phylo_full_btn_ntwk %>%
  filter(dataset_id == "full_sequences") %>%
  mutate(method = "phylo_all_infections")

# Combine GT and phylo
plot_df_allseqs <- bind_rows(
  gt_phylosims %>% mutate(method = "ground_truth"),
  phylo_fullseq %>% select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, product, method)
) %>%
  # Compute medians and IQR per rho & method
  pivot_longer(cols = c(MSM.MSMW, MSMW.MSW, product),
               names_to = "name", values_to = "value") %>%
  group_by(p_msmw_w, name, method) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    n_sims = n_distinct(sim_id),
    .groups = "drop"
  ) %>%
  # Offset x-values for visual separation
  mutate(
    p_msmw_w_offset = case_when(
      method == "ground_truth" ~ p_msmw_w - 0.01,
      TRUE ~ p_msmw_w + 0.01
    ),
    method_label = case_when(
      method == "ground_truth" ~ "Ground Truth",
      method == "phylo_all_infections" ~ "Sequences From All Infections",
      TRUE ~ method
    )
  )

message("✅ Plotting data prepared")

# ============================================================================
# 5. CREATE PLOT
# ============================================================================

message("\n🎨 Creating plot...")

comparison_plot <- plot_df_allseqs %>%
  ggplot(aes(x = p_msmw_w_offset, y = median_value, color = method_label)) +
  # Points
  geom_point() +
  # Error bars (IQR)
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.015, linewidth = 0.8) +
  # Color scale
  scale_color_manual(
    values = c(
      "Ground Truth" = "#D73027",
      "Sequences From All Infections" = "#4DDADA"
    ),
    name = "Base Population"
  ) +
  # X-axis for rho values
  scale_x_continuous(
    breaks = c(0.05, 0.25, 0.5, 0.75, 0.95),
    labels = c("0.05", "0.25", "0.50", "0.75", "0.95")
  ) +
  # Facets for measures
  facet_wrap(~factor(name,
                     levels = c("MSM.MSMW", "MSMW.MSW", "product"),
                     labels = c("Bridging Prop. 1", "Bridging Prop. 2", "Bridging Metric")),
             scales = "free_y", nrow = 1) +
  # Theme
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    panel.spacing = unit(0.8, "lines"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Ground Truth vs Phylodynamic Estimates (All Sequences)",
    x = bquote(rho ~ " - probability MSMW seek WSM partner"),
    y = "Bridging Value"
  )

message("✅ Plot created")

# ============================================================================
# 6. SAVE PLOT
# ============================================================================

message("\n💾 Saving plots...")

# Save as TIFF
ggsave(
  file.path(FIGURES_DIR, "00_allseq_3panel_comparison.tiff"),
  plot = comparison_plot,
  width = 9, height = 5, dpi = 300, bg = "white"
)

# Save as PNG
ggsave(
  file.path(FIGURES_DIR, "00_allseq_3panel_comparison.png"),
  plot = comparison_plot,
  width = 9, height = 5, dpi = 300, bg = "white"
)

# Save as RDS object for combining with other figures later
saveRDS(comparison_plot, file.path(FIGURES_DIR, "plot_3panel_allseq.rds"))

message("   ✅ 00_allseq_3panel_comparison.tiff")
message("   ✅ 00_allseq_3panel_comparison.png")
message("   ✅ plot_3panel_allseq.rds")

message("✅ All plot formats saved")

# ============================================================================
# 7. SAVE COMPARISON TABLES
# ============================================================================

message("\n📋 Saving comparison tables...")

# Table 1: Summary data for plotting
write_csv(plot_df_allseqs,
          file.path(TABLES_DIR, "allseq_comparison_plot_data.csv"))
message("   ✅ allseq_comparison_plot_data.csv")

# Table 2: Detailed differences
allseq_differences <- bind_rows(
  gt_phylosims %>% mutate(method = "ground_truth"),
  phylo_fullseq %>% select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, product, method)
) %>%
  pivot_longer(cols = c(MSM.MSMW, MSMW.MSW, product),
               names_to = "measure", values_to = "value") %>%
  pivot_wider(names_from = method, values_from = value) %>%
  mutate(
    difference = phylo_all_infections - ground_truth,
    rel_diff_to_gt = difference / ground_truth,
    rel_diff_to_est = difference / phylo_all_infections,
    measure_label = case_when(
      measure == "MSM.MSMW" ~ "MSMW with MSM",
      measure == "MSMW.MSW" ~ "MSMW with MSW",
      measure == "product" ~ "Product",
      TRUE ~ measure
    )
  ) %>%
  select(sim_id, p_msmw_w, measure_label, ground_truth, phylo_all_infections,
         difference, rel_diff_to_gt, rel_diff_to_est)

write_csv(allseq_differences,
          file.path(TABLES_DIR, "allseq_comparison_raw_differences.csv"))
message("   ✅ allseq_comparison_raw_differences.csv")

# Table 3: Difference medians
allseq_differences_summary <- allseq_differences %>%
  group_by(p_msmw_w, measure_label) %>%
  summarize(
    median_diff = median(difference, na.rm = TRUE),
    q25_diff = quantile(difference, 0.25, na.rm = TRUE),
    q75_diff = quantile(difference, 0.75, na.rm = TRUE),
    median_rel_to_gt = median(rel_diff_to_gt, na.rm = TRUE),
    q25_rel_to_gt = quantile(rel_diff_to_gt, 0.25, na.rm = TRUE),
    q75_rel_to_gt = quantile(rel_diff_to_gt, 0.75, na.rm = TRUE),
    median_rel_to_est = median(rel_diff_to_est, na.rm = TRUE),
    q25_rel_to_est = quantile(rel_diff_to_est, 0.25, na.rm = TRUE),
    q75_rel_to_est = quantile(rel_diff_to_est, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(allseq_differences_summary,
          file.path(TABLES_DIR, "allseq_comparison_differences_summary.csv"))
message("   ✅ allseq_comparison_differences_summary.csv")

# ============================================================================
# 8. SUMMARY
# ============================================================================

message("\n", strrep("=", 80))
message("✨ GT VS PHYLO (ALL SEQUENCES) ANALYSIS COMPLETE!")
message(strrep("=", 80))

message("\n📊 COMPARISON SUMMARY:")
print(allseq_differences_summary)

message("\n📁 OUTPUT LOCATIONS:")
message(strrep("-", 80))
message("Figures: ", FIGURES_DIR)
message("  ✅ 00_allseq_3panel_comparison.tiff")
message("  ✅ 00_allseq_3panel_comparison.png")
message("  ✅ plot_3panel_allseq.rds")

message("\nTables: ", TABLES_DIR)
message("  ✅ allseq_comparison_plot_data.csv")
message("  ✅ allseq_comparison_raw_differences.csv")
message("  ✅ allseq_comparison_differences_summary.csv")

message(strrep("=", 80))
message("✅ Script complete!")

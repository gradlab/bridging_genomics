# ============================================================================
# 11.19 — GROUND TRUTH VS PHYLO (COMPLETE POPULATIONS) COMPARISON
# ============================================================================
# Purpose: Compare GT to all phylo complete populations, create clean plot
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# ============================================================================
# 1. SETUP DIRECTORIES
# ============================================================================

FIGURES_DIR <- "../11_19_output/figures/gt_vs_phylo_complete_pops"
TABLES_DIR <- "../11_19_output/tables/gt_vs_phylo_complete_pops"

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

# Phylo data for complete populations (100% downsample)
phylo_complete_datasets <- phylo_full %>%
  filter(measure == "total_between_network_proportion", downsample_rate == 100) %>%
  mutate(
    MSM.MSMW = phylo_30_MSM.MSMW,
    MSMW.MSW = phylo_30_MSMW.MSW,
    method = dataset_id
  ) %>%
  select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, method)

# Ground truth between-network bridging
gt_phylosims <- gt_full %>%
  pivot_wider(names_from = method, values_from = value) %>%
  mutate(MSM.MSMW = `MSM+MSMW`, MSMW.MSW = `MSMW+MSW`) %>%
  select(-`MSM+MSMW`, -`MSMW+MSW`) %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(method = "ground_truth") %>%
  select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, method)

# Combine and compute product
complete_datasets_combined <- bind_rows(gt_phylosims, phylo_complete_datasets) %>%
  mutate(product = MSM.MSMW * MSMW.MSW) %>%
  pivot_longer(cols = c(MSM.MSMW, MSMW.MSW, product),
               names_to = "name", values_to = "value")

message("✅ Datasets prepared")

# ============================================================================
# 4. PREPARE PLOTTING DATA (PRODUCT ONLY)
# ============================================================================

message("\n📋 Preparing plotting data...")

plot_df_complete <- complete_datasets_combined %>%
  filter(name == "product") %>%
  group_by(p_msmw_w, name, method) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    n_sims = n_distinct(sim_id),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method, levels = c(
      "ground_truth", "full_sequences", "detected_only_complete",
      "symptomatic_only_complete", "symptomatic_men_complete",
      "high_activity_symptomatic_men_complete", "random_subsample_high_activity_complete"
    )),
    facet_label = paste0("rho==", p_msmw_w)
  )

message("✅ Plotting data prepared")

# ============================================================================
# 5. PREPARE COLOR PALETTE
# ============================================================================

color_palette <- c(
  "ground_truth" = "#D73027",
  "full_sequences" = "#4DDADA",
  "detected_only_complete" = "#52B88F",
  "symptomatic_only_complete" = "#3FA27A",
  "symptomatic_men_complete" = "#2C7A78",
  "high_activity_symptomatic_men_complete" = "#164B53",
  "random_subsample_high_activity_complete" = "#2c7bb6"
)

# ============================================================================
# 6. CREATE PLOT
# ============================================================================

message("\n🎨 Creating plot...")

comparison_plot <- plot_df_complete %>%
  ggplot(aes(x = method, y = median_value, color = method, fill = method)) +
  # Points
  geom_point() +
  # Error bars (IQR)
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.2, linewidth = 0.8) +
  # Color scale
  scale_color_manual(
    values = color_palette,
    labels = c(
      "ground_truth" = "Ground Truth",
      "full_sequences" = "Phylo - All Infections",
      "detected_only_complete" = "Phylo - Detected Only",
      "symptomatic_only_complete" = "Phylo - Symptomatic Only",
      "symptomatic_men_complete" = "Phylo - Symptomatic Men",
      "high_activity_symptomatic_men_complete" = "Phylo - High Activity Symptomatic Men",
      "random_subsample_high_activity_complete" = "Phylo - Random Subsample High Activity"
    ),
    name = "Method"
  ) +
  scale_fill_manual(
    values = color_palette,
    guide = "none"
  ) +
  # Facets with parsed labels for Greek rho
  facet_wrap(~facet_label, nrow = 1, labeller = label_parsed) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.spacing.x = unit(1.2, "lines"),
    panel.border = element_blank()
  ) +
  labs(
    x = "",
    y = "Bridging Value (Product)",
    title = "Complete Populations: Ground Truth vs Phylodynamic Estimates"
  )

message("✅ Plot created")

# ============================================================================
# 7. SAVE PLOT
# ============================================================================

message("\n💾 Saving plot...")

ggsave(
  file.path(FIGURES_DIR, "complete_pops_comparison.tiff"),
  plot = comparison_plot,
  width = 11, height = 6, dpi = 300, bg = "white"
)
ggsave(
  file.path(FIGURES_DIR, "complete_pops_comparison.png"),
  plot = comparison_plot,
  width = 11, height = 6, dpi = 300, bg = "white"
)

message("✅ Plot saved")

# ============================================================================
# 8. SAVE COMPARISON TABLES
# ============================================================================

message("\n📋 Saving comparison tables...")

# Table 1: Summary data for plotting (all three measures)
complete_datasets_combined_summary <- complete_datasets_combined %>%
  group_by(p_msmw_w, name, method) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    n_sims = n_distinct(sim_id),
    .groups = "drop"
  )

write_csv(complete_datasets_combined_summary,
          file.path(TABLES_DIR, "complete_pops_comparison_plot_data.csv"))
message("   ✅ complete_pops_comparison_plot_data.csv")

# Table 2: Detailed differences (product only)
complete_diffs <- complete_datasets_combined %>%
  filter(name == "product") %>%
  pivot_wider(names_from = method, values_from = value) %>%
  pivot_longer(
    cols = -c(sim_id, p_msmw_w, name, ground_truth),
    names_to = "phylo_method", 
    values_to = "phylo_value"
  ) %>%
  mutate(
    difference = phylo_value - ground_truth,
    rel_diff_to_gt = difference / ground_truth,
    rel_diff_to_est = difference / phylo_value,
    method_label = case_when(
      phylo_method == "full_sequences" ~ "Phylo - All Infections",
      phylo_method == "detected_only_complete" ~ "Phylo - Detected Only",
      phylo_method == "symptomatic_only_complete" ~ "Phylo - Symptomatic Only",
      phylo_method == "symptomatic_men_complete" ~ "Phylo - Symptomatic Men",
      phylo_method == "high_activity_symptomatic_men_complete" ~ "Phylo - High Activity Symptomatic Men",
      phylo_method == "random_subsample_high_activity_complete" ~ "Phylo - Random Subsample High Activity",
      TRUE ~ phylo_method
    )
  ) %>%
  select(sim_id, p_msmw_w, method_label, ground_truth, phylo_value,
         difference, rel_diff_to_gt, rel_diff_to_est)

write_csv(complete_diffs,
          file.path(TABLES_DIR, "complete_pops_comparison_raw_differences.csv"))
message("   ✅ complete_pops_comparison_raw_differences.csv")

# Table 3: Difference medians (product only)
complete_diffs_summary <- complete_diffs %>%
  group_by(p_msmw_w, method_label) %>%
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

write_csv(complete_diffs_summary,
          file.path(TABLES_DIR, "complete_pops_comparison_differences_summary.csv"))
message("   ✅ complete_pops_comparison_differences_summary.csv")

# ============================================================================
# 9. SUMMARY
# ============================================================================

message("\n", strrep("=", 80))
message("✨ GT VS PHYLO (COMPLETE POPULATIONS) ANALYSIS COMPLETE!")
message(strrep("=", 80))

message("\n📊 COMPARISON SUMMARY (Product measure):")
print(complete_diffs_summary)

message("\n📁 OUTPUT LOCATIONS:")
message(strrep("-", 80))
message("Figure: ", FIGURES_DIR)
message("  ✅ complete_pops_comparison.tiff")
message("  ✅ complete_pops_comparison.png")

message("\nTables: ", TABLES_DIR)
message("  ✅ complete_pops_comparison_plot_data.csv")
message("  ✅ complete_pops_comparison_raw_differences.csv")
message("  ✅ complete_pops_comparison_differences_summary.csv")

message(strrep("=", 80))
message("✅ Script complete!")

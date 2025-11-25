# ============================================================================
# 11.19 — GROUND TRUTH VS PHYLO (DOWNSAMPLES) COMPARISON
# ============================================================================
# Purpose: Create unified and side-by-side downsampling plots (no stars)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# ============================================================================
# 1. SETUP DIRECTORIES
# ============================================================================

FIGURES_DIR <- "../11_19_output/figures/gt_vs_phylo_downsamples"
TABLES_DIR <- "../11_19_output/tables/gt_vs_phylo_downsamples"

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

# Ground truth
gt_phylosims <- gt_full %>%
  pivot_wider(names_from = method, values_from = value) %>%
  mutate(MSM.MSMW = `MSM+MSMW`, MSMW.MSW = `MSMW+MSW`) %>%
  select(-`MSM+MSMW`, -`MSMW+MSW`) %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(product = MSM.MSMW * MSMW.MSW) %>%
  select(sim_id, p_msmw_w, product)

# Phylo data with downsampling
phylo_full_btn_ntwk <- phylo_full %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(
    MSM.MSMW = phylo_30_MSM.MSMW,
    MSMW.MSW = phylo_30_MSMW.MSW,
    product = MSM.MSMW * MSMW.MSW
  ) %>%
  select(sim_id, p_msmw_w, population, dataset_type, downsample_rate, replicate, product)

message("✅ Datasets prepared")

# ============================================================================
# 4. PREPARE PLOTTING DATA
# ============================================================================

message("\n📋 Preparing plotting data...")

# Pretty labels
pop_to_pretty <- c(
  "all" = "Phylo - All Infections",
  "detected_only" = "Phylo - Detected Only",
  "symptomatic_only" = "Phylo - Symptomatic Only",
  "symptomatic_men" = "Phylo - Symptomatic Men",
  "high_activity_symptomatic_men" = "Phylo - High Activity Symptomatic Men",
  "random_subsample_high_activity" = "Phylo - Random Subsample High Activity"
)

downsample_levels <- c(100, 50, 25, 10, 5, 4, 3, 2, 1)
method_levels <- unname(pop_to_pretty)

# Prepare downsample plot data
downsamples_plot_data <- phylo_full_btn_ntwk %>%
  mutate(
    population_clean = recode(population, !!!pop_to_pretty),
    downsample_rate_clean = case_when(
      dataset_type == "complete" ~ 100,
      population == "all" & dataset_type == "downsampled" ~ 100,
      TRUE ~ downsample_rate
    )
  ) %>%
  group_by(sim_id, p_msmw_w, population_clean, downsample_rate_clean) %>%
  summarise(product_sim_median = median(product, na.rm = TRUE), .groups = "drop") %>%
  group_by(p_msmw_w, population_clean, downsample_rate_clean) %>%
  summarise(
    median_value = median(product_sim_median, na.rm = TRUE),
    q25 = quantile(product_sim_median, 0.25, na.rm = TRUE),
    q75 = quantile(product_sim_median, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method_factor = factor(population_clean, levels = method_levels),
    downsample_factor = factor(downsample_rate_clean, levels = downsample_levels),
    facet_label = paste0("rho==", p_msmw_w)
  )

# GT data expanded
gt_data_expanded <- gt_phylosims %>%
  group_by(p_msmw_w) %>%
  summarise(
    median_value = median(product),
    q25 = quantile(product, 0.25),
    q75 = quantile(product, 0.75),
    .groups = "drop"
  ) %>%
  slice(rep(1:n(), each = length(downsample_levels))) %>%
  mutate(
    downsample_rate_clean = rep(downsample_levels, times = n()/length(downsample_levels)),
    downsample_factor = factor(downsample_rate_clean, levels = downsample_levels),
    facet_label = paste0("rho==", p_msmw_w)
  )

message("✅ Plotting data prepared")

# ============================================================================
# 5. BUILD BASE PLOT (UNIFIED - ALL DOWNSAMPLES)
# ============================================================================

message("\n🎨 Creating unified plot...")

plot_unified <- ggplot(
  downsamples_plot_data,
  aes(x = downsample_factor, y = median_value, color = method_factor)
) +
  # GT ribbon & line
  geom_ribbon(
    data = gt_data_expanded,
    aes(x = downsample_factor, ymin = q25, ymax = q75, group = p_msmw_w),
    alpha = 0.35, fill = "#D73027", color = NA, inherit.aes = FALSE
  ) +
  geom_line(
    data = gt_data_expanded,
    aes(x = downsample_factor, y = median_value, group = p_msmw_w),
    color = "#D73027", linewidth = 1, inherit.aes = FALSE
  ) +
  # Phylo points + IQR
  geom_point(alpha = 0.8, position = position_dodge(width = 0.6)) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.2, linewidth = 0.8, alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~facet_label, nrow = 1, labeller = label_parsed) +
  scale_x_discrete(labels = paste0(downsample_levels, "%")) +
  scale_color_manual(
    values = c(
      "Phylo - All Infections" = "#4DDADA",
      "Phylo - Detected Only" = "#52B88F",
      "Phylo - Symptomatic Only" = "#3FA27A",
      "Phylo - Symptomatic Men" = "#2C7A78",
      "Phylo - High Activity Symptomatic Men" = "#164B53",
      "Phylo - Random Subsample High Activity" = "#2c7bb6"
    ),
    name = "Method",
    limits = method_levels,
    breaks = method_levels
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = rep("blank", length(method_levels)),
        shape = rep(16, length(method_levels)),
        alpha = rep(0.8, length(method_levels)),
        size = rep(2, length(method_levels))
      )
    )
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    panel.spacing.x = unit(1.2, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.6, "lines"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Downsample Rate (%)",
    y = "Bridging Value (Product)",
    title = "All Downsampling Rates"
  )

message("✅ Unified plot created")

# ============================================================================
# 6. BUILD SIDE-BY-SIDE FUNCTION
# ============================================================================

message("\n🎨 Creating side-by-side plots...")

build_panel <- function(levels_vec, panel_title = "") {
  lev <- as.numeric(levels_vec)

  # Filter phylo data by downsample_rate_clean
  phylo_dat <- downsamples_plot_data %>%
    filter(downsample_rate_clean %in% lev) %>%
    mutate(
      downsample_factor = factor(downsample_rate_clean, levels = lev),
      facet_label = paste0("rho==", p_msmw_w)
    )

  # Filter GT data
  gt_dat <- gt_data_expanded %>%
    filter(downsample_rate_clean %in% lev) %>%
    mutate(downsample_factor = factor(downsample_rate_clean, levels = lev))

  # Build plot
  p <- ggplot(phylo_dat, aes(x = downsample_factor, y = median_value, color = method_factor)) +
    # GT ribbon & line
    geom_ribbon(
      data = gt_dat,
      aes(x = downsample_factor, ymin = q25, ymax = q75, group = p_msmw_w),
      alpha = 0.35, fill = "#D73027", color = NA, inherit.aes = FALSE
    ) +
    geom_line(
      data = gt_dat,
      aes(x = downsample_factor, y = median_value, group = p_msmw_w),
      color = "#D73027", linewidth = 1, inherit.aes = FALSE
    ) +
    # Phylo points + IQR
    geom_point(alpha = 0.8, position = position_dodge(width = 0.6)) +
    geom_errorbar(
      aes(ymin = q25, ymax = q75),
      width = 0.2, linewidth = 0.8, alpha = 0.8,
      position = position_dodge(width = 0.6)
    ) +
    facet_wrap(~facet_label, nrow = 1, labeller = label_parsed) +
    scale_x_discrete(labels = paste0(lev, "%")) +
    scale_color_manual(
      values = c(
        "Phylo - All Infections" = "#4DDADA",
        "Phylo - Detected Only" = "#52B88F",
        "Phylo - Symptomatic Only" = "#3FA27A",
        "Phylo - Symptomatic Men" = "#2C7A78",
        "Phylo - High Activity Symptomatic Men" = "#164B53",
        "Phylo - Random Subsample High Activity" = "#2c7bb6"
      ),
      name = "Method",
      limits = method_levels,
      breaks = method_levels
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          linetype = rep("blank", length(method_levels)),
          shape = rep(16, length(method_levels)),
          alpha = rep(0.8, length(method_levels)),
          size = rep(2, length(method_levels))
        )
      )
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5, "pt"),
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
      strip.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      legend.position = "right",
      panel.spacing.x = unit(1.2, "lines"),
      plot.background = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.6, "lines")
    ) +
    labs(x = "Downsample Rate (%)", y = "Bridging Value (Product)", title = panel_title)

  p
}

# Constants
upper_levels <- c(100, 50, 25, 10)
lower_levels <- c(5, 4, 3, 2, 1)

# Build side-by-side panels
left_panel <- build_panel(upper_levels, "100–10%")
right_panel <- build_panel(lower_levels, "5–1%")

# Extract legend
leg <- cowplot::get_legend(left_panel)

# Combine into side-by-side
plot_sidebyside <- cowplot::plot_grid(
  cowplot::plot_grid(
    left_panel + theme(legend.position = "none"),
    right_panel + theme(legend.position = "none"),
    ncol = 2, rel_widths = c(1, 1)
  ),
  leg, ncol = 2, rel_widths = c(1, 0.25)
)

message("✅ Side-by-side plot created")

# ============================================================================
# 7. SAVE PLOTS
# ============================================================================

message("\n💾 Saving plots...")

# Unified plot
ggsave(file.path(FIGURES_DIR, "01_downsamples_unified.tiff"),
       plot = plot_unified, width = 11, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "01_downsamples_unified.png"),
       plot = plot_unified, width = 11, height = 5, dpi = 300, bg = "white")
message("   ✅ 01_downsamples_unified")

# Side-by-side plot
ggsave(file.path(FIGURES_DIR, "02_downsamples_sidebyside.tiff"),
       plot = plot_sidebyside, width = 14, height = 6, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "02_downsamples_sidebyside.png"),
       plot = plot_sidebyside, width = 14, height = 6, dpi = 300, bg = "white")
message("   ✅ 02_downsamples_sidebyside")

# ============================================================================
# 8. SAVE TABLES
# ============================================================================

message("\n📋 Saving tables...")

# Table 1: Plot data
write_csv(downsamples_plot_data %>% 
            select(-method_factor, -downsample_factor, -facet_label),
          file.path(TABLES_DIR, "downsamples_plot_data.csv"))
message("   ✅ downsamples_plot_data.csv")

# Table 2: GT data
write_csv(gt_data_expanded %>%
            select(-downsample_factor, -facet_label),
          file.path(TABLES_DIR, "downsamples_gt_data.csv"))
message("   ✅ downsamples_gt_data.csv")

# ============================================================================
# 9. SUMMARY
# ============================================================================

message("\n", strrep("=", 80))
message("✨ DOWNSAMPLING ANALYSIS COMPLETE!")
message(strrep("=", 80))

message("\n📊 PLOTS CREATED:")
message("   • 01_downsamples_unified (all downsampling rates together)")
message("   • 02_downsamples_sidebyside (100-10% left, 5-1% right)")

message("\n📁 OUTPUT LOCATIONS:")
message(strrep("-", 80))
message("Figures: ", FIGURES_DIR)
message("  ✅ 01_downsamples_unified.tiff/.png")
message("  ✅ 02_downsamples_sidebyside.tiff/.png")

message("\nTables: ", TABLES_DIR)
message("  ✅ downsamples_plot_data.csv")
message("  ✅ downsamples_gt_data.csv")

message(strrep("=", 80))
message("✅ Script complete!")
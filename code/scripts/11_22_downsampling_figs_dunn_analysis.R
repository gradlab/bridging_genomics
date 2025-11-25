# ============================================================================
# 11.22 — DOWNSAMPLING & DUNN ANALYSIS WITH COMBINED MIDDLE FACET
# ============================================================================
# Purpose: Create all downsampling plots + Dunn heatmap, then extract rho=0.5
#          panels for combined A+B figure (Panel A narrower + colored y-axis in B)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggtext)
})

# ============================================================================
# 1. SETUP DIRECTORIES
# ============================================================================

FIGURES_DIR <- "../11_22_output/figures/downsampling_dunn_combined"
TABLES_DIR <- "../11_22_output/tables/downsampling_dunn_combined"

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
dunn_full <- read_csv("../11_19_output/tables/downsampling_kw_dunn/downsampling_dunn_results.csv")

message("✅ Data loaded")

# ============================================================================
# 3. PREPARE DOWNSAMPLING DATASETS
# ============================================================================

message("\n🔧 Preparing downsampling datasets...")

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

message("✅ Downsampling datasets prepared")

# ============================================================================
# 4. PREPARE DOWNSAMPLING PLOTTING DATA
# ============================================================================

message("\n📋 Preparing downsampling plotting data...")

# Pretty labels
pop_to_pretty <- c(
  "all" = "Phylo - All Infections",
  "detected_only" = "Phylo - Detected Only",
  "symptomatic_only" = "Phylo - Symptomatic Only",
  "symptomatic_men" = "Phylo - Symptomatic Men",
  "high_activity_symptomatic_men" = "Phylo - High Activity Symptomatic Men",
  "random_subsample_high_activity" = "Phylo - Random Subsample High Activity"
)

# Color mapping
pop_colors <- c(
  "Phylo - All Infections" = "#4DDADA",
  "Phylo - Detected Only" = "#52B88F",
  "Phylo - Symptomatic Only" = "#3FA27A",
  "Phylo - Symptomatic Men" = "#2C7A78",
  "Phylo - High Activity Symptomatic Men" = "#164B53",
  "Phylo - Random Subsample High Activity" = "#2c7bb6"
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

message("✅ Downsampling plotting data prepared")

# ============================================================================
# 5. BUILD DOWNSAMPLING PLOTS
# ============================================================================

message("\n🎨 Creating downsampling plots...")

# Unified plot
plot_unified <- ggplot(
  downsamples_plot_data,
  aes(x = downsample_factor, y = median_value, color = method_factor)
) +
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
  geom_point(alpha = 0.8, position = position_dodge(width = 0.6)) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.2, linewidth = 0.8, alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~facet_label, nrow = 1, labeller = label_parsed) +
  scale_x_discrete(labels = paste0(downsample_levels, "%")) +
  scale_color_manual(
    values = pop_colors,
    name = "Base Population",
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
    y = "Bridging Metric",
    title = "All Downsampling Rates"
  )

message("✅ Unified plot created")

# Side-by-side function
build_panel <- function(levels_vec, panel_title = "") {
  lev <- as.numeric(levels_vec)

  phylo_dat <- downsamples_plot_data %>%
    filter(downsample_rate_clean %in% lev) %>%
    mutate(
      downsample_factor = factor(downsample_rate_clean, levels = lev),
      facet_label = paste0("rho==", p_msmw_w)
    )

  gt_dat <- gt_data_expanded %>%
    filter(downsample_rate_clean %in% lev) %>%
    mutate(downsample_factor = factor(downsample_rate_clean, levels = lev))

  p <- ggplot(phylo_dat, aes(x = downsample_factor, y = median_value, color = method_factor)) +
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
    geom_point(alpha = 0.8, position = position_dodge(width = 0.6)) +
    geom_errorbar(
      aes(ymin = q25, ymax = q75),
      width = 0.2, linewidth = 0.8, alpha = 0.8,
      position = position_dodge(width = 0.6)
    ) +
    facet_wrap(~facet_label, nrow = 1, labeller = label_parsed) +
    scale_x_discrete(labels = paste0(lev, "%")) +
    scale_color_manual(
      values = pop_colors,
      name = "Base Population",
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
    labs(x = "Downsample Rate (%)", y = "Bridging Metric", title = panel_title)

  p
}

upper_levels <- c(100, 50, 25, 10)
lower_levels <- c(5, 4, 3, 2, 1)

left_panel <- build_panel(upper_levels, "100–10%")
right_panel <- build_panel(lower_levels, "5–1%")

leg <- cowplot::get_legend(left_panel)

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
# 6. PREPARE DUNN DATA
# ============================================================================

message("\n🔧 Preparing Dunn data...")

pop_order <- c("detected_only", "symptomatic_only", "symptomatic_men", 
               "high_activity_symptomatic_men", "random_subsample_high_activity")

standardize_comparison <- function(comp_string) {
  pops <- str_split(comp_string, " - ")[[1]]
  pop1 <- trimws(pops[1])
  pop2 <- trimws(pops[2])
  
  pos1 <- match(pop1, pop_order)
  pos2 <- match(pop2, pop_order)
  
  if (pos1 < pos2) {
    return(paste(pop1, pop2, sep = " - "))
  } else {
    return(paste(pop2, pop1, sep = " - "))
  }
}

dunn_data <- dunn_full %>%
  filter(!grepl("^all", `Population Comparison`, ignore.case = TRUE)) %>%
  mutate(
    comparison_standardized = map_chr(`Population Comparison`, standardize_comparison),
    comparison_order = case_when(
      comparison_standardized == "detected_only - symptomatic_only" ~ 1,
      comparison_standardized == "detected_only - symptomatic_men" ~ 2,
      comparison_standardized == "detected_only - high_activity_symptomatic_men" ~ 3,
      comparison_standardized == "detected_only - random_subsample_high_activity" ~ 4,
      comparison_standardized == "symptomatic_only - symptomatic_men" ~ 5,
      comparison_standardized == "symptomatic_only - high_activity_symptomatic_men" ~ 6,
      comparison_standardized == "symptomatic_only - random_subsample_high_activity" ~ 7,
      comparison_standardized == "symptomatic_men - high_activity_symptomatic_men" ~ 8,
      comparison_standardized == "symptomatic_men - random_subsample_high_activity" ~ 9,
      comparison_standardized == "high_activity_symptomatic_men - random_subsample_high_activity" ~ 10,
      TRUE ~ 99
    ),
    comparison_pretty = comparison_standardized %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    comparison_factor = fct_reorder(comparison_pretty, comparison_order)
  )

message("✅ Dunn data prepared")

# ============================================================================
# 7. BUILD DUNN HEATMAP
# ============================================================================

message("\n🎨 Creating Dunn heatmap...")

heatmap <- ggplot(dunn_data, 
                  aes(x = factor(`Downsample Rate (%)`, levels = c(100, 50, 25, 10, 5, 4, 3, 2, 1)),
                      y = comparison_factor,
                      fill = `Significant (FDR 0.05)`)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(
    values = c("TRUE" = "#8B008B", "FALSE" = "#C0C0E0"),
    labels = c("TRUE" = "Significant (FDR < 0.05)", "FALSE" = "Not Significant"),
    name = "Result"
  ) +
  scale_y_discrete(limits = rev(levels(dunn_data$comparison_factor))) +
  facet_wrap(~factor(paste0("ρ = ", Rho), 
                     levels = paste0("ρ = ", sort(unique(Rho)))), 
             nrow = 1) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8.5),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Significant Pairwise Comparisons Across Downsampling Rates",
    x = "Downsample Rate (%)",
    y = "Population Comparison"
  )

message("✅ Dunn heatmap created")

# ============================================================================
# 8. EXTRACT RHO=0.5 PANELS FOR COMBINED FIGURE
# ============================================================================

message("\n🎨 Extracting rho=0.5 panels...")

# Population name to color mapping
pop_name_to_color <- c(
  "Detected Only" = "#52B88F",
  "Symptomatic Only" = "#3FA27A",
  "Symptomatic Men" = "#2C7A78",
  "High Activity Symptomatic Men" = "#164B53",
  "Random Subsample High Activity" = "#2c7bb6"
)

# Extract panel A (downsampling, rho=0.5) - WITH LEGEND AND EXTENDED LINE
plot_A_data <- downsamples_plot_data %>%
  filter(p_msmw_w == 0.5)

gt_A_data <- gt_data_expanded %>%
  filter(p_msmw_w == 0.5)

plot_panel_A <- ggplot(
  plot_A_data,
  aes(x = downsample_factor, y = median_value, color = method_factor)
) +
  geom_ribbon(
    data = gt_A_data,
    aes(x = downsample_factor, ymin = q25, ymax = q75, group = p_msmw_w),
    alpha = 0.35, fill = "#D73027", color = NA, inherit.aes = FALSE
  ) +
  geom_line(
    data = gt_A_data,
    aes(x = downsample_factor, y = median_value, group = p_msmw_w),
    color = "#D73027", linewidth = 1.5, inherit.aes = FALSE
  ) +
  geom_point(alpha = 0.8, position = position_dodge(width = 0.6)) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.2, linewidth = 0.8, alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  scale_x_discrete(labels = paste0(downsample_levels, "%")) +
  scale_color_manual(
    values = pop_colors,
    name = "Base",
    limits = method_levels,
    breaks = method_levels
  ) +
  scale_y_continuous(limits = c(-0.002, NA)) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = rep("blank", length(method_levels)),
        shape = rep(16, length(method_levels)),
        alpha = rep(0.8, length(method_levels)),
        size = rep(2, length(method_levels))
      ),
      ncol = 1
    )
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.5, "lines")
  ) +
  labs(
    x = "Downsample Rate (%)",
    y = "Bridging Value (Product)"
  )

message("✅ Panel A extracted")

# Extract panel B (Dunn heatmap, rho=0.5) - WITH SPLIT-COLORED Y-AXIS
dunn_B_data <- dunn_data %>%
  filter(Rho == 0.5)

# Function to create split-colored comparison labels
create_colored_comparison <- function(comparison_text) {
  # Split on " - "
  parts <- str_split(comparison_text, " - ", n = 2)[[1]]
  
  if (length(parts) != 2) return(comparison_text)
  
  pop1 <- trimws(parts[1])
  pop2 <- trimws(parts[2])
  
  # Get colors
  color1 <- pop_name_to_color[pop1]
  color2 <- pop_name_to_color[pop2]
  
  if (is.na(color1)) color1 <- "black"
  if (is.na(color2)) color2 <- "black"
  
  # Create markdown with separate colors
  colored_text <- glue::glue(
    "<span style='color:{color1}'>{pop1}</span> - <span style='color:{color2}'>{pop2}</span>"
  )
  
  return(colored_text)
}

# Create colored labels
dunn_B_data <- dunn_B_data %>%
  mutate(
    comparison_colored = sapply(as.character(comparison_pretty), create_colored_comparison)
  )

# Create a named vector for scale_y_discrete
colored_labels <- setNames(
  dunn_B_data$comparison_colored,
  dunn_B_data$comparison_factor
)

plot_panel_B <- ggplot(dunn_B_data, 
                  aes(x = factor(`Downsample Rate (%)`, levels = c(100, 50, 25, 10, 5, 4, 3, 2, 1)),
                      y = comparison_factor,
                      fill = `Significant (FDR 0.05)`)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(
    values = c("TRUE" = "#8B008B", "FALSE" = "#C0C0E0"),
    labels = c("TRUE" = "Significant (FDR < 0.05)", "FALSE" = "Not Significant"),
    name = "Result"
  ) +
  scale_y_discrete(
    limits = rev(levels(dunn_B_data$comparison_factor)),
    labels = colored_labels
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = ggtext::element_markdown(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(1, "lines")
  ) +
  labs(
    x = "Downsample Rate (%)",
    y = "Population Comparison"
  )

message("✅ Panel B extracted with split-colored y-axis labels")

# Extract legends separately
leg_A <- cowplot::get_legend(plot_panel_A)
leg_B <- cowplot::get_legend(plot_panel_B)

# Stack legends vertically with MINIMAL spacing (no white space between them)
combined_legends <- plot_grid(
  leg_A,
  leg_B,
  nrow = 2,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "v",
  axis = "lr",
  margin = c(0, 0, 0, 0)  # Remove all margins between legends
)

# Create plots without legends
plot_A_no_leg <- plot_panel_A + theme(legend.position = "none")
plot_B_no_leg <- plot_panel_B + theme(legend.position = "none")

# MODIFIED: Panel A smaller relative to Panel B (0.7:1.3)
# This gives B more room for its long y-axis labels
plots_combined <- plot_grid(
  plot_A_no_leg,
  plot_B_no_leg,
  nrow = 1,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 14,
  label_fontface = "bold",
  rel_widths = c(0.7, 1.3)  # A is 70%, B is 130%
)

# Add stacked legends to the right with NO gap
combined_AB <- plot_grid(
  plots_combined,
  combined_legends,
  nrow = 1,
  ncol = 2,
  rel_widths = c(1, 0.35),  # ← INCREASED from 0.2 to 0.35
  axis = "tb",
  align = "v"
)

message("✅ Combined A+B created with stacked (no-gap) legends and resized panels")
# ============================================================================
# 9. SAVE PLOTS
# ============================================================================

message("\n💾 Saving all plots...")

# Unified plot
ggsave(file.path(FIGURES_DIR, "01_downsamples_unified.tiff"),
       plot = plot_unified, width = 11, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "01_downsamples_unified.png"),
       plot = plot_unified, width = 11, height = 5, dpi = 300, bg = "white")
message("   ✅ 01_downsamples_unified.tiff/.png")

# Side-by-side plot
ggsave(file.path(FIGURES_DIR, "02_downsamples_sidebyside.tiff"),
       plot = plot_sidebyside, width = 14, height = 6, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "02_downsamples_sidebyside.png"),
       plot = plot_sidebyside, width = 14, height = 6, dpi = 300, bg = "white")
message("   ✅ 02_downsamples_sidebyside.tiff/.png")

# Dunn heatmap
ggsave(file.path(FIGURES_DIR, "03_dunn_heatmap_faceted.tiff"),
       plot = heatmap, width = 13, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "03_dunn_heatmap_faceted.png"),
       plot = heatmap, width = 13, height = 5, dpi = 300, bg = "white")
message("   ✅ 03_dunn_heatmap_faceted.tiff/.png")

# Combined A+B (rho=0.5 panels)
ggsave(file.path(FIGURES_DIR, "04_combined_rho05_panels.tiff"),
       plot = combined_AB, width = 12, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "04_combined_rho05_panels.png"),
       plot = combined_AB, width = 12, height = 5, dpi = 300, bg = "white")
message("   ✅ 04_combined_rho05_panels.tiff/.png")

message("\n✅ All plots saved")

# ============================================================================
# 10. SUMMARY
# ============================================================================

message("\n", strrep("=", 80))
message("✨ DOWNSAMPLING & DUNN ANALYSIS COMPLETE!")
message(strrep("=", 80))

message("\n📊 PLOTS CREATED:")
message("   • 01_downsamples_unified (all downsampling rates)")
message("   • 02_downsamples_sidebyside (100-10% left, 5-1% right)")
message("   • 03_dunn_heatmap_faceted (all rho values)")
message("   • 04_combined_rho05_panels (A: downsampling + B: Dunn for rho=0.5)")

message("\n📁 OUTPUT LOCATIONS:")
message(strrep("-", 80))
message("Figures: ", FIGURES_DIR)
message("  ✅ 01_downsamples_unified.tiff/.png")
message("  ✅ 02_downsamples_sidebyside.tiff/.png")
message("  ✅ 03_dunn_heatmap_faceted.tiff/.png")
message("  ✅ 04_combined_rho05_panels.tiff/.png")

message(strrep("=", 80))
message("✅ Script complete!")

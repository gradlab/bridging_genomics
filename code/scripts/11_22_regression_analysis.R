# ============================================================================
# 11.22 — REGRESSION ANALYSIS: GT vs PHYLO (BASE POPULATIONS)
# ============================================================================
# Purpose: Linear regression comparing phylo estimates to ground truth
# - Pool across rho values
# - Fit regression: y = GT, x = Phylo estimate
# - Extract slope, intercept, R², correlation, and Wald test p-values
# - Create faceted scatter plot with CIs and regression lines
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(cowplot)
})

# ============================================================================
# 1. SETUP DIRECTORIES
# ============================================================================

FIGURES_DIR <- "../11_22_output/figures/regression_analysis"
TABLES_DIR <- "../11_22_output/tables/regression_tables"

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

# Phylo data for full sequences (complete dataset at 100%)
phylo_fullseq <- phylo_full %>%
  filter(measure == "total_between_network_proportion",
         dataset_id == "full_sequences") %>%
  mutate(
    MSM.MSMW = phylo_30_MSM.MSMW,
    MSMW.MSW = phylo_30_MSMW.MSW,
    product = MSM.MSMW * MSMW.MSW
  ) %>%
  select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, product)

# Phylo data for all complete populations (100% downsample)
phylo_complete_datasets <- phylo_full %>%
  filter(measure == "total_between_network_proportion", downsample_rate == 100) %>%
  mutate(
    MSM.MSMW = phylo_30_MSM.MSMW,
    MSMW.MSW = phylo_30_MSMW.MSW,
    product = MSM.MSMW * MSMW.MSW,
    population = dataset_id
  ) %>%
  select(sim_id, p_msmw_w, population, MSM.MSMW, MSMW.MSW, product)

message("✅ Datasets prepared")

# ============================================================================
# 4. PREPARE REGRESSION DATA
# ============================================================================

message("\n📊 Preparing regression data...")

# Pretty labels for populations
pop_labels <- c(
  "full_sequences" = "All Infections",
  "detected_only_complete" = "Detected Only",
  "symptomatic_only_complete" = "Symptomatic Only",
  "symptomatic_men_complete" = "Symptomatic Men",
  "high_activity_symptomatic_men_complete" = "High Activity Symptomatic Men",
  "random_subsample_high_activity_complete" = "Random Subsample High Activity"
)

# 1. Full sequences regression data
regression_data_full <- gt_phylosims %>%
  select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, product) %>%
  rename_all(~paste0("gt_", .)) %>%
  rename(sim_id = gt_sim_id, p_msmw_w = gt_p_msmw_w) %>%
  left_join(
    phylo_fullseq %>%
      select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, product) %>%
      rename_all(~paste0("phylo_", .)) %>%
      rename(sim_id = phylo_sim_id, p_msmw_w = phylo_p_msmw_w),
    by = c("sim_id", "p_msmw_w")
  ) %>%
  mutate(population = "full_sequences")

# 2. Other complete populations
regression_data_complete <- phylo_complete_datasets %>%
  select(sim_id, p_msmw_w, population, MSM.MSMW, MSMW.MSW, product) %>%
  rename_all(~paste0("phylo_", .)) %>%
  rename(
    sim_id = phylo_sim_id, 
    p_msmw_w = phylo_p_msmw_w, 
    population = phylo_population
  ) %>%
  left_join(
    gt_phylosims %>%
      select(sim_id, p_msmw_w, MSM.MSMW, MSMW.MSW, product) %>%
      rename_all(~paste0("gt_", .)) %>%
      rename(sim_id = gt_sim_id, p_msmw_w = gt_p_msmw_w),
    by = c("sim_id", "p_msmw_w")
  )

# Combine all
regression_data <- bind_rows(regression_data_full, regression_data_complete) %>%
  mutate(
    pop_label = recode(population, !!!pop_labels),
    pop_factor = factor(pop_label, levels = unname(pop_labels))
  )

message("✅ Regression data prepared: ", nrow(regression_data), " rows")

# ============================================================================
# 5. FIT LINEAR REGRESSIONS (PRODUCT MEASURE)
# ============================================================================

message("\n📈 Fitting linear regressions...")

# Focus on product measure, pool across all rho values
regression_results_base <- regression_data %>%
  select(sim_id, population, pop_label, pop_factor, gt_product, phylo_product) %>%
  rename(gt = gt_product, phylo = phylo_product) %>%
  drop_na() %>%
  group_by(population, pop_label, pop_factor) %>%
  nest() %>%
  mutate(
    # Fit linear regression
    model = map(data, ~lm(gt ~ phylo, data = .x)),
    # Extract model statistics
    glance_out = map(model, broom::glance),
    tidy_out = map(model, broom::tidy),
    # Extract correlation
    correlation = map_dbl(data, ~cor(.x$phylo, .x$gt))
  ) %>%
  ungroup()

# Extract coefficients separately
coefficients_df <- map_df(regression_results_base$tidy_out, ~{
  tidy_df <- .x
  tibble(
    intercept = tidy_df$estimate[1],
    slope = tidy_df$estimate[2],
    slope_se = tidy_df$std.error[2],
    slope_p = tidy_df$p.value[2]  # Wald test p-value
  )
})

# Combine everything
regression_results <- regression_results_base %>%
  unnest(glance_out) %>%
  select(population, pop_label, pop_factor, r.squared, adj.r.squared, correlation) %>%
  bind_cols(coefficients_df) %>%
  ungroup()

message("✅ Regressions fitted: ", nrow(regression_results), " models")

# ============================================================================
# 6. PREPARE PLOTTING DATA
# ============================================================================

message("\n🎨 Preparing plotting data...")

# Scatter plot data (product measure, pooled across rho)
scatter_data <- regression_data %>%
  select(sim_id, population, pop_label, pop_factor, 
         gt_product, phylo_product) %>%
  rename(gt = gt_product, phylo = phylo_product) %>%
  drop_na()

# Build regression lines (fit only, no CI)
regression_lines <- map_df(1:nrow(regression_results_base), ~{
  row <- regression_results_base[.x, ]
  model <- row$model[[1]]
  x_range <- seq(min(scatter_data$phylo), max(scatter_data$phylo), 
                 length.out = 100)
  pred <- predict(model, newdata = data.frame(phylo = x_range))
  
  tibble(
    phylo = x_range,
    fit = pred,
    population = row$population,
    pop_label = row$pop_label,
    pop_factor = row$pop_factor
  )
})

# Annotations for facets (R² and slope only)
facet_annotations <- regression_results %>%
  select(pop_factor, r.squared, slope) %>%
  mutate(
    label = paste0("R² = ", round(r.squared, 3), "\n", 
                   "slope = ", round(slope, 3))
  )

message("✅ Plotting data prepared")

# ============================================================================
# 7. CREATE FOUR SEPARATE PLOTS
# ============================================================================

message("\n🎨 Creating regression plots...")

# ============================================================================
# PLOT 1: All Infections Only (for combined figure with aligned axes)
# ============================================================================

scatter_data_allinfections <- scatter_data %>%
  filter(pop_label == "All Infections")

regression_lines_allinfections <- regression_lines %>%
  filter(pop_label == "All Infections")

regression_results_allinfections <- regression_results %>%
  filter(pop_label == "All Infections")

facet_annotations_allinfections <- facet_annotations %>%
  filter(pop_factor == "All Infections")

# Get axis limits from panel A (Bridging Metric/product panel)
plot_df_3panel <- read.csv("../11_22_output/tables/gt_vs_phylo_allseq/allseq_comparison_plot_data.csv")
product_data <- plot_df_3panel %>% filter(name == "product")
y_min_A <- min(product_data$q25, na.rm = TRUE)
y_max_A <- max(product_data$q75, na.rm = TRUE)
x_min_A <- min(scatter_data_allinfections$phylo, scatter_data_allinfections$gt, na.rm = TRUE)
x_max_A <- max(scatter_data_allinfections$phylo, scatter_data_allinfections$gt, na.rm = TRUE)

# Add some padding
y_range <- y_max_A - y_min_A
x_range <- x_max_A - x_min_A
y_min_A <- y_min_A - 0.05 * y_range
y_max_A <- y_max_A + 0.05 * y_range
x_min_A <- x_min_A - 0.05 * x_range
x_max_A <- x_max_A + 0.05 * x_range

plot_allinfections <- ggplot(scatter_data_allinfections, aes(x = phylo, y = gt)) +
  geom_point(aes(color = pop_factor), alpha = 0.6, size = 2.5) +
  geom_line(
    data = regression_lines_allinfections,
    aes(x = phylo, y = fit, color = pop_factor),
    linewidth = 1.2, alpha = 0.9, inherit.aes = FALSE
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#D73027", linewidth = 0.8, alpha = 0.4) +
  geom_text(
    data = facet_annotations_allinfections,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.2,
    size = 3.5, inherit.aes = FALSE,
    family = "monospace", color = "black"
  ) +
  scale_color_manual(
    values = c("All Infections" = "#4DDADA"),
    name = "Base Population"
  ) +
  scale_x_continuous(limits = c(x_min_A, x_max_A)) +
  scale_y_continuous(limits = c(y_min_A, y_max_A)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  ) +
  labs(
    x = "Bridging Metric - Phylodynamic Estimate",
    y = "Bridging Metric - Ground Truth"
  )

message("✅ All Infections plot created (axes aligned to panel A)")
# ============================================================================
# PLOT 2: 5 Base Populations (1 row, 5 columns) - NO TITLE
# ============================================================================

scatter_data_5pop <- scatter_data %>%
  filter(pop_label != "All Infections") %>%
  mutate(
    pop_factor = factor(pop_label, 
                        levels = c("Detected Only", "Symptomatic Only", 
                                  "Symptomatic Men", "High Activity Symptomatic Men",
                                  "Random Subsample High Activity"))
  )

regression_lines_5pop <- regression_lines %>%
  filter(pop_label != "All Infections") %>%
  mutate(
    pop_factor = factor(pop_label, 
                        levels = c("Detected Only", "Symptomatic Only", 
                                  "Symptomatic Men", "High Activity Symptomatic Men",
                                  "Random Subsample High Activity"))
  )

facet_annotations_5pop <- facet_annotations %>%
  filter(pop_factor != "All Infections") %>%
  mutate(
    pop_factor = factor(pop_factor, 
                        levels = c("Detected Only", "Symptomatic Only", 
                                  "Symptomatic Men", "High Activity Symptomatic Men",
                                  "Random Subsample High Activity"))
  )

plot_5populations_horizontal <- ggplot(scatter_data_5pop, aes(x = phylo, y = gt)) +
  geom_point(aes(color = pop_factor), alpha = 0.6, size = 2) +
  geom_line(
    data = regression_lines_5pop,
    aes(x = phylo, y = fit, color = pop_factor),
    linewidth = 1, alpha = 0.9, inherit.aes = FALSE
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#D73027", linewidth = 0.7, alpha = 0.4) +
  geom_text(
    data = facet_annotations_5pop,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.2,
    size = 3, inherit.aes = FALSE,
    family = "monospace", color = "black"
  ) +
  facet_wrap(~pop_factor, scales = "free", ncol = 5, nrow = 1) +
  scale_color_manual(
    values = c(
      "Detected Only" = "#52B88F",
      "Symptomatic Only" = "#3FA27A",
      "Symptomatic Men" = "#2C7A78",
      "High Activity Symptomatic Men" = "#164B53",
      "Random Subsample High Activity" = "#2c7bb6"
    ),
    name = "Base Population"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Bridging Metric - Phylodynamic Estimate",
    y = "Bridging Metric - Ground Truth"
  )

message("✅ 5-population horizontal plot created (1 row × 5 columns, no title)")

# ============================================================================
# PLOT 3: 5 Base Populations (5 rows, 1 column) - VERTICAL
# ============================================================================

plot_5populations_vertical <- ggplot(scatter_data_5pop, aes(x = phylo, y = gt)) +
  geom_point(aes(color = pop_factor), alpha = 0.6, size = 2) +
  geom_line(
    data = regression_lines_5pop,
    aes(x = phylo, y = fit, color = pop_factor),
    linewidth = 1, alpha = 0.9, inherit.aes = FALSE
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#D73027", linewidth = 0.7, alpha = 0.4) +
  geom_text(
    data = facet_annotations_5pop,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.2,
    size = 3, inherit.aes = FALSE,
    family = "monospace", color = "black"
  ) +
  facet_wrap(~pop_factor, scales = "free", ncol = 1, nrow = 5) +
  scale_color_manual(
    values = c(
      "Detected Only" = "#52B88F",
      "Symptomatic Only" = "#3FA27A",
      "Symptomatic Men" = "#2C7A78",
      "High Activity Symptomatic Men" = "#164B53",
      "Random Subsample High Activity" = "#2c7bb6"
    ),
    name = "Base Population"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Bridging Metric - Phylodynamic Estimate",
    y = "Bridging Metric - Ground Truth"
  )

message("✅ 5-population vertical plot created (5 rows × 1 column)")
# ============================================================================
# PLOT 3B: 5 Base Populations (2 rows: 3+2 centered) - 2x3 GRID
# ============================================================================

suppressPackageStartupMessages({
  library(ggtext)
})

# Create individual plots for each population
pop_order_grid <- c("Detected Only", "Symptomatic Only", "Symptomatic Men", 
                    "High Activity Symptomatic Men", "Random Subsample High Activity")

individual_plots <- map(pop_order_grid, function(pop_name) {
  
  scatter_single <- scatter_data_5pop %>%
    filter(pop_label == pop_name)
  
  regression_single <- regression_lines_5pop %>%
    filter(pop_label == pop_name)
  
  annotation_single <- facet_annotations_5pop %>%
    filter(pop_factor == pop_name)
  
  pop_color <- case_when(
    pop_name == "Detected Only" ~ "#52B88F",
    pop_name == "Symptomatic Only" ~ "#3FA27A",
    pop_name == "Symptomatic Men" ~ "#2C7A78",
    pop_name == "High Activity Symptomatic Men" ~ "#164B53",
    pop_name == "Random Subsample High Activity" ~ "#2c7bb6"
  )
  
  ggplot(scatter_single, aes(x = phylo, y = gt)) +
    geom_point(alpha = 0.6, size = 2.5, color = pop_color) +
    geom_line(
      data = regression_single,
      aes(x = phylo, y = fit),
      linewidth = 1.1, alpha = 0.9, color = pop_color
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "#D73027", linewidth = 0.7, alpha = 0.4) +
    geom_text(
      data = annotation_single,
      aes(label = label),
      x = -Inf, y = Inf,
      hjust = -0.05, vjust = 1.2,
      size = 3.2, inherit.aes = FALSE,
      family = "monospace", color = "black"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 9),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5,
        color = "black"
      )
    ) +
    labs(
      x = "Bridging Metric - Phylodynamic Estimate",
      y = "Bridging Metric - Ground Truth",
      title = pop_name
    )
})

# Top row: 3 plots
top_row <- plot_grid(
  individual_plots[[1]],  # Detected Only
  individual_plots[[2]],  # Symptomatic Only
  individual_plots[[3]],  # Symptomatic Men
  nrow = 1,
  ncol = 3,
  rel_widths = c(1, 1, 1)
)

# Bottom row: 2 plots centered (using NULL spacers)
bottom_row <- plot_grid(
  NULL,
  individual_plots[[4]],  # High Activity Symptomatic Men
  individual_plots[[5]],  # Random Subsample High Activity
  NULL,
  nrow = 1,
  ncol = 4,
  rel_widths = c(0.5, 1, 1, 0.5)
)

# Combine top and bottom rows
plot_5populations_grid <- plot_grid(
  top_row,
  bottom_row,
  nrow = 2,
  ncol = 1,
  rel_heights = c(1, 1)
)

message("✅ 5-population grid plot created (2 rows: 3 top + 2 centered bottom with grey title boxes)")
# ============================================================================
# PLOT 4: Original Full Scatter Plot (all 6 populations)
# ============================================================================

scatter_plot <- ggplot(scatter_data, aes(x = phylo, y = gt)) +
  # Points
  geom_point(aes(color = pop_factor), alpha = 0.6, size = 2) +
  # Regression lines
  geom_line(
    data = regression_lines,
    aes(x = phylo, y = fit, color = pop_factor),
    linewidth = 1, alpha = 0.9, inherit.aes = FALSE
  ) +
  # 1:1 reference line (perfect agreement)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#D73027", linewidth = 0.8, alpha = 0.4) +
  # Annotations in upper left corner
  geom_text(
    data = facet_annotations,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.2,
    size = 3.2, inherit.aes = FALSE,
    family = "monospace", color = "black"
  ) +
  # Facets by population
  facet_wrap(~pop_factor, scales = "free", ncol = 2) +
  # Color scale
  scale_color_manual(
    values = c(
      "All Infections" = "#4DDADA",
      "Detected Only" = "#52B88F",
      "Symptomatic Only" = "#3FA27A",
      "Symptomatic Men" = "#2C7A78",
      "High Activity Symptomatic Men" = "#164B53",
      "Random Subsample High Activity" = "#2c7bb6"
    ),
    name = "Base Population"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    x = "Bridging Metric - Phylodynamic Estimate",
    y = "Bridging Metric - Ground Truth",
    title = "Correlation of Ground Truth and Phylodynamic Estimates\nby Base Population Bias"
  )

message("✅ Full scatter plot created")

# ============================================================================
# PLOT 5: Combined Figure (A + B with 3-panel from first script)
# ============================================================================

message("\n🎨 Creating combined figure with GT styling...")

# Load the 3-panel plot from first script
plot_3panel <- readRDS("../11_22_output/figures/gt_vs_phylo_allseq/plot_3panel_allseq.rds")

# Remove title from 3-panel plot
plot_A_clean <- plot_3panel + 
  labs(title = NULL) +
  theme(plot.title = element_blank())

# Remove title from plot B
plot_B_clean <- plot_allinfections + 
  labs(title = NULL) +
  theme(plot.title = element_blank())

# Align plots on bottom and left axes
aligned_plots <- align_plots(plot_A_clean, plot_B_clean, align = "hv", axis = "b")

# Combine aligned plots using plot_grid
p_combined <- plot_grid(
  aligned_plots[[1]], 
  aligned_plots[[2]],
  ncol = 2, 
  rel_widths = c(1.8, 1.2),
  labels = c("A", "B"),
  label_size = 14,
  label_fontface = "bold"
)

# Add overall title
combined_figure <- plot_grid(
  ggdraw() + draw_label("Ground Truth vs Phylodynamic Estimates: All Sequences", 
                        fontface = "bold", size = 14, hjust = 0.5),
  p_combined,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

message("✅ Combined figure created with aligned x-axes (A + B)")

# ============================================================================
# 8. SAVE PLOTS
# ============================================================================

message("\n💾 Saving all plots...")

# ============================================================================
# Combined figure (A + B) - TIFF & PNG
# ============================================================================

ggsave(file.path(FIGURES_DIR, "01_combined_allinfections_regression.tiff"),
       plot = combined_figure, width = 11, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "01_combined_allinfections_regression.png"),
       plot = combined_figure, width = 11, height = 5, dpi = 300, bg = "white")
message("   ✅ 01_combined_allinfections_regression.tiff/.png")

# ============================================================================
# 5-population regression plot (horizontal) - TIFF & PNG
# ============================================================================

ggsave(file.path(FIGURES_DIR, "02_regression_5populations_horizontal.tiff"),
       plot = plot_5populations_horizontal, width = 15, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "02_regression_5populations_horizontal.png"),
       plot = plot_5populations_horizontal, width = 15, height = 5, dpi = 300, bg = "white")
message("   ✅ 02_regression_5populations_horizontal.tiff/.png")

# ============================================================================
# 5-population regression plot (vertical) - TIFF & PNG
# ============================================================================

ggsave(file.path(FIGURES_DIR, "03_regression_5populations_vertical.tiff"),
       plot = plot_5populations_vertical, width = 5, height = 14, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "03_regression_5populations_vertical.png"),
       plot = plot_5populations_vertical, width = 5, height = 14, dpi = 300, bg = "white")
message("   ✅ 03_regression_5populations_vertical.tiff/.png")

# ============================================================================
# 5-population regression plot (grid: 2 rows × 3 columns) - TIFF & PNG
# ============================================================================

ggsave(file.path(FIGURES_DIR, "03_regression_5populations_grid.tiff"),
       plot = plot_5populations_grid, width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "03_regression_5populations_grid.png"),
       plot = plot_5populations_grid, width = 14, height = 8, dpi = 300, bg = "white")
message("   ✅ 03_regression_5populations_grid.tiff/.png")

# ============================================================================
# Original full scatter plot (all 6 populations) - TIFF & PNG
# ============================================================================

ggsave(file.path(FIGURES_DIR, "04_regression_scatter_plot_all6populations.tiff"),
       plot = scatter_plot, width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "04_regression_scatter_plot_all6populations.png"),
       plot = scatter_plot, width = 12, height = 10, dpi = 300, bg = "white")
message("   ✅ 04_regression_scatter_plot_all6populations.tiff/.png")

message("\n✅ All plots saved in both TIFF and PNG formats")

# ============================================================================
# 9. SAVE TABLES
# ============================================================================

message("\n📋 Saving tables...")

# Table 1: Regression results summary (with CIs and Wald test p-value)
regression_summary_table <- regression_results %>%
  select(pop_label, slope, slope_se, slope_p, intercept, r.squared, adj.r.squared, correlation) %>%
  mutate(
    slope_ci_lower = round(slope - 1.96*slope_se, 4),
    slope_ci_upper = round(slope + 1.96*slope_se, 4),
    slope = round(slope, 4),
    slope_se = round(slope_se, 6),
    slope_p = format(slope_p, scientific = TRUE, digits = 3),
    intercept = round(intercept, 6),
    r.squared = round(r.squared, 4),
    adj.r.squared = round(adj.r.squared, 4),
    correlation = round(correlation, 4)
  ) %>%
  select(pop_label, slope, slope_se, slope_ci_lower, slope_ci_upper, slope_p, 
         intercept, r.squared, adj.r.squared, correlation) %>%
  rename(
    Population = pop_label,
    Slope = slope,
    Slope_SE = slope_se,
    `Slope 95% CI Lower` = slope_ci_lower,
    `Slope 95% CI Upper` = slope_ci_upper,
    `Slope p-value (Wald)` = slope_p,
    Intercept = intercept,
    R2 = r.squared,
    `Adj. R2` = adj.r.squared,
    `Pearson r` = correlation
  )

write_csv(regression_summary_table,
          file.path(TABLES_DIR, "regression_summary_results.csv"))
message("   ✅ regression_summary_results.csv")

# Table 2: Raw scatter data for reproduction
scatter_data_export <- scatter_data %>%
  select(sim_id, population, pop_label, gt, phylo) %>%
  rename(
    Simulation = sim_id,
    Population_ID = population,
    Population = pop_label,
    Ground_Truth = gt,
    Phylodynamic_Estimate = phylo
  )

write_csv(scatter_data_export,
          file.path(TABLES_DIR, "regression_scatter_data.csv"))
message("   ✅ regression_scatter_data.csv")

# Table 3: Full regression statistics
regression_full_table <- regression_results %>%
  select(population, pop_label, slope, intercept, slope_se, slope_p, r.squared, adj.r.squared, correlation) %>%
  mutate(
    slope = round(slope, 6),
    intercept = round(intercept, 8),
    slope_se = round(slope_se, 6),
    slope_p = format(slope_p, scientific = TRUE, digits = 3),
    r.squared = round(r.squared, 6),
    adj.r.squared = round(adj.r.squared, 6),
    correlation = round(correlation, 6)
  ) %>%
  rename(
    Population_ID = population,
    Population = pop_label,
    Slope = slope,
    Intercept = intercept,
    Slope_SE = slope_se,
    `Slope p-value (Wald)` = slope_p,
    R2 = r.squared,
    Adj_R2 = adj.r.squared,
    Pearson_r = correlation
  )

write_csv(regression_full_table,
          file.path(TABLES_DIR, "regression_full_statistics.csv"))
message("   ✅ regression_full_statistics.csv")

# ============================================================================
# 10. PRINT SUMMARY
# ============================================================================

message("\n", strrep("=", 80))
message("✨ REGRESSION ANALYSIS COMPLETE!")
message(strrep("=", 80))

message("\n📊 REGRESSION RESULTS SUMMARY:")
message(strrep("-", 80))
print(regression_summary_table)

message("\n📈 DETAILED MODEL STATISTICS:")
message(strrep("-", 80))
print(regression_full_table)

message("\n📁 OUTPUT LOCATIONS:")
message(strrep("-", 80))
message("Figures saved to: ", FIGURES_DIR)
message("  ✅ 01_combined_allinfections_regression.tiff/.png")
message("  ✅ 02_regression_5populations_horizontal.tiff/.png (1 row × 5 columns)")
message("  ✅ 03_regression_5populations_vertical.tiff/.png (5 rows × 1 column)")
message("  ✅ 04_regression_scatter_plot_all6populations.tiff/.png (2 rows × 3 columns)")

message("\nTables saved to: ", TABLES_DIR)
message("  ✅ regression_summary_results.csv")
message("  ✅ regression_scatter_data.csv")
message("  ✅ regression_full_statistics.csv")

message(strrep("=", 80))
message("✅ Script complete!")

# ============================================================================
# 11.19 — VISUALIZATION OF DUNN TEST RESULTS - WITH ORDERING
# ============================================================================

library(tidyverse)

FIGURES_DIR <- "../11_19_output/figures/downsampling_dunn_viz"
TABLES_DIR <- "../11_19_output/tables"

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD & FILTER
# ============================================================================

dunn_data <- read_csv(file.path("../11_19_output/tables/downsampling_kw_dunn", "downsampling_dunn_results.csv")) %>%
  filter(!grepl("^all", `Population Comparison`, ignore.case = TRUE))

# ============================================================================
# STANDARDIZE COMPARISON NAMES & CREATE ORDER
# ============================================================================

# Define the populations in hierarchical order
pop_order <- c("detected_only", "symptomatic_only", "symptomatic_men", 
               "high_activity_symptomatic_men", "random_subsample_high_activity")

# Function to standardize comparison names (always put lower hierarchy first)
standardize_comparison <- function(comp_string) {
  # Split by " - "
  pops <- str_split(comp_string, " - ")[[1]]
  pop1 <- trimws(pops[1])
  pop2 <- trimws(pops[2])
  
  # Get their hierarchy positions
  pos1 <- match(pop1, pop_order)
  pos2 <- match(pop2, pop_order)
  
  # Return with lower hierarchy position first
  if (pos1 < pos2) {
    return(paste(pop1, pop2, sep = " - "))
  } else {
    return(paste(pop2, pop1, sep = " - "))
  }
}

# Apply standardization and create ordering
dunn_data <- dunn_data %>%
  mutate(
    comparison_standardized = map_chr(`Population Comparison`, standardize_comparison),
    # Create a numeric order based on the standardized comparison
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
    # Pretty format: replace underscores with spaces and title case
    comparison_pretty = comparison_standardized %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  )

# Create the ordered factor
dunn_data <- dunn_data %>%
  mutate(
    comparison_factor = fct_reorder(comparison_pretty, comparison_order)
  )

message("Rows after filtering: ", nrow(dunn_data))
message("Unique comparisons: ", n_distinct(dunn_data$comparison_pretty))
message("Comparison order:\n", 
        paste(unique(dunn_data %>% 
                       arrange(comparison_order) %>% 
                       pull(comparison_pretty)), 
              collapse = "\n  "))

# ============================================================================
# PLOT
# ============================================================================

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

ggsave(
  file.path(FIGURES_DIR, "01_dunn_heatmap_faceted.tiff"),
  plot = heatmap, width = 13, height = 5, dpi = 300, bg = "white"
)
ggsave(
  file.path(FIGURES_DIR, "01_dunn_heatmap_faceted.png"),
  plot = heatmap, width = 13, height = 5, dpi = 300, bg = "white"
)

message("✅ Done!")
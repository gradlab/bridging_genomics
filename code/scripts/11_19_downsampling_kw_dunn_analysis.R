# ============================================================================
# 11.19 — KRUSKAL-WALLIS + DUNN ANALYSIS: DOWNSAMPLING EFFECTS
# ============================================================================
# Purpose: Compare bridging values across base populations within each
# downsampling rate and rho value using non-parametric tests
# - Kruskal-Wallis test for each rho × downsample_rate
# - Dunn post-hoc test for pairwise comparisons
# - Global BH correction for K-W tests
# - Within-stratum BH correction for Dunn tests
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FSA)
})

# ============================================================================
# 1. SETUP DIRECTORIES
# ============================================================================

FIGURES_DIR <- "../11_19_output/figures/downsampling_kw_dunn"
TABLES_DIR <- "../11_19_output/tables/downsampling_kw_dunn"

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

message("📂 Output directories:")
message("   Figures: ", FIGURES_DIR)
message("   Tables:  ", TABLES_DIR)

# ============================================================================
# 2. LOAD DATA
# ============================================================================

message("\n📂 Loading data...")

phylo_full <- read.csv("../output/11_7_complete_50sims/MASTER_phylogenetic_bridging_30_results.csv")

message("✅ Data loaded")

# ============================================================================
# 3. PREPARE DATA
# ============================================================================

message("\n🔧 Preparing data...")

# Extract between-network bridging data
phylo_btn_ntwk <- phylo_full %>%
  filter(measure == "total_between_network_proportion") %>%
  mutate(
    MSM.MSMW = phylo_30_MSM.MSMW,
    MSMW.MSW = phylo_30_MSMW.MSW,
    product = MSM.MSMW * MSMW.MSW
  ) %>%
  select(sim_id, p_msmw_w, population, dataset_type, downsample_rate, replicate, product)

# For each sim × population × downsample_rate, get median of product across 10 replicates
replicate_medians <- phylo_btn_ntwk %>%
  group_by(sim_id, p_msmw_w, population, downsample_rate) %>%
  summarise(
    median_product = median(product, na.rm = TRUE),
    n_replicates = n_distinct(replicate),
    .groups = "drop"
  )

# Get the 10 values (one per sim) for each population × rho × downsample_rate
comparison_data <- replicate_medians %>%
  group_by(p_msmw_w, downsample_rate, population) %>%
  summarise(
    median_products = list(median_product),
    n_sims = n(),
    .groups = "drop"
  )

message("✅ Data prepared: ", nrow(comparison_data), " population-stratum combinations")

# ============================================================================
# 4. PERFORM KRUSKAL-WALLIS TESTS
# ============================================================================

message("\n📊 Performing Kruskal-Wallis tests...")

# Get unique rho × downsample combinations
kw_strata <- replicate_medians %>%
  distinct(p_msmw_w, downsample_rate)

# Function to perform K-W test
perform_kw <- function(rho, ds_rate, data) {
  # Filter data for this stratum
  stratum_data <- data %>%
    filter(p_msmw_w == rho, downsample_rate == ds_rate)
  
  # Need at least 2 groups with data
  n_groups <- n_distinct(stratum_data$population)
  if (n_groups < 2) {
    return(tibble(
      p_msmw_w = rho,
      downsample_rate = ds_rate,
      n_populations = n_groups,
      n_observations = nrow(stratum_data),
      H_statistic = NA,
      p_value = NA,
      test_status = "insufficient_groups"
    ))
  }
  
  # Perform K-W test
  kw_result <- kruskal.test(median_product ~ population, data = stratum_data)
  
  tibble(
    p_msmw_w = rho,
    downsample_rate = ds_rate,
    n_populations = n_groups,
    n_observations = nrow(stratum_data),
    H_statistic = round(kw_result$statistic, 4),
    p_value = kw_result$p.value,
    test_status = "completed"
  )
}

# Apply K-W test to all strata
kw_results_raw <- map_df(1:nrow(kw_strata), ~{
  row <- kw_strata[.x, ]
  perform_kw(row$p_msmw_w, row$downsample_rate, replicate_medians)
})

# Apply global Benjamini-Hochberg correction across all K-W tests
kw_results <- kw_results_raw %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    sig_BH = p_adj_BH < 0.05
  ) %>%
  mutate(
    p_value = round(p_value, 6),
    p_adj_BH = round(p_adj_BH, 6)
  )

message("✅ K-W tests completed: ", sum(!is.na(kw_results$H_statistic)), " tests")
message("   Significant at FDR 0.05: ", sum(kw_results$sig_BH, na.rm = TRUE), " tests")

# ============================================================================
# 5. PERFORM DUNN POST-HOC TESTS
# ============================================================================

message("\n📊 Performing Dunn post-hoc tests...")

# Function to perform Dunn test
perform_dunn <- function(rho, ds_rate, data) {
  # Filter data for this stratum
  stratum_data <- data %>%
    filter(p_msmw_w == rho, downsample_rate == ds_rate)
  
  # Need at least 2 groups with data
  n_groups <- n_distinct(stratum_data$population)
  if (n_groups < 2) {
    return(NULL)
  }
  
  # Perform Dunn test
  dunn_result <- dunnTest(median_product ~ population, 
                          data = stratum_data,
                          method = "bonferroni") #we do not end up using these 
  
  # Extract results
  dunn_df <- as.data.frame(dunn_result$res) %>%
    as_tibble() %>%
    mutate(
      p_msmw_w = rho,
      downsample_rate = ds_rate,
      .before = everything()
    )
  
  return(dunn_df)
}

# Apply Dunn test to all strata with significant K-W results
dunn_results_raw <- map_df(1:nrow(kw_strata), ~{
  row <- kw_strata[.x, ]
  perform_dunn(row$p_msmw_w, row$downsample_rate, replicate_medians)
})

# Apply within-stratum BH correction for Dunn tests
dunn_results <- dunn_results_raw %>%
  group_by(p_msmw_w, downsample_rate) %>%
  mutate(
    p_adj_BH = p.adjust(P.unadj, method = "BH"),
    sig_BH = p_adj_BH < 0.05
  ) %>%
  ungroup() %>%
  mutate(
    P.unadj = round(P.unadj, 6),
    p_adj_BH = round(p_adj_BH, 6)
  ) %>%
  rename(
    Comparison = Comparison,
    Z_statistic = Z,
    P_unadjusted = P.unadj,
    P_adj_Bonferroni = P.adj,
    P_adj_BH = p_adj_BH,
    Sig_BH = sig_BH
  )

message("✅ Dunn tests completed: ", nrow(dunn_results), " pairwise comparisons")
message("   Significant at FDR 0.05: ", sum(dunn_results$Sig_BH, na.rm = TRUE), " comparisons")

# ============================================================================
# 6. SAVE TABLES
# ============================================================================

message("\n📋 Saving tables...")

# Table 1: Kruskal-Wallis results
kw_summary_table <- kw_results %>%
  select(p_msmw_w, downsample_rate, n_populations, n_observations, 
         H_statistic, p_value, p_adj_BH, sig_BH, test_status) %>%
  rename(
    `Rho` = p_msmw_w,
    `Downsample Rate (%)` = downsample_rate,
    `N Populations` = n_populations,
    `N Observations` = n_observations,
    `H Statistic` = H_statistic,
    `P-value` = p_value,
    `P-adj (BH)` = p_adj_BH,
    `Significant (FDR 0.05)` = sig_BH,
    `Test Status` = test_status
  )

write_csv(kw_summary_table,
          file.path(TABLES_DIR, "downsampling_kw_results.csv"))
message("   ✅ downsampling_kw_results.csv")

# Table 2: Dunn post-hoc results
if (nrow(dunn_results) > 0) {
  dunn_summary_table <- dunn_results %>%
    select(p_msmw_w, downsample_rate, Comparison, Z_statistic, 
           P_unadjusted, P_adj_Bonferroni, P_adj_BH, Sig_BH) %>%
    rename(
      `Rho` = p_msmw_w,
      `Downsample Rate (%)` = downsample_rate,
      `Population Comparison` = Comparison,
      `Z Statistic` = Z_statistic,
      `P (Unadjusted)` = P_unadjusted,
      `P (Bonferroni)` = P_adj_Bonferroni,
      `P (BH)` = P_adj_BH,
      `Significant (FDR 0.05)` = Sig_BH
    )
  
  write_csv(dunn_summary_table,
            file.path(TABLES_DIR, "downsampling_dunn_results.csv"))
  message("   ✅ downsampling_dunn_results.csv")
} else {
  message("   ⚠️  No Dunn results (no K-W tests were significant)")
}

# ============================================================================
# 7. PRINT SUMMARY
# ============================================================================

message("\n", strrep("=", 80))
message("✨ KRUSKAL-WALLIS + DUNN ANALYSIS COMPLETE!")
message(strrep("=", 80))

message("\n📊 KRUSKAL-WALLIS RESULTS SUMMARY:")
message(strrep("-", 80))
print(kw_summary_table)

if (nrow(dunn_results) > 0) {
  message("\n📊 DUNN POST-HOC RESULTS SUMMARY:")
  message(strrep("-", 80))
  print(dunn_summary_table %>% slice(1:20))  # First 20 rows
  if (nrow(dunn_summary_table) > 20) {
    message("   ... and ", nrow(dunn_summary_table) - 20, " more rows")
  }
}

message("\n📁 OUTPUT LOCATIONS:")
message(strrep("-", 80))
message("Tables saved to: ", TABLES_DIR)
message("  ✅ downsampling_kw_results.csv")
if (nrow(dunn_results) > 0) {
  message("  ✅ downsampling_dunn_results.csv")
}

message(strrep("=", 80))
message("✅ Script complete!")
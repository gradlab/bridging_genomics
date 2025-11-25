#!/usr/bin/env Rscript

# Final Tree Filter - Remove Burn-in Transmission Tips
# Takes a tree trimmed by sampling day and removes tips with burn-in transmissions

library(ape)
library(tidyverse)

# ========== FUNCTIONS ==========

parse_tip_name <- function(tip_name) {
  # Parse tip names like 'node_day_epX' 
  clean_name <- gsub("^'|'$", "", tip_name)  # Remove quotes
  parts <- strsplit(clean_name, "_")[[1]]
  
  list(
    infectee_node = as.integer(parts[1]),
    day_of_sampling = as.integer(parts[2])
  )
}

filter_burnin_transmissions <- function(tree_file, transmission_file, output_file, verbose = TRUE) {
  
  if (verbose) cat("🌳 Final Tree Filter - Removing Burn-in Transmission Tips\n")
  
  # Load tree
  if (verbose) cat(sprintf("Loading tree: %s\n", tree_file))
  tree <- read.tree(tree_file)
  if (verbose) cat(sprintf("  Input tips: %d\n", length(tree$tip.label)))
  
  # Parse tip labels
  if (verbose) cat("Parsing tip labels...\n")
  tip_df <- data.frame(
    tip_name = tree$tip.label,
    stringsAsFactors = FALSE
  )
  
  # Extract node and sampling day from each tip
  for (i in 1:nrow(tip_df)) {
    parsed <- parse_tip_name(tip_df$tip_name[i])
    tip_df$infectee_node[i] <- parsed$infectee_node
    tip_df$day_of_sampling[i] <- parsed$day_of_sampling
  }
  
  # Load transmission data
  if (verbose) cat(sprintf("Loading transmission data: %s\n", transmission_file))
  transmission_df <- read.csv(transmission_file)
  
  # Filter to burn-in transmissions only
  burnin_transmissions <- transmission_df %>%
    filter(phase != "post_tx") %>%
    filter(superseded_simultaneous == "False") %>%
    mutate(
      infectee_node = as.integer(infectee_node),
      day_of_sampling = as.integer(day_of_sampling)
    )
  
  if (verbose) cat(sprintf("  Burn-in transmissions: %d\n", nrow(burnin_transmissions)))
  
  # Find tips with burn-in transmissions
  tips_with_burnin <- tip_df %>%
    left_join(burnin_transmissions, by = c("infectee_node", "day_of_sampling")) %>%
    filter(!is.na(phase))
  
  if (verbose) {
    cat(sprintf("Tips with burn-in transmissions: %d\n", nrow(tips_with_burnin)))
    cat("Removing these tips...\n")
  }
  
  # Keep tips WITHOUT burn-in transmissions
  tips_to_keep <- tip_df %>%
    left_join(burnin_transmissions, by = c("infectee_node", "day_of_sampling")) %>%
    filter(is.na(phase)) %>%
    pull(tip_name)
  
  if (verbose) cat(sprintf("Tips to keep: %d\n", length(tips_to_keep)))
  
  # Filter tree
  filtered_tree <- keep.tip(tree, tips_to_keep)
  
  # Save result
  write.tree(filtered_tree, output_file)
  
  if (verbose) {
    cat(sprintf("✅ Final tree saved: %s\n", output_file))
    cat(sprintf("   Final tips: %d\n", length(filtered_tree$tip.label)))
    cat(sprintf("   Removed: %d burn-in transmission tips\n", 
                length(tree$tip.label) - length(filtered_tree$tip.label)))
  }
  
  invisible(filtered_tree)
}

# ========== MAIN ==========

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage: Rscript final_filter.R input_tree transmission_df output_tree\n")
  cat("  Removes tips with burn-in transmissions from pre-trimmed tree\n")
  quit(status = 1)
}

filter_burnin_transmissions(
  tree_file = args[1],
  transmission_file = args[2], 
  output_file = args[3]
)
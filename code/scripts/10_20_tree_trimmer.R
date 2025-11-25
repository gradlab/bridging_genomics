#!/usr/bin/env Rscript

# Tree Trimmer for Seq-Gen
# Trims transmission trees to post-burn-in tips and prepares for Seq-Gen

library(ape)

# ========== HELPER FUNCTIONS ==========

normalize_text <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- enc2utf8(x)
  # strip surrounding single/double quotes
  x <- gsub("^'(.*)'$", "\\1", x)
  x <- gsub('^"(.*)"$', "\\1", x)
  # normalize NBSP to space
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  # normalize any "fancy" dashes to ASCII hyphen
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2015\u2212]", "-", x, perl = TRUE)
  # collapse multiple spaces and trim
  x <- gsub("[ \t]+", " ", x, perl = TRUE)
  trimws(x)
}

parse_episode_label <- function(label) {
  # Parse episode-aware labels like 'node_day_epX' into components
  # First normalize/clean the label to remove quotes
  clean_label <- normalize_text(as.character(label))
  
  parts <- strsplit(clean_label, "_")[[1]]
  if (length(parts) >= 3 && grepl("^ep[0-9]+$", parts[length(parts)])) {
    node_id <- as.numeric(parts[1])
    day_value <- as.numeric(parts[2])
    episode <- as.numeric(sub("^ep", "", parts[length(parts)]))
    return(list(node_id = node_id, day = day_value, episode = episode))
  } else {
    return(list(node_id = NA, day = NA, episode = NA))
  }
}

infer_cutoff_day <- function(sim_dir, transmission_df_path, verbose = TRUE) {
  # Method 1: Try to load parameters file
  params_file <- file.path(sim_dir, "parameters_used.json")
  
  if (file.exists(params_file)) {
    if (verbose) cat("Loading parameters from:", params_file, "\n")
    
    # Read the file
    params_text <- readLines(params_file, warn = FALSE)
    params_json <- paste(params_text, collapse = "\n")
    
    # Simple string extraction - look for the exact patterns
    partnership_line <- grep('"partnership_burnin_days"', params_text, value = TRUE)
    transmission_line <- grep('"transmission_burnin_days"', params_text, value = TRUE)
    
    if (length(partnership_line) > 0 && length(transmission_line) > 0) {
      # Extract numbers using simple regex
      partnership_burnin <- as.numeric(gsub('.*: *([0-9]+).*', '\\1', partnership_line))
      transmission_burnin <- as.numeric(gsub('.*: *([0-9]+).*', '\\1', transmission_line))
      
      if (!is.na(partnership_burnin) && !is.na(transmission_burnin)) {
        cutoff_day <- partnership_burnin + transmission_burnin
        
        if (verbose) {
          cat(sprintf("  Partnership burnin: %d days\n", partnership_burnin))
          cat(sprintf("  Transmission burnin: %d days\n", transmission_burnin))
          cat(sprintf("  Inferred cutoff day: %d\n", cutoff_day))
        }
        
        return(cutoff_day)
      } else {
        stop(sprintf("Could not parse numeric values from parameters file.\n  Partnership line: %s\n  Transmission line: %s", 
                    partnership_line, transmission_line))
      }
    } else {
      stop(sprintf("Could not find required parameters in %s.\n  Looking for: 'partnership_burnin_days' and 'transmission_burnin_days'", 
                  params_file))
    }
  } else {
    stop(sprintf("Parameters file not found: %s", params_file))
  }
}

get_transmission_day <- function(node_id, episode, transmission_df) {
  # Find transmission day for a specific node and episode
  matching_rows <- transmission_df[
    transmission_df$infectee_node == node_id & 
    transmission_df$infectee_episode == episode, 
  ]
  
  if (nrow(matching_rows) > 0) {
    return(matching_rows$day_of_transmission[1])
  } else {
    return(NA)
  }
}

trim_and_clean_for_seqgen <- function(
  in_tree,
  sim_dir = NULL,
  transmission_df_path = NULL,
  cutoff_day = NULL,
  out_tree,
  randomize_resolution = FALSE,
  verbose = TRUE
) {
  tr <- read.tree(in_tree)

  # hard fail on non-finite lengths
  if (!is.null(tr$edge.length)) {
    bad <- which(!is.finite(tr$edge.length))
    if (length(bad)) stop(sprintf("Non-finite branch lengths at edges: %s", paste(bad, collapse = ", ")))
  }

  # Determine cutoff day
  if (is.null(cutoff_day)) {
    if (is.null(sim_dir) && !is.null(transmission_df_path)) {
      sim_dir <- dirname(transmission_df_path)
    }
    
    if (is.null(sim_dir)) {
      stop("Either cutoff_day must be specified, or sim_dir/transmission_df_path must be provided to infer cutoff")
    }
    
    if (is.null(transmission_df_path)) {
      transmission_df_path <- file.path(sim_dir, "transmission_df.csv")
    }
    
    cutoff_day <- infer_cutoff_day(sim_dir, transmission_df_path, verbose)
  }
  
  if (verbose) cat(sprintf("Using cutoff day: %d\n", cutoff_day))

  # Find tips to keep based on cutoff day
  tree_labels <- tr$tip.label
  keep_labels <- c()

    # DEBUG: Show some example tip labels and their parsed values
  if (verbose) {
    cat("Debugging tip labels:\n")
    sample_labels <- head(tree_labels, 10)
    for (label in sample_labels) {
      parsed <- parse_episode_label(label)
      cat(sprintf("  '%s' -> node:%s, day:%s, episode:%s\n", 
                  label, parsed$node_id, parsed$day, parsed$episode))
    }
    cat("\n")
  }
  
  for (label in tree_labels) {
    parsed <- parse_episode_label(label)
    if (!is.na(parsed$day) && parsed$day >= cutoff_day) {
      keep_labels <- c(keep_labels, label)
    }
  }
  
  # DEBUG: Show distribution of days
  if (verbose) {
    all_days <- c()
    for (label in tree_labels) {
      parsed <- parse_episode_label(label)
      if (!is.na(parsed$day)) {
        all_days <- c(all_days, parsed$day)
      }
    }
    if (length(all_days) > 0) {
      cat(sprintf("Day range in tree: %d to %d\n", min(all_days), max(all_days)))
      cat(sprintf("Tips with day >= %d: %d out of %d total\n", 
                  cutoff_day, sum(all_days >= cutoff_day), length(all_days)))
    }
  }

  if (length(keep_labels) == 0) {
    stop(sprintf("No tips found with day >= %d. Check your cutoff day value.", cutoff_day))
  }
  
  if (verbose) {
    cat(sprintf("Found %d post-burn-in tips (day >= %d) out of %d total tips\n", 
                length(keep_labels), cutoff_day, length(tree_labels)))
  }

  # Trim to the selected labels
  tr <- keep.tip(tr, keep_labels)

  # Resolve polytomies, collapse singles, strip internal labels
  tr <- multi2di(tr, random = randomize_resolution)
  tr <- collapse.singles(tr)
  tr$node.label <- NULL

  # Validate strictly bifurcating
  child_tab <- table(tr$edge[,1])
  if (any(child_tab != 2)) {
    bad_nodes <- names(child_tab)[child_tab != 2]
    stop(sprintf("Tree not strictly bifurcating. Offending nodes: %s", paste(bad_nodes, collapse = ", ")))
  }

  # Round-trip serialization
  tmp <- tempfile(fileext = ".newick"); write.tree(tr, file = tmp)
  ser <- readLines(tmp, warn = FALSE)
  if (!length(ser) || !grepl(";$", ser[length(ser)])) stop("Serialized tree does not end with ';'.")

  write.tree(tr, file = out_tree)

  if (verbose) {
    message(sprintf(
      "Trimmed & cleaned for Seq-Gen\n  Kept tips: %d\n  Internal nodes (final): %d\n  Cutoff day: %d\n  Output: %s",
      Ntip(tr), tr$Nnode, cutoff_day, normalizePath(out_tree)
    ))
  }
  invisible(tr)
}

# ========== COMMAND LINE INTERFACE ==========

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0 || any(c("-h", "--help") %in% args)) {
    cat("Tree Trimmer for Seq-Gen\n")
    cat("Trims transmission trees to post-burn-in tips and prepares for Seq-Gen\n\n")
    cat("Usage: Rscript tree_trimmer.R [options] input_tree sim_directory output_tree\n")
    cat("   OR: Rscript tree_trimmer.R [options] input_tree transmission_df output_tree\n\n")
    cat("Arguments:\n")
    cat("  input_tree           Input Newick tree file\n")
    cat("  sim_directory        Simulation output directory (contains parameters_used.json)\n")
    cat("  transmission_df      Transmission dataframe CSV file\n")
    cat("  output_tree          Output trimmed tree file\n\n")
    cat("Options:\n")
    cat("  --cutoff-day DAY     Override cutoff day (otherwise inferred from sim data)\n")
    cat("  --randomize          Randomize polytomy resolution (default: FALSE)\n")
    cat("  --quiet              Suppress verbose output\n")
    cat("  -h, --help           Show this help message\n\n")
    cat("Examples:\n")
    cat("  # Auto-infer cutoff from simulation directory\n")
    cat("  Rscript tree_trimmer.R tree.newick ../output/my_sim/ trimmed.newick\n\n")
    cat("  # Auto-infer cutoff from transmission file\n")
    cat("  Rscript tree_trimmer.R tree.newick transmission_df.csv trimmed.newick\n\n")
    cat("  # Override cutoff day\n")
    cat("  Rscript tree_trimmer.R --cutoff-day 15000 tree.newick sim_dir/ trimmed.newick\n")
    quit(status = 0)
  }
  
  # Parse options
  cutoff_day <- NULL
  randomize_resolution <- FALSE
  verbose <- TRUE
  
  i <- 1
  positional <- c()
  
  while (i <= length(args)) {
    arg <- args[i]
    
    if (arg == "--cutoff-day") {
      if (i == length(args)) stop("--cutoff-day requires a numeric value")
      cutoff_day <- as.numeric(args[i + 1])
      if (is.na(cutoff_day)) stop("--cutoff-day must be a number")
      i <- i + 2
    } else if (arg == "--randomize") {
      randomize_resolution <- TRUE
      i <- i + 1
    } else if (arg == "--quiet") {
      verbose <- FALSE
      i <- i + 1
    } else if (!startsWith(arg, "--")) {
      positional <- c(positional, arg)
      i <- i + 1
    } else {
      stop(sprintf("Unknown option: %s", arg))
    }
  }
  
  # Check required positional arguments
  if (length(positional) != 3) {
    stop("Exactly 3 positional arguments required: input_tree sim_data output_tree")
  }
  
  list(
    input_tree = positional[1],
    sim_data = positional[2],
    output_tree = positional[3],
    cutoff_day = cutoff_day,
    randomize_resolution = randomize_resolution,
    verbose = verbose
  )
}

# ========== MAIN ==========

main <- function() {
  args <- parse_args()
  
  # Validate input files exist
  if (!file.exists(args$input_tree)) {
    stop(sprintf("Input tree file not found: %s", args$input_tree))
  }
  
  # Determine if sim_data is directory or file
  if (dir.exists(args$sim_data)) {
    sim_dir <- normalizePath(args$sim_data, mustWork = FALSE)  # Remove trailing slash
    transmission_df_path <- file.path(sim_dir, "transmission_df.csv")
  } else if (file.exists(args$sim_data)) {
    transmission_df_path <- args$sim_data
    sim_dir <- dirname(args$sim_data)
  } else {
    stop(sprintf("Simulation data not found: %s", args$sim_data))
  }
  
  # Create output directory if needed
  output_dir <- dirname(args$output_tree)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (args$verbose) {
    cat("🌳 Tree Trimming for Seq-Gen\n")
    cat(sprintf("   Input tree: %s\n", args$input_tree))
    cat(sprintf("   Simulation data: %s\n", args$sim_data))
    cat(sprintf("   Output: %s\n", args$output_tree))
    if (!is.null(args$cutoff_day)) {
      cat(sprintf("   Cutoff day (override): %d\n", args$cutoff_day))
    } else {
      cat("   Cutoff day: auto-infer\n")
    }
  }
  
  # Run the trimming
  tryCatch({
    trim_and_clean_for_seqgen(
      in_tree = args$input_tree,
      sim_dir = sim_dir,
      transmission_df_path = transmission_df_path,
      cutoff_day = args$cutoff_day,
      out_tree = args$output_tree,
      randomize_resolution = args$randomize_resolution,
      verbose = args$verbose
    )
    
    if (args$verbose) {
      cat("✅ Tree trimming completed successfully!\n")
    }
    
  }, error = function(e) {
    cat(sprintf("❌ Error: %s\n", e$message))
    quit(status = 1)
  })
}

# Run if called as script
if (!interactive()) {
  main()
}
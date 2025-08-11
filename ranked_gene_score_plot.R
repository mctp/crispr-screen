library(tidyverse)
library(yaml)
library(ggplot2)
library(ggrepel)

config_file <- "config.yaml"

# Function to parse config.yaml
parse_config_yaml <- function(config_file) {
  # Read the YAML file as text and replace unquoted . with quoted "."
  yaml_lines <- readLines(config_file)
  yaml_lines <- gsub('^([[:space:]]*[a-zA-Z0-9_]+:[[:space:]]*)\\.$', '\\1"."', yaml_lines)
  config <- yaml::read_yaml(text = paste(yaml_lines, collapse = "\n"))
  # Replace any "." values with the current working directory
  config <- lapply(config, function(x) if (identical(x, ".")) getwd() else x)
  return(config)
}

# Function to sanitize analysis name for filenames
sanitize_filename <- function(name) {
  # Replace spaces with underscores and remove non-alphanumeric characters except underscore and hyphen
  sanitized <- gsub("[^a-zA-Z0-9_-]", "_", name)
  # Remove multiple consecutive underscores
  sanitized <- gsub("_+", "_", sanitized)
  # Remove leading/trailing underscores
  sanitized <- gsub("^_+|_+$", "", sanitized)
  return(sanitized)
}

# Function to detect data type (RRA vs MLE)
detect_data_type <- function(data) {
  colnames_lower <- tolower(colnames(data))
  
  # Check for RRA-specific columns
  has_rra_cols <- any(c("neg_lfc", "pos_lfc", "neg_fdr", "pos_fdr") %in% colnames_lower)
  
  # Check for MLE-specific columns (beta and fdr with condition prefixes)
  has_mle_cols <- any(grepl("\\|beta$|_beta$", colnames(data))) && 
                  any(grepl("\\|fdr$|_fdr$", colnames(data)))
  
  if (has_rra_cols) {
    return("RRA")
  } else if (has_mle_cols) {
    return("MLE")
  } else {
    return("unknown")
  }
}

# Function to extract conditions from MLE data
extract_mle_conditions <- function(data) {
  # Look for columns ending with |beta or _beta
  beta_cols <- colnames(data)[grepl("\\|beta$|_beta$", colnames(data))]
  
  # Extract condition names by removing the |beta or _beta suffix
  conditions <- gsub("\\|beta$|_beta$", "", beta_cols)
  
  return(conditions)
}

# Function to create ranked gene score plot for RRA data
create_ranked_gene_score_plot_rra <- function(data, comparison_name, result_type = "neg", top_n = 30) {
  
  # Select appropriate columns based on result type
  if (result_type == "neg") {
    lfc_col <- "neg_lfc"
    fdr_col <- "neg_fdr"
    title_suffix <- "Negative Selection"
    y_label <- "Negative Log2 Fold Change"
  } else {
    lfc_col <- "pos_lfc"
    fdr_col <- "pos_fdr"
    title_suffix <- "Positive Selection"
    y_label <- "Positive Log2 Fold Change"
  }
  
  # Clean column names: remove |, convert to lowercase, replace non-alphanumeric with _
  colnames(data) <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(data)))
  
  # Also update lfc_col and fdr_col to lowercase
  lfc_col <- tolower(lfc_col)
  fdr_col <- tolower(fdr_col)
  
  # Find the gene identifier column (commonly 'id' or 'gene')
  gene_id_col <- if ("id" %in% colnames(data)) {
    "id"
  } else if ("gene" %in% colnames(data)) {
    "gene"
  } else {
    colnames(data)[1] # fallback to first column
  }
  
  # Create significance colors
  data <- data %>%
    mutate(
      significance = case_when(
        .data[[fdr_col]] < 0.05 ~ "red",
        .data[[fdr_col]] >= 0.05 & .data[[fdr_col]] < 0.1 ~ "orange",
        TRUE ~ "black"
      )
    ) %>%
    arrange(.data[[lfc_col]])
  
  # Add rank for x-axis
  data <- data %>%
    mutate(rank = row_number())
  
  # Calculate point size based on FDR: sizes are powers of 2 (0.5, 1, 2, 4, 8, 16, ...)
  data <- data %>%
    mutate(
      point_size = case_when(
        .data[[fdr_col]] < 0.0001 ~ 8,
        .data[[fdr_col]] < 0.0005 ~ 6,
        .data[[fdr_col]] < 0.001 ~ 4,
        .data[[fdr_col]] < 0.005 ~ 3,
        .data[[fdr_col]] < 0.01 ~ 2,
        .data[[fdr_col]] < 0.05 ~ 1,
        .data[[fdr_col]] < 0.1 ~ 0.75,
        .data[[fdr_col]] < 0.25 ~ 0.5,
        .data[[fdr_col]] >= 0.25 ~ 0.25,
        TRUE ~ 12
      )
    )
  
  # Identify top genes for labeling (top_n by absolute log fold change)
  top_genes <- data %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    slice_head(n = min(top_n, nrow(data))) %>%
    mutate(label_color = ifelse(.data[[fdr_col]] < 0.01, "black", "gray"))
  
  # Create plot
  p <- ggplot(data, aes(x = rank, y = .data[[lfc_col]])) +
    geom_point(aes(color = significance, size = point_size), alpha = 0.7) +
    scale_color_identity() +
    scale_size_identity() +
    geom_text_repel(
      data = top_genes,
      aes(label = .data[[gene_id_col]]),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    labs(
      title = paste(comparison_name, title_suffix, sep = " - "),
      x = "Gene Rank",
      y = y_label,
      caption = "Red: FDR < 0.05, Orange: 0.05 ≤ FDR < 0.1, Black: FDR ≥ 0.1. Larger points = lower FDR."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  return(p)
}

# Function to create ranked gene score plot for MLE data
create_ranked_gene_score_plot_mle <- function(data, condition_name, analysis_name = "", top_n = 30) {
  
  # Clean column names to standardize separators
  colnames(data) <- gsub("\\|", "_", colnames(data))
  
  # Find beta and fdr columns for this condition
  beta_col <- paste0(condition_name, "_beta")
  fdr_col <- paste0(condition_name, "_fdr")
  
  # Check if columns exist
  if (!(beta_col %in% colnames(data))) {
    stop(paste("Beta column not found for condition:", condition_name, ". Expected:", beta_col))
  }
  if (!(fdr_col %in% colnames(data))) {
    stop(paste("FDR column not found for condition:", condition_name, ". Expected:", fdr_col))
  }
  
  # Find the gene identifier column
  gene_id_col <- if ("Gene" %in% colnames(data)) {
    "Gene"
  } else if ("gene" %in% colnames(data)) {
    "gene"
  } else if ("ID" %in% colnames(data)) {
    "ID"
  } else if ("id" %in% colnames(data)) {
    "id"
  } else {
    colnames(data)[1] # fallback to first column
  }
  
  # Create significance colors
  data <- data %>%
    mutate(
      significance = case_when(
        .data[[fdr_col]] < 0.05 ~ "red",
        .data[[fdr_col]] >= 0.05 & .data[[fdr_col]] < 0.1 ~ "orange",
        TRUE ~ "black"
      )
    ) %>%
    arrange(.data[[beta_col]])
  
  # Add rank for x-axis
  data <- data %>%
    mutate(rank = row_number())
  
  # Calculate point size based on FDR
  data <- data %>%
    mutate(
      point_size = case_when(
        .data[[fdr_col]] < 0.0001 ~ 8,
        .data[[fdr_col]] < 0.0005 ~ 6,
        .data[[fdr_col]] < 0.001 ~ 4,
        .data[[fdr_col]] < 0.005 ~ 3,
        .data[[fdr_col]] < 0.01 ~ 2,
        .data[[fdr_col]] < 0.05 ~ 1,
        .data[[fdr_col]] < 0.1 ~ 0.75,
        .data[[fdr_col]] < 0.25 ~ 0.5,
        .data[[fdr_col]] >= 0.25 ~ 0.25,
        TRUE ~ 12
      )
    )
  
  # Identify top genes for labeling (top_n by absolute beta value)
  top_genes <- data %>%
    arrange(desc(abs(.data[[beta_col]]))) %>%
    slice_head(n = min(top_n, nrow(data))) %>%
    mutate(label_color = ifelse(.data[[fdr_col]] < 0.01, "black", "gray"))
  
  # Create plot title
  plot_title <- if (analysis_name != "") {
    paste(analysis_name, "-", condition_name, "Beta Values")
  } else {
    paste(condition_name, "Beta Values")
  }
  
  # Create plot
  p <- ggplot(data, aes(x = rank, y = .data[[beta_col]])) +
    geom_point(aes(color = significance, size = point_size), alpha = 0.7) +
    scale_color_identity() +
    scale_size_identity() +
    geom_text_repel(
      data = top_genes,
      aes(label = .data[[gene_id_col]]),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    labs(
      title = plot_title,
      x = "Gene Rank",
      y = "Beta Value",
      caption = "Red: FDR < 0.05, Orange: 0.05 ≤ FDR < 0.1, Black: FDR ≥ 0.1. Larger points = lower FDR."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  return(p)
}

if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  # Collect command line arguments as list of paths
  args <- commandArgs(trailingOnly = TRUE)
  gene_files <- args

  # Check existence of each file
  missing_files <- gene_files[!file.exists(gene_files)]
  if (length(missing_files) == length(gene_files)) {
    stop(paste("None of the provided files exist:", paste(missing_files, collapse = ", ")))
  } else if (length(missing_files) > 0) {
    warning(paste("The following files do not exist and will be skipped:", paste(missing_files, collapse = ", ")))
    gene_files <- gene_files[file.exists(gene_files)]
  }

  # Process each file
  for (gene_file in gene_files) {
    cat("Processing:", gene_file, "\n")
    gene_data <- read.delim(gene_file, sep = "\t", header = TRUE)
    
    # Detect data type
    data_type <- detect_data_type(gene_data)
    cat("Detected data type:", data_type, "\n")
    
    # Extract comparison/analysis name from file path
    file_base <- basename(gene_file)
    comparison_name <- sub("\\.gene_summary\\.txt$", "", file_base)
    output_dir <- dirname(gene_file)
    
    if (data_type == "RRA") {
      # Create RRA plots
      neg_plot <- create_ranked_gene_score_plot_rra(gene_data, comparison_name, "neg")
      pos_plot <- create_ranked_gene_score_plot_rra(gene_data, comparison_name, "pos")
      
      # Save RRA plots
      ggsave(
        filename = file.path(output_dir, paste(comparison_name, "negative_selection.png", sep = "_")),
        plot = neg_plot,
        width = 12, height = 8, dpi = 300, bg = "white"
      )
      
      ggsave(
        filename = file.path(output_dir, paste(comparison_name, "positive_selection.png", sep = "_")),
        plot = pos_plot,
        width = 12, height = 8, dpi = 300, bg = "white"
      )
      
    } else if (data_type == "MLE") {
      # Extract conditions and create plots for each
      conditions <- extract_mle_conditions(gene_data)
      cat("Found MLE conditions:", paste(conditions, collapse = ", "), "\n")
      
      for (condition in conditions) {
        cat("Creating plot for condition:", condition, "\n")
        
        mle_plot <- create_ranked_gene_score_plot_mle(gene_data, condition, comparison_name)
        
        # Save MLE plot
        plot_filename <- paste(comparison_name, condition, "ranked_gene_plot.png", sep = "_")
        ggsave(
          filename = file.path(output_dir, plot_filename),
          plot = mle_plot,
          width = 12, height = 8, dpi = 300, bg = "white"
        )
        
        cat("Plot saved:", plot_filename, "\n")
      }
    } else {
      cat("Warning: Could not determine data type for file:", gene_file, "\n")
    }
  }

} else if (file.exists(config_file)) {
  # Parse config
  config <- parse_config_yaml(config_file)
  experiment_name <- ifelse(is.null(config$experiment_name) || config$experiment_name == "", 
                           "experiment", config$experiment_name)
  output_dir <- config$output_dir
  method <- tolower(config$method %||% "rra")
  
  # Handle both RRA and MLE analyses
  plots_created <- FALSE
  
  # Process MLE analyses if design matrices are present
  if (!is.null(config$design_matrices)) {
    design_matrices <- config$design_matrices
    cat("Found", length(design_matrices), "MLE design matrices to process\n")
    
    for (design_matrix_config in design_matrices) {
      analysis_name <- design_matrix_config$analysis_name
      if (is.null(analysis_name)) {
        cat("Warning: analysis_name missing in design matrix, skipping\n")
        next
      }
      
      safe_analysis_name <- sanitize_filename(analysis_name)
      cat("Processing MLE analysis:", analysis_name, "\n")
      
      # Construct file paths for this analysis
      mle_dir <- file.path(output_dir, paste0(experiment_name, "_mle_", safe_analysis_name))
      gene_file <- file.path(mle_dir, paste0(experiment_name, "_", safe_analysis_name, ".gene_summary.txt"))
      
      if (!file.exists(gene_file)) {
        cat("Warning: Gene summary file not found:", gene_file, "\n")
        next
      }
      
      # Read and process MLE data
      cat("Reading MLE data from:", gene_file, "\n")
      gene_data <- read.delim(gene_file, sep = "\t", header = TRUE)
      
      # Extract conditions and create plots
      conditions <- extract_mle_conditions(gene_data)
      cat("Found conditions:", paste(conditions, collapse = ", "), "\n")
      
      for (condition in conditions) {
        cat("Creating plot for condition:", condition, "in analysis:", analysis_name, "\n")
        
        mle_plot <- create_ranked_gene_score_plot_mle(gene_data, condition, analysis_name)
        
        # Save plot in the analysis-specific directory
        plot_filename <- paste(condition, "ranked_gene_plot.png", sep = "_")
        plot_filepath <- file.path(mle_dir, plot_filename)
        
        ggsave(
          filename = plot_filepath,
          plot = mle_plot,
          width = 12, height = 8, dpi = 300, bg = "white"
        )
        
        cat("Plot saved:", plot_filepath, "\n")
        plots_created <- TRUE
      }
    }
  }
  
  # Process RRA analyses if comparisons are present
  if (!is.null(config$comparisons)) {
    comparisons <- config$comparisons 
    cat("Found", length(comparisons), "RRA comparisons to process\n")
    
    # Process each comparison
    for (comparison in comparisons) {
      # Parse comparison (assumes format like "treatment_1:DMSO_1")
      comparison_parts <- strsplit(comparison, ":")[[1]]
      treatment <- comparison_parts[1]
      control <- comparison_parts[2]
      
      # Construct file path
      comparison_name <- paste(treatment, "vs", control, sep = "_")
      gene_file <- file.path(output_dir, 
                            comparison_name,
                            paste(comparison_name,".gene_summary.txt", sep = "")
                  )
      
      # Check if file exists
      if (!file.exists(gene_file)) {
        cat("Warning: Gene summary file not found:", gene_file, "\n")
        next
      }
      
      # Read gene results
      cat("Processing RRA comparison:", gene_file, "\n")
      gene_data <- read.delim(gene_file, sep = "\t", header = TRUE)
      
      # Verify this is RRA data
      data_type <- detect_data_type(gene_data)
      if (data_type != "RRA") {
        cat("Warning: Expected RRA data but detected", data_type, "for file:", gene_file, "\n")
        next
      }
      
      # Create plots for negative selection
      neg_plot <- create_ranked_gene_score_plot_rra(gene_data, comparison_name, "neg")
      
      # Create plots for positive selection
      pos_plot <- create_ranked_gene_score_plot_rra(gene_data, comparison_name, "pos")
      
      # Save plots
      ggsave(
        filename = file.path(output_dir, paste(comparison_name, "negative_selection.png", sep = "_")),
        plot = neg_plot,
        width = 12, height = 8, dpi = 300, bg = "white"
      )
      
      ggsave(
        filename = file.path(output_dir, paste(comparison_name, "positive_selection.png", sep = "_")),
        plot = pos_plot,
        width = 12, height = 8, dpi = 300, bg = "white"
      )
      
      cat("Plots saved for RRA comparison:", comparison_name, "\n")
      plots_created <- TRUE
    }
  }
  
  # Warn if no analyses were found
  if (!plots_created) {
    cat("Warning: No MLE design matrices or RRA comparisons found in config\n")
    cat("Available config keys:", paste(names(config), collapse = ", "), "\n")
  }
}

cat("Ranked gene score plot generation complete.\n")

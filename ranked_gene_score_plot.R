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

create_ranked_gene_score_plot <- function(data, comparison_name, result_type = "neg", top_n = 20) {
  
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
  
  # Identify top genes for labeling
  top_genes <- data %>%
    filter(significance %in% c("red", "orange")) %>%
    slice_head(n = top_n)
  
  # Create plot
  p <- ggplot(data, aes(x = rank, y = get(lfc_col))) +
    geom_point(aes(color = significance), alpha = 0.7) +
    scale_color_identity() +
    geom_text(
      data = top_genes,
      aes(label = .data[[gene_id_col]]),
      size = 3
    ) +
    geom_text_repel(
      data = top_genes,
      aes(label = id),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    labs(
      title = paste(comparison_name, title_suffix, sep = " - "),
      x = "Gene Rank",
      y = y_label,
      caption = "Red: FDR < 0.05, Orange: 0.05 ≤ FDR < 0.1, Black: FDR ≥ 0.1"
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

  # Use the first path to extract output_dir and comparison_name
  first_path <- gene_files[[1]]
  path_parts <- strsplit(first_path, .Platform$file.sep)[[1]]
  output_dir <- path_parts[1]

  # Extract comparison_name from file name (assumes format: comparison_name.gene_summary.txt)
  file_base <- basename(first_path)
  comparison_name <- sub("\\.gene_summary\\.txt$", "", file_base)

  # Set up comparisons as a vector of comparison_name (for each file)
  comparisons <- sapply(gene_files, function(f) sub("\\.gene_summary\\.txt$", "", basename(f)))

}else if (file.exists(config_file)) {
  # Parse config
  config <- parse_config_yaml(config_file)
  # Get experiment info
  experiment_name <- ifelse(is.null(config$experiment_name) || config$experiment_name == "", 
                           "experiment", config$experiment_name)
  output_dir <- config$output_dir
  # Get comparisons
  comparisons <- config$comparisons 
}

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
  cat("Processing:", gene_file, "\n")
  gene_data <- read.delim(gene_file, sep = "\t", header = TRUE)
  
  # Create plots for negative selection
  neg_plot <- create_ranked_gene_score_plot(gene_data, comparison_name, "neg")
  
  # Create plots for positive selection
  pos_plot <- create_ranked_gene_score_plot(gene_data, comparison_name, "pos")
  
  # Save plots
  ggsave(
    filename = file.path(output_dir, paste(comparison_name, "negative_selection.png", sep = "_")),
    plot = neg_plot,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(output_dir, paste(comparison_name, "positive_selection.png", sep = "_")),
    plot = pos_plot,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  cat("Plots saved for comparison:", comparison_name, "\n")
}

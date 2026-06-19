# install required packages if missing
required_packages <- c("tidyverse", "yaml", "conflicted")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(yaml)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

options(width = Sys.getenv("COLUMNS", unset = 80))

parse_config <- function(file) {
    if (grepl("\\.yaml$", file) || grepl("\\.yml$", file)) {
        config_content <- readLines(file)
        config_content <- gsub(":\\s*\\.$", ": \"./\"", config_content) # Replace '.' with './' only when '.' is at the end of the line
        config <- yaml::yaml.load(paste(config_content, collapse = "\n"))
        # make all keys lowercase
        names(config) <- tolower(names(config))
    } else {
        # process as bash-like key=value pairs
        lines <- readLines(file, warn = FALSE)
        for (line in lines) {
            if (grepl("=", line) && !grepl("^#", line)) {
                key_value <- strsplit(line, "=")[[1]]
                # make all keys lowercase
                key <- tolower(key_value[1])
                # clean up
                value <- gsub("\"", "", key_value[2])
                # strip comments 
                value <- sub("#.*$", "", value)
                # trim whitespace
                value <- trimws(value)
                # check if value is bash array
                if (grepl("^\\(", value) && grepl("\\)$", value)) {
                    # remove parentheses and split by spaces not within quotes
                    value <- gsub("[()]", "", value)
                    value <- strsplit(value, "(?<!\\\\)\\s+(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", perl = TRUE)[[1]]
                }
                assign(key, value, envir = .GlobalEnv)
            }
        }
        config_values <- ls(.GlobalEnv)
        # exclude temporary variables
        config_values <- config_values[!config_values %in% c("lines", "key_value", "key", "value")] 
        # finally, build the config object
        config <- lapply(config_values, function(x) get(x, envir = .GlobalEnv))
        names(config) <- config_values
    }

    # docker specific configurations
    if (is_docker()) {
        if (!is.null(config[['docker_working_dir']])) {
            config[['working_dir']] <- normalizePath(config[['docker_working_dir']], mustWork = FALSE)
        }
        if (!is.null(config[['docker_fastq_dir']])) {
            config[['fastq_dir']] <- normalizePath(config[['docker_fastq_dir']], mustWork = FALSE)
        }
        if (!is.null(config[['docker_output_dir']])) {
            config[['output_dir']] <- normalizePath(config[['docker_output_dir']], mustWork = FALSE)
        }
    }

    # if not yaml, then resolve references like $key in values
    if (!grepl("\\.yaml$", file) && !grepl("\\.yml$", file)) {
        config <- lapply(config, function(value) resolve_references(value, config))
    }

    # handle "." values explicitly (like setting working_dir=".")
    config <- lapply(config, function(value) {
        if (is.character(value)) {
            value <- sapply(value, function(v) {
                if (v == ".") {
                    return(normalizePath(".", mustWork = FALSE)) # Interpret "." as current working directory
                }
                return(v)
            }, USE.NAMES = FALSE)
        }
        return(value)
    })

    return(config)
}


resolve_references <- function(value, config) {
    if (is.character(value)) {
        for (key in names(config)) {
            if (is.character(config[[key]])) {
                for (replacement in config[[key]]) {
                    value <- gsub(paste0("\\$", key), replacement, value)
                }
            }
        }
    }
    return(value)
}

# check if running inside a Docker container
is_docker <- function() {
    file.exists("/.dockerenv")
}

# parse config.yaml or config.sh
if (file.exists("config.yaml")) {
    config <- parse_config("config.yaml")
} else if (file.exists("config.sh")) {
    config <- parse_config("config.sh")
} else {
    stop("No configuration file found. Please make config.yaml or config.sh available in current working directory.")
}

# experiment naming
experiment_name <- ifelse(!is.null(config[['experiment_name']]), config[['experiment_name']], "")
experiment_prefix <- ""
if (experiment_name != "") {
    experiment_prefix <- paste0(experiment_name, "_")
}

# detect metadata file in order of preference: 1) env var, 2) config$metadata_file, 3) auto-detect file (yaml first, then txt)
metadata_file <- Sys.getenv("METADATA_FILE")
if (metadata_file == "" || !file.exists(metadata_file)) {
  # check config object
  if (exists("config") && !is.null(config$metadata_file)) {
    metadata_file <- config$metadata_file
  } else {
    # locate file, prefer yaml over txt
    if (file.exists("sample_metadata.yaml")) {
      metadata_file <- "sample_metadata.yaml"
    } else if (file.exists("sample_metadata.txt")) {
      metadata_file <- "sample_metadata.txt"
    } else {
      stop("Neither sample_metadata.yaml nor sample_metadata.txt found.")
    }
  }
}

if (is.null(config[['metadata_file']]) || !file.exists(config[['metadata_file']])) {
    stop("Error: metadata_file is not defined or does not exist in the configuration.")
}

# parse metadata
if (grepl("\\.yaml$|\\.yml$", metadata_file)) {
    yaml_data <- yaml::read_yaml(metadata_file)
    sample_list <- sapply(yaml_data$samples, function(x) x$sample)
    meta <- data.frame(
        library = sapply(yaml_data$samples, function(x) x$library),
        sample = sample_list,
        stringsAsFactors = FALSE
    )
} else {
    # tab-delimited metadata format
    meta <- read.delim(metadata_file, sep = "\t", header = TRUE, comment.char = "#")
}

meta <- meta %>%
    mutate(sample_original = .data$sample, # preserve original sample names
           sample = gsub("-", ".", .data$sample)) # modified sample names for R compatibility

if (is.null(config[['orig_sgrna_list_file']]) || !file.exists(config[['orig_sgrna_list_file']])) {
    stop("Error: orig_sgrna_list_file is not defined or does not exist in the configuration.")
}
sgrnas <- read.delim(config[['orig_sgrna_list_file']], sep = "\t", header = TRUE)

count_method <- config[['mode']]
outdir <- config[['output_dir']]

# determine if this is a dual sgRNA library
is_dual_sgrna_library <- function(config) {
    orig_sgrna2_list_file <- if(!is.null(config[['orig_sgrna2_list_file']])) config[['orig_sgrna2_list_file']] else config[['ORIG_SGRNA2_LIST_FILE']]
    sgrna2_list_name <- if(!is.null(config[['sgrna2_list_name']])) config[['sgrna2_list_name']] else config[['SGRNA2_LIST_NAME']]
    return(!is.null(orig_sgrna2_list_file) && !is.null(sgrna2_list_name))
}

is_dual <- is_dual_sgrna_library(config)
cat("Dual sgRNA library detected:", is_dual, "\n")

sgrnas <- read.delim(config[['orig_sgrna_list_file']], sep = "\t", header = TRUE)

if ("sample_original" %in% colnames(meta)) {
    x <- meta$sample_original[1] # Example: Use the first sample as default
} else {
    warning("'sample_original' column is missing in the metadata file. Using 'sample' column instead.")
    if ("sample" %in% colnames(meta)) {
        x <- meta$sample[1] # Use the first sample as default
    } else {
        stop("Error: Neither 'sample_original' nor 'sample' column is present in the metadata file.")
    }
}
# construct file path based on library type
if (is_dual) {
    file_path <- file.path(outdir, paste0(x, "_dual_", count_method, "/", x, "_dual_", count_method, ".countsummary.txt"))
} else {
    file_path <- file.path(outdir, paste0(x, "_", count_method, "/", x, "_", count_method, ".countsummary.txt"))
}
if (is.null(x) || x == ".") {
    stop("Error: 'x' is not defined or is invalid.")
}
if (!file.exists(file_path)) {
    stop(paste("Error: File not found -", file_path))
}
df <- read.delim(file_path, sep = "\t", header = TRUE)

count_summary <- do.call(rbind, lapply(
    meta$sample_original, function(x) {
        # construct file path based on library type
        if (is_dual) {
            count_summary_file <- paste(outdir, "/", x, "_dual_", count_method, "/", x, "_dual_", count_method, ".countsummary.txt", sep = "")
        } else {
            count_summary_file <- paste(outdir, "/", x, "_", count_method, "/", x, "_", count_method, ".countsummary.txt", sep = "")
        }
        df <- read.delim(count_summary_file, sep = "\t", header = TRUE)
        df <- df[, 1:8]
        colnames(df)[1] <- "Library" # Rename the first column to "Library"
        df$Library <- gsub("_combined.*$", "", basename(df$Library)) # Modify values to extract the prefix before "_combined"
        return(df)
    }
))
count_summary$coverage <- count_summary$Mapped / count_summary$TotalsgRNAs

cat(capture.output(print(count_summary)), sep = "\n")

##
## prepare output folder
##

plots_dir <- file.path(outdir, "count_distributions")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

read_matrix_with_prefix <- function(outdir, prefix, suffix) {
    f <- file.path(outdir, paste0(prefix, suffix))
    if (file.exists(f)) return(read.delim(f, sep = "\t", header = TRUE))
    # fallbackprefix "experiment_"
    f2 <- file.path(outdir, paste0("experiment_", suffix))
    if (file.exists(f2)) return(read.delim(f2, sep = "\t", header = TRUE))
    stop(paste("Error: matrix file not found:", f, "or", f2))
}

# set first two columns to sgRNA/gene for plotting consistency
normalize_cols <- function(df) {
    cols <- colnames(df)
    if ("sgrna_id" %in% cols) {
        cols[cols == "sgrna_id"] <- "sgRNA"
    }
    if ("sgrna_target" %in% cols) {
        cols[cols == "sgrna_target"] <- "gene"
    }
    # If still unnamed, coerce first two columns
    if (!("sgRNA" %in% cols && "gene" %in% cols)) {
        cols[1:2] <- c("sgRNA", "gene")
    }
    colnames(df) <- cols
    df
}

count_matrix_file_suffix <- "sgrna_count_matrix.txt"
cpm_matrix_file_suffix   <- "sgrna_cpm_matrix.txt"
sgrna_count.df <- read_matrix_with_prefix(outdir, experiment_prefix, count_matrix_file_suffix) %>% normalize_cols()
sgrna_cpm.df   <- read_matrix_with_prefix(outdir, experiment_prefix, cpm_matrix_file_suffix) %>% normalize_cols()

meta.2 <- left_join(meta, count_summary, by = c("sample" = "Label"))

cat("Experiment Name:", experiment_name, "\n")
cat("Experiment Prefix:", experiment_prefix, "\n")

# pivot longer
sgrna_count_long.df <- sgrna_count.df %>%
    pivot_longer(cols = -c(sgRNA, gene), names_to = "sample", values_to = "count")
sgrna_cpm_long.df <- sgrna_cpm.df %>%
    pivot_longer(cols = -c(sgRNA, gene), names_to = "sample", values_to = "cpm")

# plot count and cpm histograms
plot_count_cpm_hists <- function(count_long, cpm_long, label_prefix = "") {
    # per-sample list derived from data, robust to naming
    sample_list <- unique(count_long$sample)

    for (sample_name in sample_list) {
        # count
        p_log <- ggplot(
            count_long %>% filter(sample == sample_name) %>% mutate(log_count = log10(count + 1)),
            aes(x = log_count)
        ) +
            geom_histogram(fill = "#aaaaaa", bins = 30) +
            xlab("Log10(Count + 1)") + ylab("Frequency") +
            ggtitle(paste("Log10(Count + 1) Histogram for", sample_name)) +
            theme_minimal()

        ggsave(
            filename = file.path(plots_dir, paste0(label_prefix, sample_name, "_", count_method, "_log_histogram.png")),
            plot = p_log, width = 8, height = 5, bg = "white"
        )

        p_non_log <- ggplot(
            count_long %>% filter(sample == sample_name), aes(x = count)
        ) +
            geom_histogram(fill = "#aaaaaa", bins = 30) +
            xlab("Count") + ylab("Frequency") +
            ggtitle(paste("Count Histogram for", sample_name)) +
            theme_minimal()

        ggsave(
            filename = file.path(plots_dir, paste0(label_prefix, sample_name, "_", count_method, "_histogram.png")),
            plot = p_non_log, width = 8, height = 5, bg = "white"
        )

        # CPM
        p_log_cpm <- ggplot(
            cpm_long %>% filter(sample == sample_name) %>% mutate(log_cpm = log10(cpm + 1)),
            aes(x = log_cpm)
        ) +
            geom_histogram(fill = "#aaaaaa", bins = 30) +
            xlab("Log10(CPM + 1)") + ylab("Frequency") +
            ggtitle(paste("Log10(CPM + 1) Histogram for", sample_name)) +
            theme_minimal()

        ggsave(
            filename = file.path(plots_dir, paste0(label_prefix, sample_name, "_", count_method, "_log_cpm_histogram.png")),
            plot = p_log_cpm, width = 8, height = 5, bg = "white"
        )

        p_non_log_cpm <- ggplot(
            cpm_long %>% filter(sample == sample_name), aes(x = cpm)
        ) +
            geom_histogram(fill = "#aaaaaa", bins = 30) +
            xlab("CPM") + ylab("Frequency") +
            ggtitle(paste("CPM Histogram for", sample_name)) +
            theme_minimal()

        ggsave(
            filename = file.path(plots_dir, paste0(label_prefix, sample_name, "_", count_method, "_cpm_histogram.png")),
            plot = p_non_log_cpm, width = 8, height = 5, bg = "white"
        )
    }

    ### multi sample plots
    n_samp <- length(unique(count_long$sample))

    # count
    p_log_all <- ggplot(
        count_long %>% mutate(log_count = log10(count + 1)), aes(x = log_count)
    ) +
        geom_histogram(fill = "#aaaaaa", bins = 30) +
        xlab("Log10(Count + 1)") + ylab("Frequency") +
        ggtitle("Count Histograms (Log)") + theme_minimal() +
        facet_wrap(~ sample, ncol = 1)

    ggsave(
        filename = file.path(plots_dir, paste0(label_prefix, count_method, "_log_histogram_all_samples.png")),
        plot = p_log_all, width = 3, height = 1 * n_samp, bg = "white"
    )

    p_non_log_all <- ggplot(count_long, aes(x = count)) +
        geom_histogram(fill = "#aaaaaa", bins = 30) +
        xlab("Count") + ylab("Frequency") + ggtitle("Count Histograms") + theme_minimal() +
        facet_wrap(~ sample, ncol = 1)

    ggsave(
        filename = file.path(plots_dir, paste0(label_prefix, count_method, "_histogram_all_samples.png")),
        plot = p_non_log_all, width = 3, height = 1 * n_samp, bg = "white"
    )

    # CPM
    p_log_cpm_all <- ggplot(
        cpm_long %>% mutate(log_cpm = log10(cpm + 1)), aes(x = log_cpm)
    ) +
        geom_histogram(fill = "#aaaaaa", bins = 30) +
        xlab("Log10(CPM + 1)") + ylab("Frequency") + ggtitle("CPM Histograms (Log)") + theme_minimal() +
        facet_wrap(~ sample, ncol = 1)

    ggsave(
        filename = file.path(plots_dir, paste0(label_prefix, count_method, "_log_cpm_histogram_all_samples.png")),
        plot = p_log_cpm_all, width = 3, height = 1 * n_samp, bg = "white"
    )

    p_non_log_cpm_all <- ggplot(cpm_long, aes(x = cpm)) +
        geom_histogram(fill = "#aaaaaa", bins = 30) +
        xlab("CPM") + ylab("Frequency") + ggtitle("CPM Histograms") + theme_minimal() +
        facet_wrap(~ sample, ncol = 1)

    ggsave(
        filename = file.path(plots_dir, paste0(label_prefix, count_method, "_cpm_histogram_all_samples.png")),
        plot = p_non_log_cpm_all, width = 3, height = 1 * n_samp, bg = "white"
    )
}

# plots for standard matrices (non-dual)
plot_count_cpm_hists(sgrna_count_long.df, sgrna_cpm_long.df, label_prefix = "")

# if dual sgrna library, process pg_count/pg_cpm matrices
if (is_dual) {
    pg_count_suffix <- "sgrna_pg_count_matrix.txt"
    pg_cpm_suffix   <- "sgrna_pg_cpm_matrix.txt"
    pg_count_path <- file.path(outdir, paste0(experiment_prefix, pg_count_suffix))
    pg_cpm_path   <- file.path(outdir, paste0(experiment_prefix, pg_cpm_suffix))

    if (!file.exists(pg_count_path)) {
        # Try fallback prefix
        pg_count_path2 <- file.path(outdir, paste0("experiment_", pg_count_suffix))
        if (file.exists(pg_count_path2)) pg_count_path <- pg_count_path2
    }
    if (!file.exists(pg_cpm_path)) {
        pg_cpm_path2 <- file.path(outdir, paste0("experiment_", pg_cpm_suffix))
        if (file.exists(pg_cpm_path2)) pg_cpm_path <- pg_cpm_path2
    }

    if (file.exists(pg_count_path) && file.exists(pg_cpm_path)) {
        pg_count.df <- read.delim(pg_count_path, sep = "\t", header = TRUE) %>% normalize_cols()
        pg_cpm.df   <- read.delim(pg_cpm_path, sep = "\t", header = TRUE) %>% normalize_cols()

        pg_count_long.df <- pg_count.df %>%
            pivot_longer(cols = -c(sgRNA, gene), names_to = "sample", values_to = "count")
        pg_cpm_long.df <- pg_cpm.df %>%
            pivot_longer(cols = -c(sgRNA, gene), names_to = "sample", values_to = "cpm")

        # "pg_" prefix for plots
        plot_count_cpm_hists(pg_count_long.df, pg_cpm_long.df, label_prefix = "pg_")
    } else {
        warning("Dual library detected but pg_count/pg_cpm matrices not found; skipping pg plots.")
    }
}


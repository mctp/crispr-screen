library(tidyverse)
library(yaml)

options(width = Sys.getenv("COLUMNS", unset = 80))
# Function to parse config file and set environment variables
library(conflicted)

parse_config <- function(file) {
    if (grepl("\\.yaml$", file) || grepl("\\.yml$", file)) {
        config_content <- readLines(file)
        config_content <- gsub(":\\s*\\.$", ": \"./\"", config_content) # Replace '.' with './' only when '.' is at the end of the line
        config <- yaml::yaml.load(paste(config_content, collapse = "\n"))
        names(config) <- tolower(names(config)) # Convert keys to lowercase
    } else {
        lines <- readLines(file, warn = FALSE)
        for (line in lines) {
            if (grepl("=", line) && !grepl("^#", line)) {
                key_value <- strsplit(line, "=")[[1]]
                key <- tolower(key_value[1]) # Convert key to lowercase
                value <- gsub("\"", "", key_value[2])
                value <- sub("#.*$", "", value) # Remove comments
                value <- trimws(value) # Remove leading and trailing whitespace
                # Check if the value is a bash array
                if (grepl("^\\(", value) && grepl("\\)$", value)) {
                    value <- gsub("[()]", "", value)
                    value <- strsplit(value, "(?<!\\\\)\\s+(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", perl = TRUE)[[1]]
                }
                assign(key, value, envir = .GlobalEnv)
            }
        }
        config_values <- ls(.GlobalEnv)
        config_values <- config_values[!config_values %in% c("lines", "key_value", "key", "value")] # Exclude temporary variables
        config <- lapply(config_values, function(x) get(x, envir = .GlobalEnv))
        names(config) <- config_values
    }

    # Check for Docker status and set directories accordingly
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

    if (!grepl("\\.yaml$", file) && !grepl("\\.yml$", file)) {
        config <- lapply(config, function(value) resolve_references(value, config))
    }

    # Handle "." values explicitly
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

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

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

# Check if running inside a Docker container
is_docker <- function() {
    file.exists("/.dockerenv")
}

# Parse the config.yaml file if it exists, otherwise try config.sh
if (file.exists("config.yaml")) {
    config <- parse_config("config.yaml")
} else if (file.exists("config.sh")) {
    config <- parse_config("config.sh")
} else {
    stop("No configuration file found. Please make config.yaml or config.sh available in current working directory.")
}


experiment_name <- ifelse(!is.null(config[['experiment_name']]), config[['experiment_name']], "")
experiment_prefix <- ""
if (experiment_name != "") {
    experiment_prefix <- paste0(experiment_name, "_")
}

if (is.null(config[['metadata_file']]) || !file.exists(config[['metadata_file']])) {
    stop("Error: metadata_file is not defined or does not exist in the configuration.")
}
meta <- read.delim(config[['metadata_file']], sep = "\t", header = TRUE, comment.char = "#")
meta <- meta %>%
    mutate(sample_original = sample, # Preserve original sample names
           sample = gsub("-", ".", sample)) # Create a modified version for other operations

if (is.null(config[['orig_sgrna_list_file']]) || !file.exists(config[['orig_sgrna_list_file']])) {
    stop("Error: orig_sgrna_list_file is not defined or does not exist in the configuration.")
}
sgrnas <- read.delim(config[['orig_sgrna_list_file']], sep = "\t", header = TRUE)

count_method <- config[['mode']]
outdir <- config[['output_dir']]

sgrnas <- read.delim(config[['orig_sgrna_list_file']], sep = "\t", header = TRUE)

# Define 'x' based on the metadata file
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
file_path <- file.path(outdir, paste0(meta$sample_original[1], "_", count_method, "/", meta$sample_original[1], "_", count_method, ".countsummary.txt"))
file_path <- file.path(outdir, paste0(x, "_", count_method, "/", x, "_", count_method, ".countsummary.txt"))
if (is.null(x) || x == ".") {
    stop("Error: 'x' is not defined or is invalid.")
}
if (!file.exists(file_path)) {
    stop(paste("Error: File not found -", file_path))
}
df <- read.delim(file_path, sep = "\t", header = TRUE)

count_summary <- do.call(rbind, lapply(
    meta$sample_original, function(x) {
        df <- read.delim(paste(outdir, "/", x, "_", count_method, "/", x, "_", count_method, ".countsummary.txt", sep = ""), sep = "\t", header = TRUE)
        df <- df[, 1:8]
        colnames(df)[1] <- "Library" # Rename the first column to "Library"
        df$Library <- gsub("_combined.*$", "", basename(df$Library)) # Modify values to extract the prefix before "_combined"
        return(df)
    }
))
count_summary$coverage <- count_summary$Mapped / count_summary$TotalsgRNAs

cat(capture.output(print(count_summary)), sep = "\n")
# Load counts and CPM matrices
count_matrix_file <- file.path(outdir, paste0(experiment_prefix, "sgrna_count_matrix.txt"))
cpm_matrix_file <- file.path(outdir, paste0(experiment_prefix, "sgrna_cpm_matrix.txt"))

if (!file.exists(count_matrix_file) || !file.exists(cpm_matrix_file)) {
    # Check with "experiment_" prefix if files are not found
    experiment_prefix <- "experiment_"
    count_matrix_file <- file.path(outdir, paste0(experiment_prefix, "sgrna_count_matrix.txt"))
    cpm_matrix_file <- file.path(outdir, paste0(experiment_prefix, "sgrna_cpm_matrix.txt"))
}

if (!file.exists(count_matrix_file)) {
    stop(paste("Error: Count matrix file not found -", count_matrix_file))
}
if (!file.exists(cpm_matrix_file)) {
    stop(paste("Error: CPM matrix file not found -", cpm_matrix_file))
}
sgrna_cpm.df <- read.delim(cpm_matrix_file, sep = "\t", header = TRUE)
meta.2 <- left_join(meta,count_summary,by=c("sample"="Label"))

# Print experiment name and prefix for debugging
cat("Experiment Name:", experiment_name, "\n")
cat("Experiment Prefix:", experiment_prefix, "\n")
# Load counts and CPM matrices
count_matrix_file <- file.path(outdir, paste0(experiment_prefix, "sgrna_count_matrix.txt"))
cpm_matrix_file <- file.path(outdir, paste0(experiment_prefix, "sgrna_cpm_matrix.txt"))

sgrna_count.df <- read.delim(count_matrix_file, sep = "\t", header = TRUE)
colnames(sgrna_count.df)[1:2] <- c("sgRNA", "gene")

sgrna_cpm.df <- read.delim(cpm_matrix_file, sep = "\t", header = TRUE)
colnames(sgrna_cpm.df)[1:2] <- c("sgRNA", "gene")

# Pivot long the sgrna_count.df
sgrna_count_long.df <- sgrna_count.df %>%
    pivot_longer(
        cols = -c(sgRNA, gene),
        names_to = "sample",
        values_to = "count"
    )

# Pivot long the sgrna_cpm.df
sgrna_cpm_long.df <- sgrna_cpm.df %>%
    pivot_longer(
        cols = -c(sgRNA, gene),
        names_to = "sample",
        values_to = "cpm"
    )


# Create and save histogram plots for each sample
for (sample_name in unique(meta$sample)) {

    ## Count histograms

    # Log scale histogram including zeros
    p_log <- ggplot(
        sgrna_count_long.df %>% 
        filter(sample == sample_name) %>% 
        mutate(log_count = log10(count + 1)),
        aes(x = log_count)
    ) + 
    geom_histogram(fill = "#aaaaaa", bins = 30) +
    xlab("Log10(Count + 1)") +
    ylab("Frequency") +
    ggtitle(paste("Log10(Count + 1) Histogram for", sample_name)) +
    theme_minimal()
    
    # Save the log scale histogram
    ggsave(
        filename = file.path(outdir, paste0(sample_name, "_", count_method, "_log_histogram.png")),
        plot = p_log,
        width = 8,
        height = 5,
        bg = "white"
    )
    
    # Non-log scale histogram
    p_non_log <- ggplot(
        sgrna_count_long.df %>% 
        filter(sample == sample_name),
        aes(x = count)
    ) + 
    geom_histogram(fill = "#aaaaaa", bins = 30) +
    xlab("Count") +
    ylab("Frequency") +
    ggtitle(paste("Count Histogram for", sample_name)) +
    theme_minimal()
    
    # Save the non-log scale histogram
    ggsave(
        filename = file.path(outdir, paste0(sample_name, "_", count_method, "_histogram.png")),
        plot = p_non_log,
        width = 8,
        height = 5,
        bg = "white"
    )

    ## CPM histograms

    # Log scale histogram including zeros
    p_log_cpm <- ggplot(
        sgrna_cpm_long.df %>% 
        filter(sample == sample_name) %>% 
        mutate(log_cpm = log10(cpm + 1)),
        aes(x = log_cpm)
    ) + 
    geom_histogram(fill = "#aaaaaa", bins = 30) +
    xlab("Log10(CPM + 1)") +
    ylab("Frequency") +
    ggtitle(paste("Log10(CPM + 1) Histogram for", sample_name)) +
    theme_minimal()
    
    # Save the log scale histogram
    ggsave(
        filename = file.path(outdir, paste0(sample_name, "_", count_method, "_log_cpm_histogram.png")),
        plot = p_log_cpm,
        width = 8,
        height = 5,
        bg = "white"
    )
    
    # Non-log scale histogram
    p_non_log_cpm <- ggplot(
        sgrna_cpm_long.df %>% 
        filter(sample == sample_name),
        aes(x = cpm)
    ) + 
    geom_histogram(fill = "#aaaaaa", bins = 30) +
    xlab("CPM") +
    ylab("Frequency") +
    ggtitle(paste("CPM Histogram for", sample_name)) +
    theme_minimal()
    
    # Save the non-log scale histogram
    ggsave(
        filename = file.path(outdir, paste0(sample_name, "_", count_method, "_cpm_histogram.png")),
        plot = p_non_log_cpm,
        width = 8,
        height = 5,
        bg = "white"
    )
}

# Log scale count histogram for all samples
p_log_all <- ggplot(
    sgrna_count_long.df %>% 
    mutate(log_count = log10(count + 1)),
    aes(x = log_count)
) + 
geom_histogram(fill = "#aaaaaa", bins = 30) +
xlab("Log10(Count + 1)") +
ylab("Frequency") +
ggtitle("Count Histograms (Log)") +
theme_minimal() +
facet_wrap(~ sample, ncol = 1)

ggsave(
    filename = file.path(outdir, paste0(count_method, "_log_histogram_all_samples.png")),
    plot = p_log_all,
    width = 3,
    height = 1 * length(unique(meta$sample)),
    bg = "white"
)

# Non-log scale count histogram for all samples
p_non_log_all <- ggplot(
    sgrna_count_long.df,
    aes(x = count)
) + 
geom_histogram(fill = "#aaaaaa", bins = 30) +
xlab("Count") +
ylab("Frequency") +
ggtitle("Count Histograms") +
theme_minimal() +
facet_wrap(~ sample, ncol = 1)

ggsave(
    filename = file.path(outdir, paste0(count_method, "_histogram_all_samples.png")),
    plot = p_non_log_all,
    width = 3,
    height = 1 * length(unique(meta$sample)),
    bg = "white"
)

# Log scale CPM histogram for all samples
p_log_cpm_all <- ggplot(
    sgrna_cpm_long.df %>% 
    mutate(log_cpm = log10(cpm + 1)),
    aes(x = log_cpm)
) + 
geom_histogram(fill = "#aaaaaa", bins = 30) +
xlab("Log10(CPM + 1)") +
ylab("Frequency") +
ggtitle("CPM Histograms (Log)") +
theme_minimal() +
facet_wrap(~ sample, ncol = 1)

ggsave(
    filename = file.path(outdir, paste0(count_method, "_log_cpm_histogram_all_samples.png")),
    plot = p_log_cpm_all,
    width = 3,
    height = 1 * length(unique(meta$sample)),
    bg = "white"
)

# Non-log scale CPM histogram for all samples
p_non_log_cpm_all <- ggplot(
    sgrna_cpm_long.df,
    aes(x = cpm)
) + 
geom_histogram(fill = "#aaaaaa", bins = 30) +
xlab("CPM") +
ylab("Frequency") +
ggtitle("CPM Histograms") +
theme_minimal() +
facet_wrap(~ sample, ncol = 1)

ggsave(
    filename = file.path(outdir, paste0(count_method, "_cpm_histogram_all_samples.png")),
    plot = p_non_log_cpm_all,
    width = 3,
    height = 1 * length(unique(meta$sample)),
    bg = "white"
)


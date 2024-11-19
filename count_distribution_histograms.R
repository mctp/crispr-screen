library(tidyverse)

options(width = Sys.getenv("COLUMNS", unset = 80))

# Function to parse config.sh and set environment variables
parse_config <- function(file) {
    lines <- readLines(file, warn = FALSE)
    for (line in lines) {
        if (grepl("=", line) && !grepl("^#", line)) {
            key_value <- strsplit(line, "=")[[1]]
            key <- key_value[1]
            value <- gsub("\"", "", key_value[2])
            value <- sub("#.*$", "", value) # Remove comments
            value <- trimws(value) # Remove leading and trailing whitespace
            # Check if the value is a bash array
            if (grepl("^\\(", value) && grepl("\\)$", value)) {
                value <- gsub("[()]", "", value)
                # Split the string 'value' by whitespace that is not preceded by a backslash and is outside of quoted substrings.
                # The regex breakdown:
                # - (?<!\\\\): Negative lookbehind to ensure the whitespace is not preceded by a backslash.
                # - \\s+: Matches one or more whitespace characters.
                # - (?= ... ): Positive lookahead to ensure the following condition is met:
                #   - (?:[^\"]*\"[^\"]*\")*: Non-capturing group that matches pairs of quotes with any characters except quotes in between.
                #   - [^\"]*$: Matches any characters except quotes until the end of the string.
                # The 'perl = TRUE' argument enables Perl-compatible regex for advanced features like lookaheads and lookbehinds.
                value <- strsplit(value, "(?<!\\\\)\\s+(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", perl = TRUE)[[1]]
            }
            assign(key, value, envir = .GlobalEnv)
        }
    }
    config_values <- ls(.GlobalEnv)
    config_values <- config_values[!config_values %in% c("lines", "key_value", "key", "value")] # Exclude temporary variables
    config <- lapply(config_values, function(x) get(x, envir = .GlobalEnv))
    names(config) <- config_values

    # Check for Docker status and set directories accordingly
    if (is_docker()) {
        if (!is.null(config[['DOCKER_WORKING_DIR']])) {
            config[['WORKING_DIR']] <- config[['DOCKER_WORKING_DIR']]
        }
        if (!is.null(config[['DOCKER_FASTQ_DIR']])) {
            config[['FASTQ_DIR']] <- config[['DOCKER_FASTQ_DIR']]
        }
        if (!is.null(config[['DOCKER_OUTPUT_DIR']])) {
            config[['OUTPUT_DIR']] <- config[['DOCKER_OUTPUT_DIR']]
        }
    }

    config <- lapply(config, function(value) resolve_references(value, config))
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

# Check if running inside a Docker container
is_docker <- function() {
    file.exists("/.dockerenv")
}

# Parse the config.sh file
config <- parse_config("config.sh")


experiment_name <- ifelse(!is.null(config[['experiment_name']]), config[['experiment_name']], "")
experiment_prefix <- ""
if (experiment_name != "") {
    experiment_prefix <- paste0(experiment_name, "_")
}

meta <- read.delim(config[['METADATA_FILE']], sep = "\t", header = TRUE, comment.char = "#")
meta <- meta %>%
    mutate(sample_original = sample) %>%
    mutate(sample = gsub("-", ".", sample))

sgrnas <- read.delim(config[['ORIG_SGRNA_LIST_FILE']],sep="\t",header=TRUE)

count_method <- config[['MODE']]
outdir <- config[['OUTPUT_DIR']]

count_summary <- bind_rows(lapply(
    meta$sample_original, function(x) {
        df <- read.delim(paste(outdir, "/", x, "_", count_method, "/", x, "_", count_method, ".countsummary.txt", sep = ""), sep = "\t", header = TRUE)
        return(df[, 1:8])
    }))
count_summary$coverage <- count_summary$Mapped / count_summary$TotalsgRNAs

cat(capture.output(print(count_summary)), sep = "\n")

meta.2 <- left_join(meta,count_summary,by=c("sample"="Label"))

# Load counts and CPM matrices
sgrna_count.df <- read.delim(file.path(outdir, "sgrna_count_matrix.txt"), sep = "\t", header = TRUE)
colnames(sgrna_count.df)[1:2] <- c("sgRNA", "gene")
sgrna_cpm.df <- read.delim(file.path(outdir, "sgrna_cpm_matrix.txt"), sep = "\t", header = TRUE)
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


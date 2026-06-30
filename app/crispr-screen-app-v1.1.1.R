#!/usr/bin/env Rscript

# A shiny app to explore CRISPR screen results.
# Brian Magnuson (bmagnuso@umich.edu) -- Michigan Center for Translational Pathology

library(shiny)
library(tidyverse)
library(yaml)
library(ggrepel)
library(colourpicker)
library(shinyFiles)

# Use Cairo for bitmap rendering whenever available.
# The default (non-Cairo) PNG device on Linux fails with "failed to load
# default encoding" when re-opened with a different height (e.g. after
# adding a plot). Cairo is more robust and doesn't require afm font files.
if (capabilities("cairo")) {
    options(bitmapType = "cairo")
}

PLOTLY_AVAILABLE <- requireNamespace("plotly", quietly = TRUE)
if (PLOTLY_AVAILABLE) library(plotly)

DT_AVAILABLE <- requireNamespace("DT", quietly = TRUE)
if (DT_AVAILABLE) library(DT)

OPENXLSX2_AVAILABLE <- requireNamespace("openxlsx2", quietly = TRUE)
OPENXLSX_AVAILABLE  <- requireNamespace("openxlsx",  quietly = TRUE)
GGEXTRA_AVAILABLE   <- requireNamespace("ggExtra",   quietly = TRUE)

# Global UI constants
PICKER_WIDTH <- "60px"
PALETTE_WIDTH <- "95px"
DEFAULT_FDR_COLOR_1 <- "#FF0000"
DEFAULT_FDR_COLOR_2 <- "#FF7F00"

`%||%` <- function(lhs, rhs) {
    if (!is.null(lhs) && length(lhs) > 0) lhs else rhs
}

is_docker <- function() file.exists("/.dockerenv")

parse_config <- function(file) {
    if (grepl("\\.yaml$", file) || grepl("\\.yml$", file)) {
        config_content <- readLines(file)
        config_content <- gsub(":\\s*\\.$", ": \"./\"", config_content)
        config <- yaml::yaml.load(paste(config_content, collapse = "\n"))
        names(config) <- tolower(names(config))
    } else {
        local_env <- new.env(parent = emptyenv())
        lines <- readLines(file, warn = FALSE)
        for (line in lines) {
            if (grepl("=", line) && !grepl("^#", line)) {
                parts <- strsplit(line, "=")[[1]]
                key   <- trimws(parts[1])
                value <- trimws(parts[2])
                assign(key, value, envir = local_env)
            }
        }
        config_values <- ls(local_env)
        config <- lapply(config_values, function(x) get(x, envir = local_env))
        names(config) <- config_values
    }
    if (is_docker()) {
        if (!is.null(config[['docker_working_dir']])) config[['working_dir']] <- normalizePath(config[['docker_working_dir']], mustWork = FALSE)
        if (!is.null(config[['docker_fastq_dir']])) config[['fastq_dir']] <- normalizePath(config[['docker_fastq_dir']], mustWork = FALSE)
        if (!is.null(config[['docker_output_dir']])) config[['output_dir']] <- normalizePath(config[['docker_output_dir']], mustWork = FALSE)
    }
    return(config)
}

extract_comparison_string <- function(comparison_entry) {
    if (is.character(comparison_entry)) return(comparison_entry)
    else if (is.list(comparison_entry) && !is.null(comparison_entry$comparison)) return(comparison_entry$comparison)
    else stop("Invalid comparison format in config")
}

extract_plot_config <- function(comparison_entry) {
    if (is.list(comparison_entry) && !is.null(comparison_entry$plot_config)) return(comparison_entry$plot_config)
    return(NULL)
}

parse_comparison_with_replicates <- function(comparison_str) {
    parts <- strsplit(comparison_str, ":")[[1]]
    if (length(parts) != 2) stop(paste("Invalid comparison format:", comparison_str))
    treatment_samples <- trimws(strsplit(parts[1], ",")[[1]])
    control_samples <- trimws(strsplit(parts[2], ",")[[1]])
    treatment_name <- generate_group_name(treatment_samples)
    control_name <- generate_group_name(control_samples)
    list(treatment_samples = treatment_samples, control_samples = control_samples,
         comparison_name = paste0(treatment_name, "_vs_", control_name))
}

generate_group_name <- function(sample_list) {
    if (length(sample_list) == 1) return(sample_list[1])
    base_names <- character(length(sample_list))
    for (i in seq_along(sample_list)) {
        sample <- sample_list[i]
        if (grepl("_\\d+$", sample)) base_names[i] <- sub("_\\d+$", "", sample)
        else base_names[i] <- NA_character_
    }
    if (all(!is.na(base_names)) && length(unique(base_names)) == 1) return(paste0(base_names[1], "_reps"))
    return(paste0(sample_list[1], "_reps"))
}

find_rra_results <- function(output_dir, experiment_name, comparisons) {
    result_files <- list()
    for (comparison_entry in comparisons) {
        comparison_str <- extract_comparison_string(comparison_entry)
        plot_config <- extract_plot_config(comparison_entry)
        parsed <- parse_comparison_with_replicates(comparison_str)
        comparison_name <- parsed$comparison_name
        rra_dir <- file.path(output_dir, paste0(experiment_name, "_rra_", comparison_name))
        gene_file <- file.path(rra_dir, paste0(comparison_name, ".gene_summary.txt"))
        if (file.exists(gene_file)) {
            result_files[[comparison_str]] <- list(
                comparison_name = comparison_name, gene_file = gene_file,
                output_dir = rra_dir, original_comparison = comparison_str, plot_config = plot_config)
            cat("Found RRA results for:", comparison_str, "->", comparison_name, "\n")
        } else {
            cat("Warning: No RRA results found for:", comparison_str, "at", gene_file, "\n")
        }
    }
    return(result_files)
}

load_gene_list <- function(file_path) {
    if (!file.exists(file_path)) { cat("Warning: Gene list file not found:", file_path, "\n"); return(character(0)) }
    genes <- readLines(file_path, warn = FALSE)
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    # Strip common single-word column headers (gene, symbol, name, etc.)
    if (length(genes) > 0 && tolower(genes[1]) %in% c("gene","genes","symbol","name","geneid","gene_id","hgnc_symbol","gene_symbol","genename"))
        genes <- genes[-1]
    genes
}

# Returns a named list of system gene group definitions derived from config.
# Keys: "neg_ctrl", "essential", "nonessential" (only present when defined/applicable).
# neg_ctrl is ALWAYS included (so the edit link is visible) but ctrl_id may be NULL.
get_system_gene_groups <- function(config, config_dir, neg_ctrl_override = NULL) {
    groups <- list()
    ctrl_id <- if (!is.null(neg_ctrl_override) && nzchar(trimws(neg_ctrl_override %||% "")))
                   trimws(neg_ctrl_override)
               else if (!is.null(config$control_gene_id) && nzchar(trimws(config$control_gene_id %||% "")))
                   trimws(config$control_gene_id)
               else NULL
    groups[["neg_ctrl"]] <- list(key="neg_ctrl", group_name="Negative Controls",
                                  default_color="#E41A1C", ctrl_id=ctrl_id,
                                  genes=ctrl_id, is_custom=TRUE)
    resolve_file <- function(f) {
        if (is.null(f) || !nzchar(trimws(f %||% ""))) return(NULL)
        fp <- trimws(f)
        if (!is.null(config_dir) && !file.exists(fp)) fp <- file.path(config_dir, fp)
        fp
    }
    fp <- resolve_file(config$core_essential_genes_file)
    if (!is.null(fp))
        groups[["essential"]] <- list(key="essential", group_name="Core Essential",
                                       default_color="#4DAF4A", file=fp, is_custom=FALSE)
    fp <- resolve_file(config$non_essential_genes_file)
    if (!is.null(fp))
        groups[["nonessential"]] <- list(key="nonessential", group_name="Non-Essential",
                                          default_color="#984EA3", file=fp, is_custom=FALSE)
    groups
}

SYS_GROUP_DEFAULTS  <- c(neg_ctrl="#1F77B4", essential="#9467BD", nonessential="#2CA02C")
GENESET_COLOR_CYCLE <- c("#17BECF", "#E377C2", "#8C564B", "#BCBD22", "#56B4E9", "#009E73")

# Render one system-group row (checkbox + name/edit-link + color picker).
# tab_prefix: "counts" or "comp"
make_sys_group_row <- function(tab_prefix, sg) {
    key    <- sg$key
    pid    <- paste0(tab_prefix, "_sys_", key)
    name_col <- if (key == "neg_ctrl") {
        id_text <- if (!is.null(sg$ctrl_id) && nzchar(sg$ctrl_id)) sg$ctrl_id else "not configured"
        tagList(
            tags$span(sg$group_name, style="font-weight:600;font-size:0.9em;"),
            tags$span(paste0("(", id_text, ")"), style="font-size:0.78em;color:#888;margin-left:4px;"),
            actionLink(paste0(tab_prefix, "_neg_ctrl_edit"), "\u270e",
                       style="font-size:0.9em;margin-left:3px;color:#666;", title="Set control gene ID")
        )
    } else {
        tags$span(sg$group_name, style="font-weight:600;font-size:0.9em;")
    }
    div(style="border:1px solid #ddd;padding:5px;margin-bottom:5px;border-radius:3px;",
        fluidRow(
            column(3, checkboxInput(pid, label=NULL, value=FALSE)),
            column(9, name_col)
        ),
        fluidRow(
            column(6, div(class="palette-dropdown",
                selectInput(paste0(pid, "_palette"), NULL,
                    choices=c("Basic"="limited","Web-safe"="websafe","Full"="square"),
                    selected="square", width="100%"))),
            column(6, uiOutput(paste0(pid, "_color_picker")))
        )
    )
}

render_plot_with_fallback <- function(ggplot_obj, use_plotly = TRUE, tooltip = "text") {
    if (PLOTLY_AVAILABLE && use_plotly) {
        ggplot_obj <- ggplot_obj + guides(size = "none")
        p <- ggplotly(ggplot_obj, tooltip = tooltip) %>% layout(hovermode = "closest")
        # ggplotly names each trace by the combination of all discrete aesthetics,
        # e.g. "(FDR < 0.05,NA)" or "(FDR < 0.05,Group1)". Clean these up so the
        # legend shows one clean entry per significance level or gene group.
        seen <- character(0)
        p$x$data <- lapply(p$x$data, function(tr) {
            nm <- tr$name
            if (is.null(nm) || !nzchar(nm)) return(tr)
            # Pattern: "(fill_value,color_value)"
            if (grepl("^\\(.*,.*\\)$", nm)) {
                parts <- regmatches(nm, regexec("^\\((.+),(.+)\\)$", nm))[[1]]
                if (length(parts) == 3) {
                    fill_val  <- parts[2]
                    color_val <- parts[3]
                    if (color_val == "NA") {
                        # Background point — label by significance level
                        tr$name <- fill_val
                    } else {
                        # Highlighted point — label by gene group, deduplicate
                        if (color_val %in% seen) {
                            tr$showlegend <- FALSE
                        } else {
                            seen <<- c(seen, color_val)
                            tr$name <- color_val
                        }
                    }
                }
            }
            tr
        })
        return(p)
    }
    return(ggplot_obj)
}

# ── QC Helper Functions ────────────────────────────────────────────────────────

parse_histogram_metadata <- function(filepath) {
    fname <- tools::file_path_sans_ext(basename(filepath))
    # All-samples combined files: {mode}_log_cpm_histogram_all_samples.png etc.
    if (grepl("_histogram_all_samples$", fname, ignore.case = TRUE)) {
        prefix <- sub("_histogram_all_samples$", "", fname, ignore.case = TRUE)
        if (grepl("_log_cpm$", prefix, ignore.case = TRUE)) type <- "Log CPM"
        else if (grepl("_cpm$", prefix, ignore.case = TRUE)) type <- "CPM"
        else if (grepl("_log$", prefix, ignore.case = TRUE)) type <- "Log Count"
        else type <- "Count"
        return(list(sample = "All Samples", type = type, filepath = filepath))
    }
    # Individual: {sample}_{mode}_log_cpm_histogram.png etc.
    if (grepl("_log_cpm_histogram$", fname, ignore.case = TRUE)) {
        type <- "Log CPM"; fname <- sub("_log_cpm_histogram$", "", fname, ignore.case = TRUE)
    } else if (grepl("_cpm_histogram$", fname, ignore.case = TRUE)) {
        type <- "CPM"; fname <- sub("_cpm_histogram$", "", fname, ignore.case = TRUE)
    } else if (grepl("_log_histogram$", fname, ignore.case = TRUE)) {
        type <- "Log Count"; fname <- sub("_log_histogram$", "", fname, ignore.case = TRUE)
    } else {
        type <- "Count"; fname <- sub("_histogram$", "", fname, ignore.case = TRUE)
    }
    # Strip trailing mode suffix (fastq, bam, etc.)
    fname <- sub("_(fastq|bam)$", "", fname, ignore.case = TRUE)
    list(sample = fname, type = type, filepath = filepath)
}

find_qc_files <- function(output_dir) {
    if (!dir.exists(output_dir)) return(list(countsummary = c(), histograms = c(), pca = c(), barplots = c()))

    # Histograms may be directly in output_dir (older pipeline) or in count_distributions/ subdir (newer pipeline)
    histograms <- list.files(output_dir, pattern = "_histogram\\.png$", full.names = TRUE)
    count_dist_dir <- file.path(output_dir, "count_distributions")
    if (dir.exists(count_dist_dir))
        histograms <- c(histograms, list.files(count_dist_dir, pattern = "_histogram\\.png$", full.names = TRUE))
    barplots   <- list.files(output_dir, pattern = "_barplot\\.png$",   full.names = TRUE)

    # Countsummary files live one level deep inside per-sample subdirs
    subdirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
    countsummary <- unlist(lapply(subdirs, function(d)
        list.files(d, pattern = "\\.countsummary\\.txt$", full.names = TRUE)))
    # Fallback: also check output_dir itself
    countsummary <- c(countsummary,
                      list.files(output_dir, pattern = "\\.countsummary\\.txt$", full.names = TRUE))
    countsummary <- unique(countsummary)

    # PCA: check known PCA_variance subdir first; otherwise scan immediate subdirs only
    pca_subdir <- file.path(output_dir, "PCA_variance")
    if (dir.exists(pca_subdir)) {
        pca <- list.files(pca_subdir, pattern = "pca.*\\.png$", full.names = TRUE, ignore.case = TRUE)
    } else {
        pca <- c(
            list.files(output_dir, pattern = "pca.*\\.png$", full.names = TRUE, ignore.case = TRUE),
            unlist(lapply(subdirs, function(d)
                list.files(d, pattern = "pca.*\\.png$", full.names = TRUE, ignore.case = TRUE)))
        )
        pca <- unique(pca)
    }

    list(countsummary = countsummary, histograms = histograms, pca = pca, barplots = barplots)
}

load_countsummary_data <- function(files) {
    if (length(files) == 0) return(NULL)
    dfs <- lapply(files, function(f) tryCatch(
        read.delim(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) { cat("Warning: Could not read countsummary:", f, "\n"); NULL }
    ))
    dfs <- Filter(Negate(is.null), dfs)
    if (length(dfs) == 0) return(NULL)
    tryCatch(bind_rows(dfs), error = function(e) dfs[[1]])
}

# ── Experiment discovery helpers ─────────────────────────────────────────────

# Derive a human-readable label from a config.yaml file path.
# Expected convention: <experiment>/crispr-screen[-modifier]/config.yaml
#   "my_exp/crispr-screen-a/config.yaml"  →  "my exp (a)"
#   "my_exp/crispr-screen/config.yaml"    →  "my exp"
# Non-standard paths: last 2-3 dirs joined with " - ".
config_path_to_label <- function(path) {
    path <- gsub("\\\\", "/", path)
    parts <- strsplit(path, "/")[[1]]
    parts <- parts[nzchar(parts)]
    parts <- parts[!grepl("^config\\.ya?ml$", parts, ignore.case = TRUE)]
    if (length(parts) == 0) return(basename(path))

    cs_idx <- grep("^crispr-screen(-.*)?$", parts, ignore.case = TRUE)
    if (length(cs_idx) > 0) {
        i    <- cs_idx[length(cs_idx)]
        mod  <- sub("^crispr-screen-?", "", parts[i], ignore.case = TRUE)
        name <- if (i > 1) gsub("[_-]", " ", parts[i - 1]) else "experiment"
        if (nzchar(trimws(mod))) paste0(name, " (", trimws(mod), ")") else name
    } else {
        n <- min(length(parts), 3)
        paste(gsub("[_-]", " ", tail(parts, n)), collapse = " / ")
    }
}

# Targeted config discovery exploiting the known 3-level structure:
#   <root>/<experiment>/<crispr-screen[-modifier]>/config.yaml
# Uses only two shallow (non-recursive) list.dirs calls per root, then a
# direct file.exists check — no deep traversal, no Windows MAX_PATH issues.
safe_ld <- function(path) tryCatch(
    suppressWarnings(list.dirs(path, recursive = FALSE, full.names = TRUE)),
    error = function(e) character(0)
)

discover_configs <- function(search_dirs) {
    configs <- character(0)
    for (root in search_dirs) {
        if (!isTRUE(dir.exists(root))) next
        for (exp_dir in safe_ld(root)) {
            cs_dirs <- safe_ld(exp_dir)
            cs_dirs <- cs_dirs[grepl("^crispr-screen", basename(cs_dirs), ignore.case = TRUE)]
            for (cs_dir in cs_dirs) {
                cfg <- file.path(cs_dir, "config.yaml")
                if (file.exists(cfg)) configs <- c(configs, cfg)
            }
        }
    }
    unique(configs)
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
    tags$head(
        tags$style(HTML(paste0(
            "
            .sidebar { font-size: 90%; }
            .sidebar .form-group { margin-bottom: 8px; }
            .sidebar h4 { margin-top: 8px; margin-bottom: 8px; font-size: 1.1em; }
            .sidebar label { margin-bottom: 3px; }
            .sidebar .control-label { font-size: 0.95em; }
            .section-header {
                cursor: pointer; padding: 8px 10px; margin: 5px 0;
                background-color: #f8f9fa; border: 1px solid #dee2e6;
                border-radius: 4px; user-select: none; display: block;
                text-decoration: none; color: #212529; font-weight: 500;
            }
            .section-header:hover { background-color: #e9ecef; text-decoration: none; color: #212529; }
            .section-header .arrow { float: right; transition: transform 0.2s; }
            .section-content { padding: 10px 5px; border-left: 2px solid #dee2e6; margin-left: 5px; margin-bottom: 10px; }
            .counts-plot-header { padding: 3px 8px !important; margin: 2px 0 !important; font-size: 0.9em; }
            .counts-plot-content { padding: 3px 4px 0 4px !important; margin-left: 0 !important; margin-bottom: 2px !important; border-left: none !important; }
            .counts-plot-content .form-group { margin-bottom: 3px !important; }
            .counts-plot-content .control-label { font-size: 0.85em !important; margin-bottom: 1px !important; }
            .counts-move-btn { cursor:pointer; color:#bbb; font-size:0.9em; line-height:1; padding:0 3px; user-select:none; }
            .counts-move-btn:hover { color:#555; }
            .sortable-enabled .counts-move-btn { display:none !important; }
            .sortable-enabled .counts-plot-header { cursor:grab !important; }
            .sortable-enabled .counts-plot-header:active { cursor:grabbing !important; }
            .palette-dropdown .selectize-control { margin-bottom: 0; }
            .palette-dropdown .selectize-control.single .selectize-input { min-height: 32px; height: 32px; padding-top: 4px; padding-bottom: 4px; }
            .palette-dropdown .selectize-control.single .selectize-input.input-active { min-height: 32px; height: 32px; }
            .palette-dropdown .selectize-control.single { display: inline-flex; align-items: center; }
            .palette-dropdown .selectize-input > input { display: none; }
            .palette-dropdown .selectize-control.single .selectize-input { width: ", PALETTE_WIDTH, " !important; max-width: ", PALETTE_WIDTH, " !important; min-width: ", PALETTE_WIDTH, " !important; }
            .config-panel { background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 4px; padding: 10px 15px; margin-bottom: 12px; }
            .qc-image-container { text-align: center; }
            .qc-image-container img { max-width: 100%; height: auto; cursor: zoom-in; }
            .modal-backdrop { z-index: 10000; }
            .modal { z-index: 10001; }
            .modal-header .close { font-size: 2.5em; font-weight: 700; opacity: 0.7; padding: 5px 12px; }
            .modal-header .close:hover { opacity: 1; }
            .modal-lg { max-width: min(95vw, 1400px) !important; }
            .qc-overlay-scroll { overflow: auto; max-height: 80vh; }
            .qc-overlay-img { display: block; width: 100%; height: auto; }
            .qc-zoom-btn.active { background-color: #337ab7 !important; color: #fff !important; border-color: #2e6da4 !important; }
        "))),
        tags$script(HTML(paste0("
            $(document).on('shiny:connected', function() {
                $(document).on('click', '.section-header', function() {
                    var targetId = $(this).data('target');
                    var content = $('#' + targetId);
                    var arrow = $(this).find('.arrow');
                    content.slideToggle(200);
                    if (arrow.text() === '▼') arrow.text('▶'); else arrow.text('▼');
                });
                $('#config_path_input').on('keypress', function(e) {
                    if (e.which === 13) { e.preventDefault(); $('#load_config_path').click(); }
                });
            });
            $(document).on('click', '#copy_config_path', function() {
                var val = $('#config_path_text').val();
                if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(val);
                else { var $temp = $('<textarea>'); $('body').append($temp); $temp.val(val).select(); document.execCommand('copy'); $temp.remove(); }
            });
            $(document).ready(function() {
                setTimeout(function() {
                    var picker2 = $('#fdr_color_2');
                    var expectedColor = '", DEFAULT_FDR_COLOR_2, "';
                    if (picker2.length && picker2.val() === '#000000') {
                        picker2.val(expectedColor); picker2.css('background-color', expectedColor); picker2.css('color', '#000');
                    }
                }, 500);
            });
        "))),
        tags$script(src = "Sortable.min.js"),
        tags$script(HTML("
            var _countsSortable = null;

            function initCountsSortable() {
                var el = document.getElementById('counts_plot_configs');
                if (!el) return;
                if (typeof Sortable === 'undefined') {
                    Shiny.setInputValue('counts_sortable_ok', false, {priority: 'event'});
                    return;
                }
                if (_countsSortable) { try { _countsSortable.destroy(); } catch(e) {} }
                try {
                    _countsSortable = new Sortable(el, {
                        handle: '.counts-plot-header',
                        animation: 150,
                        onStart: function() { document.body.style.cursor = 'grabbing'; },
                        onEnd: function() {
                            document.body.style.cursor = '';
                            var items = el.querySelectorAll('[data-plot-id]');
                            var order = Array.from(items).map(function(d) {
                                return d.getAttribute('data-plot-id');
                            });
                            Shiny.setInputValue('counts_plot_order', order, {priority: 'event'});
                        }
                    });
                    el.classList.add('sortable-enabled');
                    Shiny.setInputValue('counts_sortable_ok', true, {priority: 'event'});
                } catch(e) {
                    Shiny.setInputValue('counts_sortable_ok', false, {priority: 'event'});
                }
            }

            // Re-initialize after each renderUI update of the plot config panel
            $(document).on('shiny:value', function(e) {
                if (e.name === 'counts_plot_configs') {
                    setTimeout(initCountsSortable, 60);
                }
            });

            // Initial availability check: wait 2s for CDN script to load, then report
            $(document).on('shiny:connected', function() {
                setTimeout(function() {
                    if (typeof Sortable === 'undefined') {
                        Shiny.setInputValue('counts_sortable_ok', false, {priority: 'event'});
                    }
                }, 2000);
            });

            // Rename modal: auto-focus + select text when modal finishes opening
            $(document).on('shown.bs.modal', function() {
                var inp = document.getElementById('counts_rename_input');
                if (inp) { inp.focus(); inp.select(); }
            });
            // Rename modal: Enter confirms, Escape cancels
            $(document).on('keydown', '#counts_rename_input', function(e) {
                if (e.which === 13) { e.preventDefault(); $('#counts_rename_confirm').click(); }
                if (e.which === 27) { e.preventDefault(); $(this).closest('.modal').modal('hide'); }
            });
        "))
    ),

    titlePanel("CRISPR Screen Explorer"),

    tags$div(
        style = "position: fixed; top: 15px; right: 20px; z-index: 9999; display: flex; align-items: center; gap: 10px;",
        uiOutput("loaded_experiment_label"),
        actionButton("reload_data", "Reload Data", class = "btn-warning btn-sm",
                     title = "Re-scan output files using the currently loaded config"),
        tags$span("v1.1", style = "font-size: 0.75em; color: #aaa;")
    ),

    # ── Main tabs ──────────────────────────────────────────────────────────────
    tabsetPanel(
        id = "main_tabs",
        selected = "Configuration",

        # ── Configuration tab ─────────────────────────────────────────────────
        tabPanel("Configuration",
            br(),
            div(class = "config-panel",
                fluidRow(
                    column(10,
                        tags$label("Select Experiment:", `for` = "experiment_select",
                                   style = "font-weight:600; margin-bottom:4px; display:block;"),
                        selectInput("experiment_select", NULL,
                                    choices = c("Scanning..." = ""),
                                    width = "100%")
                    ),
                    column(2,
                        tags$label("\u00a0", style = "display:block; margin-bottom:4px;"),
                        actionButton("rescan_experiments", "Rescan",
                                     class = "btn-default btn-sm", style = "width:100%;",
                                     title = "Search for config.yaml files on disk")
                    )
                ),
                tags$hr(style = "margin: 10px 0;"),
                tags$p(style = "font-size:0.85em; color:#666; margin-bottom:8px;",
                       "Or specify a config manually:"),
                fluidRow(
                    column(2, shinyFiles::shinyFilesButton("config_select", "Browse...",
                                                           title = "Select a config.yaml", multiple = FALSE)),
                    column(2, actionButton("reload_default", "Reload Default", class = "btn-default btn-sm")),
                    column(5, textInput("config_path_input", NULL, value = "",
                                        placeholder = "Paste path to config.yaml or use Browse")),
                    column(2, actionButton("load_config_path", "Load Config",
                                           class = "btn-primary btn-sm", style = "width:100%;")),
                    column(1, checkboxInput("show_debug", "Debug", value = TRUE))
                ),
                uiOutput("config_error_message")
            ),
            uiOutput("config_loaded_summary")
        ), # end Configuration tab

        # ── Comparisons tab ───────────────────────────────────────────────────
        tabPanel("Comparisons",
            br(),
            sidebarLayout(
                sidebarPanel(
                    width = 3, class = "sidebar",

                    selectInput("comparison", "Select Comparison:", choices = NULL),
                    uiOutput("global_sig_warning"),

                    radioButtons("selection_type", "Selection Type:",
                                choices = c("Negative" = "neg", "Positive" = "pos"),
                                selected = "neg", inline = TRUE),
                    radioButtons("plot_metric", "Plot metric",
                                choices = c("LFC" = "lfc", "Score" = "score"),
                                selected = "lfc", inline = TRUE),
                    radioButtons("rank_method", "Rank genes by",
                                choices = c("LFC" = "lfc", "RRA" = "rra"),
                                selected = "lfc", inline = TRUE),
                    checkboxInput("reverse_rank_order", "Reverse rank order", value = FALSE),
                    hr(style = "margin: 10px 0;"),

                    if (PLOTLY_AVAILABLE) checkboxInput("use_interactive_plots", "Interactive Plots", value = TRUE)
                    else tags$p(style = "font-size:0.9em;color:#999;margin:8px 0;", "Interactive plots not available."),

                    tags$a(class = "section-header", `data-target` = "highlight_content",
                           "Highlight Groups", tags$span(class = "arrow", "▼")),
                    tags$div(id = "highlight_content", class = "section-content", style = "display:block;",
                        checkboxInput("show_labels", "Label points", value = TRUE),
                        conditionalPanel(
                            condition = "input.show_labels",
                            if (PLOTLY_AVAILABLE)
                                tags$p(style = "font-size:0.8em;color:#888;margin:2px 0 4px 0;",
                                       "\u2139 Labels only appear in static plot mode."),
                            sliderInput("label_scale", "Label scale:", min = 0.5, max = 2.5, value = 1, step = 0.1)
                        ),
                        hr(style = "margin:6px 0;"),
                        uiOutput("gene_set_checkboxes")
                    ),

                    tags$a(class = "section-header", `data-target` = "significant_content",
                           "Top Genes", tags$span(class = "arrow", "▶")),
                    tags$div(id = "significant_content", class = "section-content", style = "display:none;",
                        checkboxInput("highlight_significant", "Label top genes", value = FALSE),
                        conditionalPanel(
                            condition = "input.highlight_significant && input.use_interactive_plots",
                            tags$p(style = "font-size:0.8em; color:#888; margin: 2px 0 6px 0;",
                                   "\u2139 Gene labels only appear in static plot mode.")
                        ),
                        numericInput("sig_top_n", "Max labeled genes:", value = 10, min = 1, max = 1000, step = 1),
                        checkboxInput("use_plot_thresholds", "Customize thresholds & colors", value = FALSE),
                        uiOutput("sig_warning"),
                        conditionalPanel(
                            condition = "input.use_plot_thresholds",
                            tags$p(style = "font-size:0.8em; color:#666; margin: 2px 0 4px 0;",
                                   "Set independent thresholds and colors for gene labeling:"),
                            fluidRow(
                                column(6,
                                    checkboxInput("sig_use_thresh1", "Enable threshold 1", value = TRUE),
                                    textInput("sig_fdr_threshold_1", "FDR <", value = "0.05", width = "80px"),
                                    div(class = "palette-dropdown",
                                        selectInput("sig_palette_1", "", choices = c("Basic"="limited","Web-safe"="websafe","Full"="square"), selected = "limited", width = PALETTE_WIDTH)),
                                    uiOutput("sig_color_picker_1")),
                                column(6,
                                    checkboxInput("sig_use_thresh2", "Enable threshold 2", value = FALSE),
                                    textInput("sig_fdr_threshold_2", "FDR <", value = "0.1", width = "80px"),
                                    div(class = "palette-dropdown",
                                        selectInput("sig_palette_2", "", choices = c("Basic"="limited","Web-safe"="websafe","Full"="square"), selected = "limited", width = PALETTE_WIDTH)),
                                    uiOutput("sig_color_picker_2"))
                            )
                        )
                    ),

                    tags$a(class = "section-header", `data-target` = "custom_content",
                           "Custom Gene List", tags$span(class = "arrow", "▶")),
                    tags$div(id = "custom_content", class = "section-content", style = "display:none;",
                        fileInput("custom_gene_file", "Upload Gene List File", accept = c(".txt",".csv",".tsv")),
                        textInput("custom_list_name", "List Name:", value = "Custom"),
                        fluidRow(
                            column(5, selectInput("custom_palette_type", "Color palette:",
                                                  choices = c("Basic"="limited","Web-safe"="websafe","Full"="square"),
                                                  selected = "limited", width = "100%")),
                            column(7, uiOutput("custom_color_picker"))
                        ),
                        conditionalPanel(
                            condition = "input.custom_palette_type == 'websafe'",
                            selectInput("websafe_order", "Web-Safe Ordering:",
                                       choices = c("RGB grid (R fast)"="rgbgrid","HSV (by hue)"="hsv"), selected = "rgbgrid")),
                        textAreaInput("custom_genes", "Gene List:", placeholder = "GENE1 GENE2\nGENE3, GENE4", rows = 3),
                        fluidRow(
                            column(4, actionButton("apply_custom", "Apply", class = "btn-primary btn-sm", style = "width:100%;")),
                            column(4, actionButton("clear_custom", "Clear Custom", class = "btn-default btn-sm", style = "width:100%;")),
                            column(4, checkboxInput("show_custom", "Show on plot", value = TRUE))
                        ),
                        uiOutput("custom_warning")
                    ),

                    tags$a(class = "section-header", `data-target` = "plot_content",
                           "Plot Options", tags$span(class = "arrow", "▶")),
                    tags$div(id = "plot_content", class = "section-content", style = "display:none;",
                        tags$style(HTML(paste0(
                            ".compact-row .shiny-input-container{margin-bottom:0;}\n",
                            ".compact-row input:not([type=\"checkbox\"]), .compact-row select {height: 32px;}\n",
                            ".compact-row .checkbox {display: inline-block; width: auto; margin: 0;}\n",
                            ".compact-row .checkbox label {display: inline-flex; align-items: center; gap: 6px; margin: 0; white-space: nowrap;}\n",
                            ".compact-row .checkbox input[type=\"checkbox\"] {position: static; margin: 0;}\n"
                        ))),
                        textInput("volcano_title", "Volcano Title:", value = "", placeholder = "Auto"),
                        textInput("ranked_title", "Ranked Title:", value = "", placeholder = "Auto"),
                        sliderInput("point_size", "Base Point Size:", min = 0.5, max = 4, value = 2, step = 0.25),
                        fluidRow(column(12,
                            tags$div(tags$strong("FDR Threshold 1:")),
                            tags$div(class = "compact-row", style = "display:flex;gap:8px;align-items:center;flex-wrap:wrap;",
                                textInput("fdr_threshold_1", label = NULL, value = "0.05", placeholder = "0.05", width = "60px"),
                                uiOutput("fdr_color_picker_1"),
                                div(class = "palette-dropdown", selectInput("fdr_palette_1", label = NULL,
                                    choices = c("Basic"="limited","Web-safe"="websafe","Full"="square"), selected = "limited", width = PALETTE_WIDTH)),
                                checkboxInput("show_fdr_line_1", label = "Show line", value = TRUE))
                        )),
                        fluidRow(column(12,
                            tags$div(tags$strong("FDR Threshold 2:")),
                            tags$div(class = "compact-row", style = "display:flex;gap:8px;align-items:center;flex-wrap:wrap;",
                                textInput("fdr_threshold_2", label = NULL, value = "0.1", placeholder = "0.1", width = "60px"),
                                uiOutput("fdr_color_picker_2"),
                                div(class = "palette-dropdown", selectInput("fdr_palette_2", label = NULL,
                                    choices = c("Basic"="limited","Web-safe"="websafe","Full"="square"), selected = "limited", width = PALETTE_WIDTH)),
                                checkboxInput("show_fdr_line_2", label = "Show line", value = TRUE))
                        )),
                        checkboxInput("scale_point_size", "Scale point size by FDR", value = TRUE),
                        checkboxInput("color_by_fdr", "Color points by FDR", value = TRUE),
                        sliderInput("text_scale", "Text Scale:", min = 0.5, max = 2.5, value = 1, step = 0.1),
                        sliderInput("title_scale", "Title Scale:", min = 0.5, max = 4, value = 1, step = 0.1)
                    ),

                    tags$a(class = "section-header", `data-target` = "download_content",
                           "Download Options", tags$span(class = "arrow", "▶")),
                    tags$div(id = "download_content", class = "section-content", style = "display:none;",
                        fluidRow(
                            column(6, selectInput("download_format", "Format:", choices = c("PNG"="png","PDF"="pdf"), selected = "png")),
                            column(6, selectInput("download_dpi", "DPI:", choices = c("72"=72,"150"=150,"300"=300,"600"=600,"1200"=1200), selected = "300"))
                        ),
                        fluidRow(
                            column(6, sliderInput("download_width", "Width (in):", min = 4, max = 20, value = 10, step = 0.5)),
                            column(6, sliderInput("download_height", "Height (in):", min = 4, max = 20, value = 8, step = 0.5))
                        ),
                        downloadButton("download_volcano", "Download Volcano", style = "width:100%;margin-bottom:5px;"),
                        downloadButton("download_ranked", "Download Ranked", style = "width:100%;")
                    )
                ),

                mainPanel(width = 9,
                    fluidRow(column(12,
                        h3("Ranked Gene Plot"),
                        if (PLOTLY_AVAILABLE) conditionalPanel(condition = "input.use_interactive_plots",
                            plotly::plotlyOutput("ranked_plot", height = "500px")),
                        if (PLOTLY_AVAILABLE) conditionalPanel(condition = "!input.use_interactive_plots",
                            plotOutput("ranked_plot_static", height = "500px"))
                        else plotOutput("ranked_plot", height = "500px")
                    )),
                    hr(),
                    fluidRow(column(12,
                        h3("Volcano Plot"),
                        if (PLOTLY_AVAILABLE) conditionalPanel(condition = "input.use_interactive_plots",
                            plotly::plotlyOutput("volcano_plot", height = "500px")),
                        if (PLOTLY_AVAILABLE) conditionalPanel(condition = "!input.use_interactive_plots",
                            plotOutput("volcano_plot_static", height = "500px"))
                        else plotOutput("volcano_plot", height = "500px")
                    ))
                )
            )
        ), # end Comparisons tab

        # ── QC tab ────────────────────────────────────────────────────────────
        tabPanel("QC",
            br(),
            uiOutput("qc_no_config_msg"),

            # Row 1: MAGeCK Stats (left) | Count Distributions (right)
            fluidRow(
                column(6,
                    wellPanel(
                        h4("MAGeCK Stats", style = "margin-top:0;"),
                        uiOutput("qc_stats_status"),
                        div(style = "overflow-x: auto;",
                            if (DT_AVAILABLE) DT::dataTableOutput("qc_countsummary_dt")
                            else tableOutput("qc_countsummary_table")
                        )
                    )
                ),
                column(6,
                    wellPanel(
                        h4("Count Distributions", style = "margin-top:0;"),
                        fluidRow(
                            column(7, selectInput("qc_hist_sample", "Sample:", choices = NULL)),
                            column(5, radioButtons("qc_hist_type", "Type:",
                                choices = c("Count"="Count","Log Count"="Log Count",
                                            "CPM"="CPM","Log CPM"="Log CPM"),
                                selected = "Log CPM"))
                        ),
                        uiOutput("qc_hist_status"),
                        div(class = "qc-image-container",
                            onclick = "Shiny.setInputValue('qc_image_click', {type:'hist', ts:Date.now()}, {priority:'event'})",
                            imageOutput("qc_hist_image", width = "100%", height = "auto"))
                    )
                )
            ),

            # Row 2: PCA Results (left) | Reorientation Barplots (right)
            fluidRow(
                column(6,
                    wellPanel(
                        h4("PCA Results", style = "margin-top:0;"),
                        selectInput("qc_pca_file", "Select plot:", choices = NULL),
                        uiOutput("qc_pca_status"),
                        div(class = "qc-image-container",
                            onclick = "Shiny.setInputValue('qc_image_click', {type:'pca', ts:Date.now()}, {priority:'event'})",
                            imageOutput("qc_pca_image", width = "100%", height = "auto"))
                    )
                ),
                column(6,
                    wellPanel(
                        h4("Reorientation Barplots", style = "margin-top:0;"),
                        selectInput("qc_barplot_sample", "Select sample:", choices = NULL),
                        uiOutput("qc_barplot_status"),
                        div(class = "qc-image-container",
                            onclick = "Shiny.setInputValue('qc_image_click', {type:'barplot', ts:Date.now()}, {priority:'event'})",
                            imageOutput("qc_barplot_image", width = "100%", height = "auto"))
                    )
                )
            )
        ) # end QC tab

        # ── Counts tab ────────────────────────────────────────────────────────
        ,tabPanel("Counts",
            br(),
            uiOutput("counts_no_config_msg"),
            conditionalPanel(
                condition = "output.counts_has_data",
                sidebarLayout(
                    sidebarPanel(width = 3,
                        selectInput("counts_matrix_type", "Matrix:",
                                    choices = c("Counts"="count", "CPM"="cpm"), selected = "count"),
                        uiOutput("counts_matrix_status"),
                        if (PLOTLY_AVAILABLE)
                            checkboxInput("use_interactive_counts", "Interactive plots", value = TRUE)
                        else
                            tags$p(style="font-size:0.9em;color:#999;margin:8px 0;", "Interactive plots not available."),
                        hr(),
                        tags$a(class = "section-header", `data-target` = "counts_highlight_content",
                               "Highlight Groups", tags$span(class = "arrow", "▼")),
                        tags$div(id = "counts_highlight_content", class = "section-content", style = "display:block;",
                            uiOutput("counts_gene_set_checkboxes")
                        ),
                        hr(),
                        fluidRow(
                            column(6, numericInput("counts_ncols",   "Columns:",      value=2,   min=1,   max=8,   step=1)),
                            column(6, numericInput("counts_plot_px", "Height (px):",  value=350, min=150, max=800, step=50))
                        ),
                        fluidRow(
                            column(6, numericInput("counts_plot_w",  "Width (px):",   value=350, min=150, max=800, step=50)),
                            column(6, div(style="margin-top:24px;",
                                checkboxInput("counts_aspect_ratio1", "Square panel (static)", value=TRUE)))
                        ),
                        uiOutput("counts_sortable_status"),
                        uiOutput("counts_plot_configs"),
                        fluidRow(style="margin:0;margin-top:4px;margin-bottom:4px;",
                            column(7, style="padding-right:2px;",
                                actionButton("counts_add_plot", "+ Add Plot",
                                             class="btn-sm btn-default", style="width:100%;")),
                            column(5, style="padding-left:2px;",
                                actionButton("counts_reset_names", "Reset names",
                                             class="btn-sm btn-default", style="width:100%;"))
                        ),
                        selectInput("counts_title_fmt", "Plot title:",
                            choices = c(
                                "X vs Y (plot name)" = "full",
                                "X vs Y"             = "xy",
                                "Plot name only"     = "name"
                            ), selected = "full"),
                        selectInput("counts_level", "Level:",
                            choices = c("sgRNA" = "sgrna", "Gene" = "gene"),
                            selected = "sgrna"),
                        conditionalPanel(condition = "input.counts_level == 'gene'",
                            selectInput("counts_gene_summary", "Summarize by:",
                                choices = c("Mean"="mean","Median"="median","Max"="max"),
                                selected = "mean"),
                            fluidRow(
                                column(6, numericInput("counts_gene_min", "Exclude if count <",
                                                       value=0, min=0, step=1)),
                                column(6, selectInput("counts_gene_filter", "in sample(s):",
                                    choices = c("—"="none",
                                                "Either"="either",
                                                "Both"="both"),
                                    selected = "none"))
                            )
                        ),
                        fluidRow(
                            column(6, selectInput("counts_transform", "Transform:",
                                choices = c("None"="none",
                                            "log\u2081\u2080"="log10",
                                            "log\u2082"="log2"),
                                selected = "none")),
                            column(6, conditionalPanel(
                                condition = "input.counts_transform != 'none'",
                                numericInput("counts_pseudocount", "Pseudocount:",
                                             value = 1, min = 0, step = 0.1)))
                        ),
                        uiOutput("counts_log_warning"),
                        hr(),
                        h5("Identity Line (y = x)", style="margin-bottom:6px;"),
                        checkboxInput("counts_show_identity", "Show identity line", value = TRUE),
                        conditionalPanel(condition = "input.counts_show_identity",
                            fluidRow(
                                column(5, div(class="palette-dropdown",
                                    selectInput("counts_id_palette", "",
                                        choices=c("Basic"="limited","Web-safe"="websafe","Full"="square"),
                                        selected="limited", width="100%"))),
                                column(7, uiOutput("counts_id_color_picker"))
                            ),
                            selectInput("counts_id_lty", "Line style:",
                                choices=c("Dashed"="dashed","Solid"="solid","Dotted"="dotted",
                                          "Long dash"="longdash","Dot-dash"="dotdash"),
                                selected="dashed"),
                            sliderInput("counts_id_lwd", "Line width:", min=0.3, max=4, value=1, step=0.1)
                        ),
                        hr(),
                        h5("Scale & Appearance", style="margin-bottom:6px;"),
                        sliderInput("counts_point_size",  "Point size:",  min=0.1, max=3,   value=0.6, step=0.1),
                        sliderInput("counts_point_alpha", "Opacity:",     min=0.05, max=1,  value=0.35, step=0.05),
                        if (GGEXTRA_AVAILABLE)
                            checkboxInput("counts_show_hist", "Marginal histograms (static only)", value = FALSE)
                        else
                            tags$p(style="font-size:0.85em;color:#999;margin:4px 0;", "Install ggExtra for marginal histograms."),
                        hr(),
                        h5("Download", style="margin-bottom:6px;"),
                        fluidRow(
                            column(6, selectInput("counts_dl_fmt", "Format:",
                                choices=c("PNG"="png","PDF"="pdf"), selected="png")),
                            column(6, selectInput("counts_dl_dpi", "DPI:",
                                choices=c("150"=150,"300"=300,"600"=600), selected="300"))
                        ),
                        fluidRow(
                            column(6, sliderInput("counts_dl_w", "Width (in):",  min=2, max=30, value=7, step=0.5)),
                            column(6, sliderInput("counts_dl_h", "Height (in):", min=2, max=30, value=7, step=0.5))
                        ),
                        downloadButton("download_counts_scatter", "Download Plot",
                                       style="width:100%;", class="btn-primary btn-sm"),
                        hr(),
                        h6("Download Data", style="margin-bottom:4px;"),
                        fluidRow(
                            column(6, downloadButton("download_counts_tsv", "TSV",
                                style="width:100%;", class="btn-default btn-sm")),
                            column(6,
                                if (OPENXLSX2_AVAILABLE || OPENXLSX_AVAILABLE)
                                    tagList(
                                        downloadButton("download_counts_xlsx", "Excel",
                                            style="width:100%;", class="btn-default btn-sm"),
                                        tags$small(
                                            if (OPENXLSX2_AVAILABLE) "openxlsx2" else "openxlsx",
                                            style="color:#888; display:block; text-align:center; margin-top:2px;")
                                    )
                                else
                                    tags$small("Excel unavailable (install openxlsx2)",
                                        style="color:#888; display:block; margin-top:4px;")
                            )
                        )
                    ),
                    mainPanel(width = 9,
                        if (PLOTLY_AVAILABLE) conditionalPanel(
                            condition = "input.use_interactive_counts",
                            uiOutput("counts_plotly_ui")),
                        if (PLOTLY_AVAILABLE) conditionalPanel(
                            condition = "!input.use_interactive_counts",
                            uiOutput("counts_static_ui"))
                        else uiOutput("counts_static_ui")
                    )
                )
            )
        ) # end Counts tab
    ) # end tabsetPanel
) # end fluidPage


# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
    cat("\n========== SERVER FUNCTION STARTING ==========\n")

    steps <- c("00","33","66","99","CC","FF")
    websafe_colors_all <- as.vector(outer(steps, steps, FUN = function(r, g) paste0(r, g)))
    websafe_colors_all <- unlist(lapply(steps, function(b) paste0("#", paste0(websafe_colors_all, b))))
    websafe_colors_rgbgrid <- unlist(lapply(steps, function(b)
        unlist(lapply(steps, function(g) unlist(lapply(steps, function(r) paste0("#", r, g, b)))))))
    hsv_vals <- grDevices::rgb2hsv(col2rgb(websafe_colors_rgbgrid), maxColorValue = 255)
    ord_hsv <- order(hsv_vals["h",], hsv_vals["s",], hsv_vals["v",])
    websafe_colors_hsv <- websafe_colors_rgbgrid[ord_hsv]
    websafe_colors <- websafe_colors_rgbgrid

    log_debug <- function(msg) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n"))

    rv <- reactiveValues(
        config = NULL, config_dir = NULL, gene_data = NULL, comparisons = NULL,
        result_files = NULL, custom_gene_list = NULL, gene_sets_enabled = list(),
        neg_ctrl_override = NULL,
        warning_message = NULL, sig_warning = NULL, config_error = NULL, config_warnings = NULL,
        has_highlight_sets = FALSE,
        sig_color_1 = DEFAULT_FDR_COLOR_1, sig_color_2 = DEFAULT_FDR_COLOR_2,
        fdr_color_1 = DEFAULT_FDR_COLOR_1, fdr_color_2 = DEFAULT_FDR_COLOR_2,
        # QC state
        output_dir = NULL, experiment_name = NULL,
        qc_countsummary = NULL, qc_histogram_meta = NULL,
        qc_pca_files = NULL, qc_barplot_files = NULL,
        # Counts matrix state
        count_matrix_file = NULL, cpm_matrix_file = NULL,
        counts_id_color = "#FF0000", counts_hl_color = "#E41A1C",
        counts_plot_ids = "p1", counts_plot_id_counter = 1L,
        counts_plot_names = list(p1 = "Plot 1"), counts_renaming_pid = NULL,
        # Experiment discovery
        discovered_configs = NULL
    )

    clear_state <- function() {
        rv$config <- NULL; rv$gene_data <- NULL; rv$result_files <- NULL; rv$comparisons <- NULL
        updateSelectInput(session, "comparison", choices = NULL)
    }

    normalize_comparisons <- function(comparisons) {
        if (is.null(comparisons) || length(comparisons) == 0) return(NULL)
        normalized <- list()
        for (item in comparisons) {
            if (is.character(item)) normalized <- c(normalized, list(item))
            else if (is.list(item) && "comparison" %in% names(item)) normalized <- c(normalized, list(item))
            else if (is.list(item) && length(item) == 1 && is.character(item[[1]])) normalized <- c(normalized, list(item[[1]]))
        }
        if (length(normalized) == 0) return(NULL)
        return(normalized)
    }

    load_config_from_path <- function(config_path) {
        cat("LOAD_CONFIG: Starting with path:", config_path, "\n")
        rv$config_error <- NULL
        rv$config_warnings <- NULL

        # ── Hard failures: cannot proceed at all ──────────────────────────────
        if (is.null(config_path) || !nzchar(config_path)) { rv$config_error <- "No config file selected."; clear_state(); return(FALSE) }
        if (!file.exists(config_path)) { rv$config_error <- paste0("Config file not found: ", config_path); clear_state(); return(FALSE) }

        rv$config_dir <- dirname(config_path)
        tryCatch({ rv$config <- parse_config(config_path); cat("LOAD_CONFIG: Config parsed\n") },
                 error = function(e) { rv$config_error <- paste0("Error parsing config: ", e$message); rv$config <- NULL })
        if (is.null(rv$config)) { if (is.null(rv$config_error)) rv$config_error <- "Failed to parse config."; clear_state(); return(FALSE) }

        output_dir <- rv$config$output_dir %||% "output"
        experiment_name <- rv$config$experiment_name %||% "experiment"
        output_dir_relative <- file.path(rv$config_dir, output_dir)

        if (file.exists(output_dir_relative)) {
            output_dir <- output_dir_relative
        } else if (!file.exists(output_dir)) {
            rv$config_error <- paste0("Output directory not found. Tried: ", output_dir_relative, " and ", output_dir)
            clear_state(); return(FALSE)
        }

        # Output dir and experiment name are valid — store now so QC always loads
        rv$output_dir     <- output_dir
        rv$experiment_name <- experiment_name

        # ── Soft failures: warn but don't block QC ────────────────────────────
        warns <- character(0)

        comparisons_raw <- rv$config$comparisons
        if (is.null(comparisons_raw)) {
            warns <- c(warns, "Config missing 'comparisons' field — Comparisons tab unavailable.")
            rv$result_files <- NULL
            updateSelectInput(session, "comparison", choices = NULL)
        } else {
            comparisons <- normalize_comparisons(comparisons_raw)
            if (is.null(comparisons) || length(comparisons) == 0) {
                warns <- c(warns, "Config 'comparisons' is empty or invalid — Comparisons tab unavailable.")
                rv$result_files <- NULL
                updateSelectInput(session, "comparison", choices = NULL)
            } else {
                tryCatch({
                    rv$result_files <- find_rra_results(output_dir, experiment_name, comparisons)
                    if (is.null(rv$result_files) || length(rv$result_files) == 0) {
                        warns <- c(warns, paste0("No RRA result files found in: ", output_dir, " — Comparisons tab unavailable."))
                        updateSelectInput(session, "comparison", choices = NULL)
                    } else {
                        updateSelectInput(session, "comparison", choices = names(rv$result_files))
                    }
                }, error = function(e) {
                    warns <<- c(warns, paste0("Error finding result files: ", e$message))
                    updateSelectInput(session, "comparison", choices = NULL)
                })
            }
        }

        if (length(warns) > 0) rv$config_warnings <- warns

        updateTextInput(session, "config_path_input", value = config_path)
        cat("LOAD_CONFIG: SUCCESS\n")
        return(TRUE)
    }

    # File browser setup
    vols <- shinyFiles::getVolumes()()
    roots <- vols
    default_root_path <- Sys.getenv("APP_ROOT_PATH", unset = "/data")
    if (file.exists(default_root_path)) roots <- c(roots, CRISPR = default_root_path)
    roots <- c(roots, WD = getwd(), Home = normalizePath("~", winslash = "/", mustWork = FALSE))

    if ("CRISPR" %in% names(roots))
        shinyFiles::shinyFileChoose(input, id = "config_select", roots = roots, filetypes = c("yaml","yml"), defaultRoot = "CRISPR", defaultPath = ".")
    else if ("WD" %in% names(roots))
        shinyFiles::shinyFileChoose(input, id = "config_select", roots = roots, filetypes = c("yaml","yml"), defaultRoot = "WD", defaultPath = ".")
    else
        shinyFiles::shinyFileChoose(input, id = "config_select", roots = roots, filetypes = c("yaml","yml"))

    # ── Experiment discovery ──────────────────────────────────────────────────

    search_roots <- unique(c(default_root_path, dirname(getwd()), getwd()))

    scan_and_update_experiments <- function() {
        cat("SCAN: Searching for config files under:", paste(search_roots, collapse = ", "), "\n")
        paths <- discover_configs(search_roots)
        cat("SCAN: Found", length(paths), "config file(s)\n")
        if (length(paths) == 0) {
            rv$discovered_configs <- character(0)
            updateSelectInput(session, "experiment_select",
                              choices = c("No experiments found" = ""))
            return(invisible(NULL))
        }
        labels <- sapply(paths, config_path_to_label, USE.NAMES = FALSE)
        # Deduplicate labels by appending a counter if needed
        dupes <- which(duplicated(labels) | duplicated(labels, fromLast = TRUE))
        if (length(dupes) > 0) {
            for (lbl in unique(labels[dupes])) {
                idx <- which(labels == lbl)
                labels[idx] <- paste0(lbl, " [", seq_along(idx), "]")
            }
        }
        choices <- c("-- select an experiment --" = "", setNames(paths, labels))
        rv$discovered_configs <- setNames(paths, labels)
        updateSelectInput(session, "experiment_select", choices = choices, selected = "")
    }

    observeEvent(input$rescan_experiments, {
        scan_and_update_experiments()
    })

    observeEvent(input$experiment_select, {
        path <- input$experiment_select
        if (!is.null(path) && nzchar(path)) {
            updateTextInput(session, "config_path_input", value = path)
            load_config_from_path(path)
        }
    }, ignoreInit = TRUE)

    output$config_error_message <- renderUI({
        items <- list()
        if (!is.null(rv$config_error))
            items <- c(items, list(tags$div(
                style = "margin-top:8px;padding:8px;background-color:#f8d7da;color:#721c24;border:1px solid #f5c6cb;border-radius:4px;font-size:0.9em;",
                tags$strong("Error: "), rv$config_error)))
        if (!is.null(rv$config_warnings))
            items <- c(items, list(tags$div(
                style = "margin-top:8px;padding:8px;background-color:#fff3cd;color:#856404;border:1px solid #ffc107;border-radius:4px;font-size:0.9em;",
                tags$strong("Warning: "),
                if (length(rv$config_warnings) == 1) rv$config_warnings
                else tags$ul(style="margin:4px 0 0 0;padding-left:18px;",
                    lapply(rv$config_warnings, tags$li)))))
        if (length(items) > 0) do.call(tagList, items)
    })

    observeEvent(input$config_select, {
        paths <- shinyFiles::parseFilePaths(roots, input$config_select)
        if (nrow(paths) > 0) {
            cfg <- as.character(paths$datapath[1])
            updateTextInput(session, "config_path_input", value = cfg)
            load_config_from_path(cfg)
        }
    }, ignoreInit = TRUE)

    observeEvent(input$show_debug, { if (input$show_debug) cat("[DEBUG] Debugging on\n") })

    observeEvent(input$load_config_path, {
        if (!is.null(input$config_path_input) && nzchar(input$config_path_input)) load_config_from_path(input$config_path_input)
        else rv$config_error <- "Please enter a config path."
    })

    observeEvent(input$reload_default, {
        default_config <- file.path(default_root_path, "config.yaml")
        if (file.exists(default_config)) load_config_from_path(default_config)
        else if (file.exists("config.yaml")) load_config_from_path("config.yaml")
        else { rv$config_error <- "No default config found."; clear_state() }
    })

    observeEvent(input$reload_data, {
        scan_and_update_experiments()
        path <- input$config_path_input
        if (!is.null(path) && nzchar(path)) load_config_from_path(path)
        else rv$config_error <- "No config loaded. Use Browse or paste a path and click Load Config first."
    })

    output$loaded_experiment_label <- renderUI({
        req(rv$config, rv$output_dir)
        label <- config_path_to_label(input$config_path_input)
        tags$span(
            style = "font-size: 0.85em; color: #555; white-space: nowrap;",
            tags$strong("Loaded: "), label
        )
    })

    output$config_loaded_summary <- renderUI({
        if (!is.null(rv$config) && !is.null(rv$output_dir)) {
            n_comp <- length(rv$result_files)
            tags$div(style = "margin-top:10px;padding:8px;background-color:#d4edda;color:#155724;border:1px solid #c3e6cb;border-radius:4px;font-size:0.9em;",
                tags$strong("Config loaded. "),
                paste0("Experiment: ", rv$experiment_name %||% "unknown", " | ",
                       n_comp, " comparison(s) found | Output dir: ", rv$output_dir))
        }
    })

    isolate({
        scan_and_update_experiments()
        default_config <- file.path(default_root_path, "config.yaml")
        if (file.exists(default_config)) load_config_from_path(default_config)
        else if (file.exists("config.yaml")) load_config_from_path("config.yaml")
        # Auto-load errors are not user-initiated — clear them so the UI starts clean
        rv$config_error <- NULL
    })

    # ── QC file scanning ──────────────────────────────────────────────────────

    observeEvent(rv$output_dir, {
        req(rv$output_dir)
        out_dir <- rv$output_dir
        addResourcePath("qc_files", normalizePath(out_dir, mustWork = FALSE))
        cat("QC SCAN: Scanning", out_dir, "\n")

        qc_files <- find_qc_files(out_dir)

        # Countsummary
        rv$qc_countsummary <- load_countsummary_data(qc_files$countsummary)
        cat("QC SCAN: Found", length(qc_files$countsummary), "countsummary files\n")

        # Histograms
        if (length(qc_files$histograms) > 0) {
            hist_meta <- lapply(qc_files$histograms, parse_histogram_metadata)
            rv$qc_histogram_meta <- hist_meta
            samples <- sort(unique(sapply(hist_meta, function(x) x$sample)))
            updateSelectInput(session, "qc_hist_sample", choices = samples,
                              selected = if (length(samples) > 0) samples[1] else NULL)
            cat("QC SCAN: Found", length(hist_meta), "histogram files,", length(samples), "samples\n")
        } else {
            rv$qc_histogram_meta <- NULL
            updateSelectInput(session, "qc_hist_sample", choices = NULL)
        }

        # PCA plots
        if (length(qc_files$pca) > 0) {
            pca_labels <- gsub("_", " ", tools::file_path_sans_ext(basename(qc_files$pca)))
            pca_choices <- setNames(qc_files$pca, pca_labels)
            rv$qc_pca_files <- pca_choices
            updateSelectInput(session, "qc_pca_file", choices = pca_choices)
            cat("QC SCAN: Found", length(pca_choices), "PCA plot files\n")
        } else {
            rv$qc_pca_files <- NULL
            updateSelectInput(session, "qc_pca_file", choices = NULL)
        }

        # Barplots
        if (length(qc_files$barplots) > 0) {
            bar_labels <- sub("_barplot\\.png$", "", basename(qc_files$barplots))
            bar_choices <- setNames(qc_files$barplots, bar_labels)
            rv$qc_barplot_files <- bar_choices
            updateSelectInput(session, "qc_barplot_sample", choices = bar_choices)
            cat("QC SCAN: Found", length(bar_choices), "barplot files\n")
        } else {
            rv$qc_barplot_files <- NULL
            updateSelectInput(session, "qc_barplot_sample", choices = NULL)
        }

        # Count / CPM matrix files
        exp_name <- rv$experiment_name %||% "experiment"
        cfg <- rv$config %||% list()
        find_matrix_file <- function(cfg_key, default_suffix) {
            explicit <- cfg[[cfg_key]]
            if (!is.null(explicit) && nzchar(explicit)) {
                for (p in c(explicit, file.path(out_dir, explicit), file.path(rv$config_dir %||% ".", explicit)))
                    if (file.exists(p)) return(normalizePath(p))
            }
            default <- file.path(out_dir, paste0(exp_name, default_suffix))
            if (file.exists(default)) normalizePath(default) else NULL
        }
        rv$count_matrix_file <- find_matrix_file("count_matrix_file", "_sgrna_count_matrix.txt")
        rv$cpm_matrix_file   <- find_matrix_file("cpm_matrix_file",   "_sgrna_cpm_matrix.txt")
        cat("COUNTS: count matrix:", rv$count_matrix_file %||% "not found", "\n")
        cat("COUNTS: CPM matrix:  ", rv$cpm_matrix_file   %||% "not found", "\n")
    })

    # Build a browser-accessible URL for a file under the registered qc_files resource path.
    qc_file_url <- function(filepath) {
        base <- normalizePath(rv$output_dir, mustWork = FALSE)
        file <- normalizePath(filepath,      mustWork = FALSE)
        base_esc <- gsub("([][{()*+?.\\\\^$|])", "\\\\\\1", base)
        rel  <- sub(paste0("^", base_esc, "[/\\\\]?"), "", file)
        rel  <- gsub("\\\\", "/", rel)
        paste0("/qc_files/", paste(vapply(strsplit(rel, "/")[[1]], URLencode, character(1), repeated = TRUE), collapse = "/"))
    }

    # ── QC image click → full-resolution overlay modal ───────────────────────

    observeEvent(input$qc_image_click, {
        req(rv$output_dir)
        type <- input$qc_image_click$type

        filepath <- switch(type,
            hist = {
                req(rv$qc_histogram_meta, input$qc_hist_sample, input$qc_hist_type)
                m <- Filter(function(x) x$sample == input$qc_hist_sample &&
                                        x$type   == input$qc_hist_type, rv$qc_histogram_meta)
                if (length(m) > 0 && file.exists(m[[1]]$filepath)) m[[1]]$filepath else NULL
            },
            pca     = { req(input$qc_pca_file);     if (file.exists(input$qc_pca_file))     input$qc_pca_file     else NULL },
            barplot = { req(input$qc_barplot_sample); if (file.exists(input$qc_barplot_sample)) input$qc_barplot_sample else NULL }
        )
        req(filepath)

        img_url  <- qc_file_url(filepath)
        pdf_path <- paste0(tools::file_path_sans_ext(filepath), ".pdf")
        pdf_btn  <- if (file.exists(pdf_path))
            tags$a(href = qc_file_url(pdf_path), target = "_blank", rel = "noopener",
                   class = "btn btn-default btn-sm", style = "margin-right:8px;",
                   "\U1F4C4 Open PDF")
        else NULL

        dl_btn <- tags$a(
            href = img_url, download = basename(filepath),
            class = "btn btn-primary btn-sm", style = "margin-right:8px;",
            "\u2B07 Download"
        )
        zoom_btns <- tags$div(
            class = "btn-group btn-group-sm", style = "margin-right:12px;",
            tags$button(class="btn btn-default qc-zoom-btn active",
                onclick="qcSetZoom('fit',this)", title="Scale to fit window", "Fit"),
            tags$button(class="btn btn-default qc-zoom-btn",
                onclick="qcSetZoom(50,this)",  title="50% of native resolution",  "50%"),
            tags$button(class="btn btn-default qc-zoom-btn",
                onclick="qcSetZoom(100,this)", title="Native resolution (1 image pixel = 1 screen pixel)", "100%"),
        )
        showModal(modalDialog(
            title = basename(filepath),
            tags$div(class = "qc-overlay-scroll",
                tags$img(src = img_url, class = "qc-overlay-img")),
            tags$script(HTML("
                function qcSetZoom(pct, btn) {
                    var img = document.querySelector('.qc-overlay-img');
                    function apply() {
                        if (pct === 'fit') {
                            img.style.width = '100%';
                            img.style.maxWidth = '';
                        } else {
                            var nw = img.naturalWidth;
                            img.style.width = (nw > 0 ? Math.round(nw * pct / 100) : pct) + 'px';
                            img.style.maxWidth = 'none';
                        }
                        img.style.height = 'auto';
                    }
                    if (img.complete && img.naturalWidth > 0) apply();
                    else img.addEventListener('load', apply, {once: true});
                    document.querySelectorAll('.qc-zoom-btn').forEach(function(b) {
                        b.classList.remove('active');
                    });
                    btn.classList.add('active');
                }
            ")),
            footer = tagList(zoom_btns, dl_btn, pdf_btn, modalButton("Close")),
            size = "l", easyClose = TRUE
        ))
    })

    # ── QC renders ────────────────────────────────────────────────────────────

    output$qc_no_config_msg <- renderUI({
        if (is.null(rv$output_dir))
            tags$div(style = "padding:15px;background:#fff3cd;border:1px solid #ffc107;border-radius:4px;margin-bottom:10px;",
                     tags$strong("No experiment loaded."), " Load a config.yaml using the controls above.")
    })

    output$qc_stats_status <- renderUI({
        if (is.null(rv$qc_countsummary))
            tags$p(style = "color:#888;", "No countsummary.txt files found in the output directory.")
    })

    if (DT_AVAILABLE) {
        output$qc_countsummary_dt <- DT::renderDataTable({
            req(rv$qc_countsummary)
            df <- rv$qc_countsummary
            # Select key columns (drop NegSelQC columns)
            cn <- tolower(names(df))
            key_cols <- c(label="label", reads="reads", mapped="mapped",
                          percentage="percentage", zerocounts="zerocounts",
                          giniindex="giniindex")
            found <- sapply(key_cols, function(k) which(cn == k)[1])
            found <- found[!is.na(found)]
            if (length(found) > 0) df <- df[, found, drop = FALSE]

            # Compute Coverage = mapped reads / total sgRNAs in library
            total_sgrna <- NA_real_
            sgrna_file  <- rv$config$orig_sgrna_list_file
            if (!is.null(sgrna_file) && nzchar(sgrna_file)) {
                full_path <- sgrna_file
                if (!file.exists(full_path) && !is.null(rv$config_dir))
                    full_path <- file.path(rv$config_dir, sgrna_file)
                if (file.exists(full_path))
                    total_sgrna <- length(readLines(full_path, warn = FALSE)) - 1L
            }
            mapped_idx <- which(tolower(names(df)) == "mapped")[1]
            if (!is.na(total_sgrna) && !is.na(mapped_idx)) {
                coverage <- round(df[[mapped_idx]] / total_sgrna, 1)
                df <- cbind(df[, seq_len(mapped_idx),          drop = FALSE],
                            Coverage = coverage,
                            df[, seq(mapped_idx + 1L, ncol(df)), drop = FALSE])
            }

            # Resolve actual column names (case may vary from MAGeCK version)
            col_name <- function(key) names(df)[which(tolower(names(df)) == key)[1]]
            pct_col   <- col_name("percentage")
            gini_col  <- col_name("giniindex")
            zero_col  <- col_name("zerocounts")

            dt <- DT::datatable(df, rownames = FALSE,
                options = list(
                    paging    = FALSE,
                    searching = FALSE,
                    dom       = "t",
                    scrollX   = TRUE),
                caption = paste0("MAGeCK Count Summary — ", rv$experiment_name %||% ""))

            # Color rules: green = good, yellow = caution, red = bad
            green  <- "#d4edda"; yellow <- "#fff3cd"; red <- "#f8d7da"
            if (!is.na(pct_col))
                dt <- dt %>% DT::formatStyle(pct_col,
                    backgroundColor = DT::styleInterval(c(50, 70), c(red, yellow, green)))
            blue <- "#cce5ff"
            if ("Coverage" %in% names(df))
                dt <- dt %>% DT::formatStyle("Coverage",
                    backgroundColor = DT::styleInterval(c(100, 300, 1000), c(red, yellow, blue, green)))
            if (!is.na(gini_col))
                dt <- dt %>% DT::formatStyle(gini_col,
                    backgroundColor = DT::styleInterval(c(0.1, 0.2), c(green, yellow, red)))
            if (!is.na(zero_col))
                dt <- dt %>% DT::formatStyle(zero_col,
                    backgroundColor = DT::styleInterval(c(1, 500), c(green, yellow, red)))
            dt
        }, server = FALSE)
    } else {
        output$qc_countsummary_table <- renderTable({
            req(rv$qc_countsummary)
            df <- rv$qc_countsummary
            cn <- tolower(names(df))
            key_cols <- c("label","reads","mapped","percentage","zerocounts","giniindex")
            found_idx <- which(cn %in% key_cols)
            if (length(found_idx) > 0) df[, found_idx, drop = FALSE] else df
        }, striped = TRUE, hover = TRUE, bordered = TRUE)
    }

    # Shared renderImage helper: returns a list suitable for renderImage().
    # If filepath is NULL or does not exist, generates a small placeholder PNG.
    qc_image_result <- function(filepath, alt = "") {
        if (!is.null(filepath) && file.exists(filepath))
            return(list(src = filepath, contentType = "image/png", width = "100%",
                        alt = alt, deleteFile = FALSE))
        tmp <- tempfile(fileext = ".png")
        png(tmp, width = 400, height = 100)
        par(mar = c(0,0,0,0)); plot.new()
        text(0.5, 0.5, if (nzchar(alt)) paste0("Not found: ", alt) else "Image not found", cex = 1.2)
        dev.off()
        list(src = tmp, contentType = "image/png", width = "100%", alt = "Not found", deleteFile = TRUE)
    }

    output$qc_hist_status <- renderUI({
        if (is.null(rv$qc_histogram_meta) || length(rv$qc_histogram_meta) == 0)
            tags$p(style = "color:#888;", "No histogram PNG files found in the output directory.")
    })
    output$qc_hist_image <- renderImage({
        req(rv$qc_histogram_meta, input$qc_hist_sample, input$qc_hist_type)
        m <- Filter(function(x) x$sample == input$qc_hist_sample &&
                                x$type   == input$qc_hist_type, rv$qc_histogram_meta)
        qc_image_result(if (length(m) > 0) m[[1]]$filepath else NULL,
                        paste(input$qc_hist_sample, input$qc_hist_type))
    }, deleteFile = FALSE)

    output$qc_pca_status <- renderUI({
        if (is.null(rv$qc_pca_files) || length(rv$qc_pca_files) == 0)
            tags$p(style = "color:#888;", "No PCA PNG files found in the output directory.")
    })
    output$qc_pca_image <- renderImage({
        req(rv$qc_pca_files, input$qc_pca_file)
        qc_image_result(input$qc_pca_file, basename(input$qc_pca_file))
    }, deleteFile = FALSE)

    output$qc_barplot_status <- renderUI({
        if (is.null(rv$qc_barplot_files) || length(rv$qc_barplot_files) == 0)
            tags$p(style = "color:#888;", "No barplot PNG files found in the output directory.")
    })
    output$qc_barplot_image <- renderImage({
        req(rv$qc_barplot_files, input$qc_barplot_sample)
        qc_image_result(input$qc_barplot_sample, basename(input$qc_barplot_sample))
    }, deleteFile = FALSE)

    # ── Comparisons tab — existing server logic ───────────────────────────────

    observeEvent(input$comparison, {
        req(input$comparison, rv$result_files)
        result_info <- rv$result_files[[input$comparison]]
        if (!is.null(result_info) && file.exists(result_info$gene_file)) {
            rv$gene_data <- read.delim(result_info$gene_file, sep = "\t", header = TRUE)
            has_sets <- !is.null(result_info$plot_config) &&
                !is.null(result_info$plot_config$highlight_genes) &&
                length(result_info$plot_config$highlight_genes) > 0
            rv$has_highlight_sets <- has_sets
            if (!has_sets && (is.null(rv$custom_gene_list) || length(rv$custom_gene_list$genes) == 0))
                updateCheckboxInput(session, "highlight_significant", value = TRUE)
            if (has_sets) {
                for (i in seq_along(result_info$plot_config$highlight_genes)) {
                    group_info <- result_info$plot_config$highlight_genes[[i]]
                    group_name <- group_info$group_name %||% paste0("Group ", i)
                    if (is.null(rv$gene_sets_enabled[[group_name]])) rv$gene_sets_enabled[[group_name]] <- TRUE
                }
            }
        }
    })

    observeEvent(input$clear_custom, {
        rv$custom_gene_list <- NULL
        updateTextAreaInput(session, "custom_genes", value = "")
        rv$warning_message <- NULL
        if (!is.null(input$comparison) && !is.null(rv$result_files)) {
            result_info <- rv$result_files[[input$comparison]]
            has_sets <- !is.null(result_info$plot_config) && !is.null(result_info$plot_config$highlight_genes) && length(result_info$plot_config$highlight_genes) > 0
            if (!has_sets) updateCheckboxInput(session, "highlight_significant", value = TRUE)
        }
    })

    output$gene_set_checkboxes <- renderUI({
        req(rv$config)
        sys_groups <- get_system_gene_groups(rv$config, rv$config_dir, rv$neg_ctrl_override)
        sys_ui <- lapply(Filter(function(sg) sg$key != "neg_ctrl", sys_groups), function(sg) make_sys_group_row("comp", sg))

        pred_ui <- NULL
        if (!is.null(input$comparison) && !is.null(rv$result_files)) {
            gene_sets <- rv$result_files[[input$comparison]]$plot_config$highlight_genes
            if (!is.null(gene_sets) && length(gene_sets) > 0) {
                pred_controls <- lapply(seq_along(gene_sets), function(i) {
                    group_name <- gene_sets[[i]]$group_name %||% paste0("Group ", i)
                    div(style="border:1px solid #ddd;padding:5px;margin-bottom:5px;border-radius:3px;",
                        fluidRow(
                            column(3, checkboxInput(paste0("geneset_",i), label=NULL,
                                                    value=isTRUE(isolate(input[[paste0("geneset_",i)]]) %||% TRUE))),
                            column(9, textInput(paste0("geneset_name_",i), NULL,
                                               value=group_name, placeholder="Name"))
                        ),
                        fluidRow(
                            column(6, div(class="palette-dropdown",
                                selectInput(paste0("geneset_palette_",i), NULL,
                                    choices=c("Basic"="limited","Web-safe"="websafe","Full"="square"),
                                    selected="limited"))),
                            column(6, uiOutput(paste0("geneset_color_picker_",i)))
                        )
                    )
                })
                pred_ui <- tagList(
                    hr(style="margin:6px 0;"),
                    tags$p(tags$strong("From config:"), style="margin:4px 0;font-size:0.9em;"),
                    do.call(tagList, pred_controls)
                )
            }
        }
        do.call(tagList, c(sys_ui, list(pred_ui)))
    })

    output$fdr_color_picker_1 <- renderUI({
        palette_type <- if (is.null(input$fdr_palette_1)) "limited" else input$fdr_palette_1
        current_color <- isolate(if (!is.null(input$fdr_color_1)) input$fdr_color_1 else DEFAULT_FDR_COLOR_1)
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput("fdr_color_1", NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=FALSE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
        } else
            colourpicker::colourInput("fdr_color_1", NULL, value=current_color, palette=palette_type, returnName=FALSE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
    })

    output$fdr_color_picker_2 <- renderUI({
        palette_type <- if (is.null(input$fdr_palette_2)) "limited" else input$fdr_palette_2
        current_color <- isolate(if (!is.null(input$fdr_color_2)) input$fdr_color_2 else DEFAULT_FDR_COLOR_2)
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput("fdr_color_2", NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=FALSE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
        } else
            colourpicker::colourInput("fdr_color_2", NULL, value=current_color, palette=palette_type, returnName=FALSE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
    })

    outputOptions(output, "fdr_color_picker_1", suspendWhenHidden = FALSE)
    outputOptions(output, "fdr_color_picker_2", suspendWhenHidden = FALSE)

    get_fdr_color_1 <- reactive({ color <- input$fdr_color_1; if (is.null(color)) DEFAULT_FDR_COLOR_1 else color })
    get_fdr_color_2 <- reactive({
        color <- input$fdr_color_2
        if (is.null(color)) return(DEFAULT_FDR_COLOR_2)
        if (color == "#000000") return(DEFAULT_FDR_COLOR_2)
        return(color)
    })

    observe({
        req(!is.null(input$fdr_color_2))
        if (input$fdr_color_2 == "#000000") {
            js_code <- sprintf("
                var picker = $('#fdr_color_2');
                if (picker.length) { picker.val('%s'); picker.css('background-color','%s'); picker.trigger('change'); }
            ", DEFAULT_FDR_COLOR_2, DEFAULT_FDR_COLOR_2)
            shinyjs::runjs(js_code)
        }
    })

    output$sig_color_picker_1 <- renderUI({
        palette_type <- if (is.null(input$sig_palette_1)) "limited" else input$sig_palette_1
        current_color <- if (!is.null(input$sig_color_1)) input$sig_color_1 else rv$sig_color_1
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput("sig_color_1", NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=TRUE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
        } else
            colourpicker::colourInput("sig_color_1", NULL, value=current_color, palette=palette_type, returnName=TRUE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
    })

    output$sig_color_picker_2 <- renderUI({
        palette_type <- if (is.null(input$sig_palette_2)) "limited" else input$sig_palette_2
        current_color <- if (!is.null(input$sig_color_2)) input$sig_color_2 else rv$sig_color_2
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput("sig_color_2", NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=TRUE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
        } else
            colourpicker::colourInput("sig_color_2", NULL, value=current_color, palette=palette_type, returnName=TRUE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
    })

    observeEvent(input$sig_color_1, { if (!is.null(input$sig_color_1)) rv$sig_color_1 <- input$sig_color_1 })
    observeEvent(input$sig_color_2, { if (!is.null(input$sig_color_2)) rv$sig_color_2 <- input$sig_color_2 })

    observeEvent(input$use_plot_thresholds, {
        if (is.null(input$use_plot_thresholds)) return()
        if (!input$use_plot_thresholds) {
            rv$sig_color_1 <- get_fdr_color_1(); rv$sig_color_2 <- get_fdr_color_2()
        } else {
            # Reset threshold checkboxes to defaults each time the customize panel is opened
            updateCheckboxInput(session, "sig_use_thresh1", value = TRUE)
            updateCheckboxInput(session, "sig_use_thresh2", value = FALSE)
        }
    })

    output$custom_color_picker <- renderUI({
        palette_type <- if (is.null(input$custom_palette_type)) "square" else input$custom_palette_type
        current_color <- if (is.null(input$custom_color)) "#00FF00" else input$custom_color
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourInput("custom_color", NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=TRUE, showColour="both", width=PICKER_WIDTH)
        } else if (palette_type == "limited")
            colourInput("custom_color", NULL, value=current_color, palette="limited", returnName=TRUE, showColour="both", width=PICKER_WIDTH)
        else
            colourInput("custom_color", NULL, value=current_color, palette="square", returnName=TRUE, showColour="both", width=PICKER_WIDTH)
    })

    outputOptions(output, "custom_color_picker", suspendWhenHidden = FALSE)

    observe({
        req(input$comparison, rv$result_files)
        result_info <- rv$result_files[[input$comparison]]
        if (!is.null(result_info$plot_config) && !is.null(result_info$plot_config$highlight_genes)) {
            gene_sets <- result_info$plot_config$highlight_genes
            lapply(seq_along(gene_sets), function(i) {
                output[[paste0("geneset_color_picker_",i)]] <- renderUI({
                    palette_type <- input[[paste0("geneset_palette_",i)]]; if (is.null(palette_type)) palette_type <- "square"
                    current_color <- input[[paste0("geneset_color_",i)]]; if (is.null(current_color)) current_color <- gene_sets[[i]]$color %||% GENESET_COLOR_CYCLE[[(i-1L)%%length(GENESET_COLOR_CYCLE)+1L]]
                    if (palette_type == "websafe") {
                        colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
                        colourInput(paste0("geneset_color_",i), NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=TRUE, showColour="both", width=PICKER_WIDTH)
                    } else if (palette_type == "limited")
                        colourInput(paste0("geneset_color_",i), NULL, value=current_color, palette="limited", returnName=TRUE, showColour="both", width=PICKER_WIDTH)
                    else
                        colourInput(paste0("geneset_color_",i), NULL, value=current_color, palette="square", returnName=TRUE, showColour="both", width=PICKER_WIDTH)
                })
            })
        }
    })

    # ── System group color pickers — Comparisons tab (neg_ctrl excluded) ────
    for (.key in setdiff(names(SYS_GROUP_DEFAULTS), "neg_ctrl")) {
        local({
            k       <- .key
            def_col <- SYS_GROUP_DEFAULTS[[k]]
            prefix  <- paste0("comp_sys_", k)
            output[[paste0(prefix, "_color_picker")]] <- renderUI({
                pal <- input[[paste0(prefix, "_palette")]] %||% "square"
                cur <- input[[paste0(prefix, "_color")]]  %||% def_col
                if (pal == "websafe") {
                    cvec <- if (identical(input$websafe_order, "hsv")) websafe_colors_hsv else websafe_colors
                    colourInput(paste0(prefix, "_color"), NULL, value=cur, palette="limited",
                                allowedCols=cvec, returnName=TRUE, showColour="both", width=PICKER_WIDTH)
                } else
                    colourInput(paste0(prefix, "_color"), NULL, value=cur, palette=pal,
                                returnName=TRUE, showColour="both", width=PICKER_WIDTH)
            })
            outputOptions(output, paste0(prefix, "_color_picker"), suspendWhenHidden=FALSE)
        })
    }

    # ── Negative control ID modal ─────────────────────────────────────────────
    show_neg_ctrl_modal <- function() {
        showModal(modalDialog(
            title = "Set Negative Control Gene ID",
            tags$p("Enter the gene/target ID for negative control sgRNAs:",
                   style="font-size:0.9em;color:#555;margin-bottom:8px;"),
            textInput("neg_ctrl_id_input", NULL,
                      value = rv$neg_ctrl_override %||% rv$config$control_gene_id %||% "",
                      placeholder = "e.g. Non-Targeting_Control"),
            footer = tagList(
                modalButton("Cancel"),
                actionButton("neg_ctrl_id_save", "Save", class="btn-primary btn-sm")
            ),
            size = "s", easyClose = TRUE
        ))
    }
    observeEvent(input$counts_neg_ctrl_edit, { show_neg_ctrl_modal() })
    observeEvent(input$neg_ctrl_id_save, {
        rv$neg_ctrl_override <- trimws(input$neg_ctrl_id_input %||% "")
        removeModal()
    })

    observeEvent(input$custom_gene_file, {
        req(input$custom_gene_file)
        file_content <- readLines(input$custom_gene_file$datapath, warn = FALSE)
        genes <- trimws(file_content); genes <- genes[genes != ""]
        if (length(genes) > 0) updateTextAreaInput(session, "custom_genes", value = paste(genes, collapse = "\n"))
    })

    observeEvent(input$apply_custom, {
        req(input$custom_genes)
        raw <- gsub("\r\n", "\n", input$custom_genes); raw <- gsub("\r", "\n", raw)
        lines <- unlist(strsplit(raw, "\n")); lines <- lines[trimws(lines) != ""]
        raw_norm <- raw
        try_split <- function(text, split_func) {
            toks <- split_func(text); toks <- trimws(toks); toks <- toks[toks != ""]
            bad <- any(grepl("[[[:space:]],;:]", toks))
            list(tokens = toks, ok = !bad)
        }
        candidates_try <- list(
            newline   = try_split(raw_norm, function(t) unlist(strsplit(t, "\n"))),
            comma     = try_split(raw_norm, function(t) unlist(strsplit(t, ","))),
            semicolon = try_split(raw_norm, function(t) unlist(strsplit(t, ";"))),
            colon     = try_split(raw_norm, function(t) unlist(strsplit(t, ":"))),
            whitespace= try_split(raw_norm, function(t) unlist(strsplit(t, "\\s+")))
        )
        chosen <- NULL; tokens <- character(0)
        for (d in c("newline","comma","semicolon","colon","whitespace")) {
            if (length(candidates_try[[d]]$tokens) > 0 && candidates_try[[d]]$ok) { chosen <- d; tokens <- candidates_try[[d]]$tokens; break }
        }
        if (is.null(chosen)) { chosen <- "whitespace"; tokens <- candidates_try$whitespace$tokens }
        tokens <- gsub("^[\"'`,;:[[:space:]]]+|[\"'`,;:[[:space:]]]+$", "", tokens)
        tokens <- tokens[tokens != ""]
        all_tokens <- tokens
        lines_norm <- trimws(unlist(strsplit(raw_norm, "\n"))); lines_norm <- lines_norm[lines_norm != ""]
        line_last <- if (length(lines_norm) > 0) {
            ll <- sapply(lines_norm, function(l) {
                if (grepl("[,;:]", l)) { parts <- trimws(unlist(strsplit(l, "[,;:]+"))); parts <- parts[parts != ""]; if (length(parts) == 0) NA_character_ else tail(parts,1) } else l
            }); ll <- trimws(ll); ll[!is.na(ll) & ll != ""]
        } else character(0)
        if (length(all_tokens) == 0 && length(line_last) == 0) { rv$warning_message <- "No tokens parsed from input."; rv$custom_gene_list <- NULL; return() }
        candidates <- unique(c(line_last, all_tokens))
        candidates <- candidates[!grepl("^[0-9]+$", candidates)]
        candidates <- candidates[!grepl("^[ACGTNacgtn]{15,}$", candidates)]
        candidates <- candidates[!grepl("^.*_[0-9]+$", candidates)]
        if (length(candidates) == 0) { rv$warning_message <- paste0("After filtering, no gene-like tokens left."); rv$custom_gene_list <- NULL; return() }
        genes <- gsub("[^[:print:]]", "", candidates); genes <- trimws(genes); genes <- genes[genes != ""]
        if (is.null(rv$gene_data)) {
            rv$custom_gene_list <- list(name = input$custom_list_name %||% "Custom", genes = unique(genes), color = input$custom_color)
            rv$warning_message <- paste0("Loaded ", length(unique(genes)), " candidate genes; no dataset loaded for validation.")
            return()
        }
        colnames_clean <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(rv$gene_data)))
        gene_data_temp <- rv$gene_data; colnames(gene_data_temp) <- colnames_clean
        candidate_cols <- unique(c(intersect(c("id","gene"), colnames_clean), grep("gene|symbol|name", colnames_clean, value=TRUE)))
        if (length(candidate_cols) == 0) candidate_cols <- colnames_clean
        best_col <- NULL; best_hits <- -1; best_avail <- NULL
        for (col in candidate_cols) { avail <- gene_data_temp[[col]]; hits <- sum(genes %in% avail, na.rm=TRUE); if (hits > best_hits) { best_hits <- hits; best_col <- col; best_avail <- avail } }
        gene_id_col <- if (!is.null(best_col)) best_col else colnames_clean[1]
        available_genes <- if (!is.null(best_avail)) best_avail else gene_data_temp[[gene_id_col]]
        matched_mask <- genes %in% available_genes
        matched_tokens <- genes[matched_mask]
        matched_genes <- if (length(matched_tokens) > 0) unique(available_genes[match(matched_tokens, available_genes)]) else character(0)
        unmatched_genes <- genes[!matched_mask]
        if (length(matched_genes) == 0) {
            rv$warning_message <- paste0("None of ", length(genes), " candidate tokens matched genes in dataset.")
            rv$custom_gene_list <- NULL
        } else {
            rv$custom_gene_list <- list(name = input$custom_list_name, genes = unique(matched_genes), color = input$custom_color)
            rv$warning_message <- if (length(unmatched_genes) > 0)
                paste0("Matched ", length(matched_genes), " genes. ", length(unmatched_genes), " unmatched: ", paste(head(unmatched_genes, 20), collapse=" "), if (length(unmatched_genes) > 20) paste0(" (+", length(unmatched_genes)-20, " more)"))
            else paste0("Successfully loaded ", length(matched_genes), " unique genes.")
        }
    })

    output$custom_warning <- renderUI({
        if (!is.null(rv$warning_message)) {
            warning_color <- if (grepl("Successfully", rv$warning_message)) "green" else if (grepl("None|No tokens", rv$warning_message)) "red" else "orange"
            div(style = paste0("margin-top:10px;padding:8px;border-radius:4px;background-color:", warning_color, "22;border:1px solid ", warning_color, ";"),
                p(style = paste0("margin:0;color:", warning_color, ";font-size:12px;"), rv$warning_message))
        }
    })

    get_valid_fdr_threshold <- function(value, default) {
        if (is.null(value) || value == "") return(default)
        num_val <- suppressWarnings(as.numeric(value))
        if (is.na(num_val) || num_val < 0 || num_val > 1) return(default)
        return(num_val)
    }

    get_significant_genes <- function(gene_data, result_type, thresh1, thresh2, use_thresh1, use_thresh2, top_n = 10, rank_by_rra = FALSE, reverse_rank = FALSE) {
        if (result_type == "neg") { fdr_col <- "neg_fdr"; lfc_col <- "neg_lfc"; rank_col <- "neg_rank" }
        else { fdr_col <- "pos_fdr"; lfc_col <- "pos_lfc"; rank_col <- "pos_rank" }
        colnames(gene_data) <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(gene_data)))
        fdr_col <- tolower(fdr_col); lfc_col <- tolower(lfc_col); rank_col <- tolower(rank_col)
        gene_id_col <- if ("id" %in% colnames(gene_data)) "id" else if ("gene" %in% colnames(gene_data)) "gene" else colnames(gene_data)[1]
        sort_data <- function(df) {
            if (rank_by_rra && rank_col %in% colnames(df))
                df %>% arrange(.data[[fdr_col]], if (reverse_rank) desc(.data[[rank_col]]) else .data[[rank_col]])
            else if (result_type == "neg")
                df %>% arrange(.data[[fdr_col]], if (reverse_rank) desc(.data[[lfc_col]]) else .data[[lfc_col]])
            else
                df %>% arrange(.data[[fdr_col]], if (reverse_rank) .data[[lfc_col]] else desc(.data[[lfc_col]]))
        }
        sig_groups <- list()
        if (!use_thresh1 && !use_thresh2) {
            # No FDR filter — label top N by rank regardless of significance
            sort_no_fdr <- function(df) {
                if (rank_by_rra && rank_col %in% colnames(df))
                    df %>% arrange(if (reverse_rank) desc(.data[[rank_col]]) else .data[[rank_col]])
                else if (result_type == "neg")
                    df %>% arrange(if (reverse_rank) desc(.data[[lfc_col]]) else .data[[lfc_col]])
                else
                    df %>% arrange(if (reverse_rank) .data[[lfc_col]] else desc(.data[[lfc_col]]))
            }
            genes_top <- gene_data %>% filter(!is.na(.data[[lfc_col]])) %>% sort_no_fdr() %>% head(top_n) %>% pull(.data[[gene_id_col]])
            if (length(genes_top) > 0) sig_groups <- list(list(genes = genes_top, threshold = NULL, slot = 1))
        } else {
            if (use_thresh1) {
                genes_t1 <- gene_data %>% filter(!is.na(.data[[fdr_col]]) & .data[[fdr_col]] < thresh1) %>%
                    sort_data() %>% head(top_n) %>% pull(.data[[gene_id_col]])
                if (length(genes_t1) > 0) sig_groups <- c(sig_groups, list(list(genes = genes_t1, threshold = thresh1, slot = 1)))
            }
            if (use_thresh2) {
                genes_t2 <- gene_data %>% filter(!is.na(.data[[fdr_col]]) & .data[[fdr_col]] >= thresh1 & .data[[fdr_col]] < thresh2) %>%
                    sort_data() %>% head(top_n) %>% pull(.data[[gene_id_col]])
                if (length(genes_t2) > 0) sig_groups <- c(sig_groups, list(list(genes = genes_t2, threshold = thresh2, slot = 2)))
            }
        }
        return(sig_groups)
    }

    significant_genes_data <- reactive({
        req(input$highlight_significant, input$selection_type, rv$gene_data)
        use_custom <- isTRUE(input$use_plot_thresholds)
        if (!use_custom) {
            # Default (unchecked): mirror the Plot Options FDR thresholds and colors exactly
            thresh1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05)
            thresh2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
            use_t1  <- TRUE
            use_t2  <- FALSE
            color1  <- get_fdr_color_1()
            color2  <- get_fdr_color_2()
        } else {
            # Customized (checked): use the independent threshold inputs and per-threshold enable toggles
            thresh1 <- get_valid_fdr_threshold(input$sig_fdr_threshold_1, 0.05)
            thresh2 <- get_valid_fdr_threshold(input$sig_fdr_threshold_2, 0.1)
            use_t1  <- if (is.null(input$sig_use_thresh1)) TRUE else input$sig_use_thresh1
            use_t2  <- if (is.null(input$sig_use_thresh2)) FALSE else input$sig_use_thresh2
            color1  <- if (!is.null(input$sig_color_1)) input$sig_color_1 else rv$sig_color_1
            color2  <- if (!is.null(input$sig_color_2)) input$sig_color_2 else rv$sig_color_2
        }
        top_n <- if (is.null(input$sig_top_n) || is.na(input$sig_top_n)) 10 else as.integer(input$sig_top_n)
        sig_groups <- get_significant_genes(rv$gene_data, input$selection_type, thresh1, thresh2, use_t1, use_t2, top_n,
                                            isTRUE(input$rank_method == "rra"), isTRUE(input$reverse_rank_order))
        list(sig_groups = sig_groups, thresh2 = thresh2, color1 = color1, color2 = color2)
    })

    observe({
        req(input$highlight_significant)
        if (input$highlight_significant) {
            sig_data <- significant_genes_data(); highlight_config <- get_highlight_config()
            both_fdr_off <- isTRUE(input$use_plot_thresholds) && !isTRUE(input$sig_use_thresh1) && !isTRUE(input$sig_use_thresh2)
            if ((is.null(highlight_config) || length(highlight_config) == 0) && length(sig_data$sig_groups) == 0)
                rv$sig_warning <- if (both_fdr_off) "No gene data available." else
                    paste0("No genes with FDR < ", sig_data$thresh2, " for ", if (input$selection_type == "neg") "negative" else "positive", " selection.")
            else rv$sig_warning <- NULL
        } else rv$sig_warning <- NULL
    })

    get_highlight_config <- reactive({
        req(input$comparison, rv$result_files)
        rank_method <- input$rank_method; reverse_rank <- input$reverse_rank_order
        result_info <- rv$result_files[[input$comparison]]
        highlight_config <- list()

        # System groups (essential, non-essential only — neg_ctrl is sgRNA-level only)
        sys_groups <- get_system_gene_groups(rv$config, rv$config_dir, rv$neg_ctrl_override)
        for (sg in Filter(function(sg) sg$key != "neg_ctrl", sys_groups)) {
            pid <- paste0("comp_sys_", sg$key)
            if (!isTRUE(input[[pid]])) next
            entry <- sg
            col <- input[[paste0(pid, "_color")]]; if (!is.null(col)) entry$color <- col
            highlight_config <- c(highlight_config, list(entry))
        }

        if (!is.null(result_info$plot_config) && !is.null(result_info$plot_config$highlight_genes)) {
            gene_sets <- result_info$plot_config$highlight_genes
            for (i in seq_along(gene_sets)) {
                checkbox_id <- paste0("geneset_",i)
                if (!is.null(input[[checkbox_id]]) && input[[checkbox_id]]) {
                    custom_name <- input[[paste0("geneset_name_",i)]]
                    custom_color <- input[[paste0("geneset_color_",i)]]
                    modified_set <- gene_sets[[i]]
                    if (!is.null(custom_name) && custom_name != "") modified_set$group_name <- custom_name
                    if (!is.null(custom_color)) modified_set$color <- custom_color
                    highlight_config <- c(highlight_config, list(modified_set))
                }
            }
        }
        if (!is.null(rv$custom_gene_list) && !is.null(input$show_custom) && input$show_custom)
            highlight_config <- c(highlight_config, list(list(group_name=rv$custom_gene_list$name, color=rv$custom_gene_list$color, genes=rv$custom_gene_list$genes, is_custom=TRUE)))
        if (!is.null(input$highlight_significant) && input$highlight_significant && !is.null(rv$gene_data)) {
            sig_data <- significant_genes_data()
            for (i in seq_along(sig_data$sig_groups)) {
                group <- sig_data$sig_groups[[i]]
                highlight_config <- c(highlight_config, list(list(group_name=if(is.null(group$threshold)) "Top N" else paste0("Top N (FDR < ", group$threshold, ")"), color=if(isTRUE(group$slot==1)) sig_data$color1 else sig_data$color2, genes=group$genes, is_custom=TRUE)))
            }
        }
        if (length(highlight_config) > 0) highlight_config else NULL
    })

    output$global_sig_warning <- renderUI({
        if (isTRUE(input$highlight_significant) && !is.null(rv$sig_warning) && nzchar(rv$sig_warning))
            tags$div(style = "margin-top:6px;padding:6px;background-color:#fff3cd;color:#856404;border:1px solid #ffeeba;border-radius:4px;font-size:0.9em;", rv$sig_warning)
    })
    outputOptions(output, "global_sig_warning", suspendWhenHidden = FALSE)

    # ── Counts tab ────────────────────────────────────────────────────────────

    observeEvent(input$counts_id_color, { if (!is.null(input$counts_id_color)) rv$counts_id_color <- input$counts_id_color })
    observeEvent(input$counts_hl_color,  { if (!is.null(input$counts_hl_color))  rv$counts_hl_color  <- input$counts_hl_color  })

    output$counts_id_color_picker <- renderUI({
        palette_type  <- input$counts_id_palette %||% "limited"
        current_color <- rv$counts_id_color
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput("counts_id_color", NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=FALSE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
        } else
            colourpicker::colourInput("counts_id_color", NULL, value=current_color, palette=palette_type, returnName=FALSE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
    })

    output$counts_hl_color_picker <- renderUI({
        palette_type  <- input$counts_hl_palette %||% "limited"
        current_color <- rv$counts_hl_color
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput("counts_hl_color", NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=FALSE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
        } else
            colourpicker::colourInput("counts_hl_color", NULL, value=current_color, palette=palette_type, returnName=FALSE, showColour="both", allowTransparent=FALSE, width=PICKER_WIDTH)
    })
    outputOptions(output, "counts_hl_color_picker", suspendWhenHidden = FALSE)

    # ── System group color pickers — Counts tab ───────────────────────────────
    for (.key in names(SYS_GROUP_DEFAULTS)) {
        local({
            k       <- .key
            def_col <- SYS_GROUP_DEFAULTS[[k]]
            prefix  <- paste0("counts_sys_", k)
            output[[paste0(prefix, "_color_picker")]] <- renderUI({
                pal <- input[[paste0(prefix, "_palette")]] %||% "square"
                cur <- input[[paste0(prefix, "_color")]]  %||% def_col
                if (pal == "websafe") {
                    cvec <- if (identical(input$websafe_order, "hsv")) websafe_colors_hsv else websafe_colors
                    colourInput(paste0(prefix, "_color"), NULL, value=cur, palette="limited",
                                allowedCols=cvec, returnName=TRUE, showColour="both", width=PICKER_WIDTH)
                } else
                    colourInput(paste0(prefix, "_color"), NULL, value=cur, palette=pal,
                                returnName=TRUE, showColour="both", width=PICKER_WIDTH)
            })
            outputOptions(output, paste0(prefix, "_color_picker"), suspendWhenHidden=FALSE)
        })
    }

    # Collect unique predefined gene sets across all comparisons (deduped by group_name)
    counts_predefined_gene_sets <- reactive({
        if (is.null(rv$result_files)) return(list())
        all_sets <- list(); seen <- character(0)
        for (result_info in rv$result_files) {
            gene_sets <- result_info$plot_config$highlight_genes
            if (is.null(gene_sets)) next
            for (gs in gene_sets) {
                gname <- gs$group_name %||% "Unknown"
                if (!(gname %in% seen)) { all_sets <- c(all_sets, list(gs)); seen <- c(seen, gname) }
            }
        }
        all_sets
    })

    # Register color pickers for each predefined gene set in the Counts tab
    observe({
        gene_sets <- counts_predefined_gene_sets()
        lapply(seq_along(gene_sets), function(i) {
            output[[paste0("counts_geneset_color_picker_", i)]] <- renderUI({
                palette_type  <- input[[paste0("counts_geneset_palette_", i)]] %||% "square"
                current_color <- input[[paste0("counts_geneset_color_",   i)]] %||% (gene_sets[[i]]$color %||% GENESET_COLOR_CYCLE[[(i-1L)%%length(GENESET_COLOR_CYCLE)+1L]])
                if (palette_type == "websafe") {
                    colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
                    colourInput(paste0("counts_geneset_color_",i), NULL, value=current_color, palette="limited", allowedCols=colors_vec, returnName=TRUE, showColour="both", width=PICKER_WIDTH)
                } else if (palette_type == "limited")
                    colourInput(paste0("counts_geneset_color_",i), NULL, value=current_color, palette="limited", returnName=TRUE, showColour="both", width=PICKER_WIDTH)
                else
                    colourInput(paste0("counts_geneset_color_",i), NULL, value=current_color, palette="square", returnName=TRUE, showColour="both", width=PICKER_WIDTH)
            })
        })
    })

    output$counts_gene_set_checkboxes <- renderUI({
        sys_groups <- get_system_gene_groups(rv$config, rv$config_dir, rv$neg_ctrl_override)
        gene_sets  <- counts_predefined_gene_sets()

        sys_ui <- lapply(sys_groups, function(sg) make_sys_group_row("counts", sg))

        pred_ui <- if (length(gene_sets) > 0) {
            set_controls <- lapply(seq_along(gene_sets), function(i) {
                group_name <- gene_sets[[i]]$group_name %||% paste0("Group ", i)
                div(style="border:1px solid #ddd;padding:5px;margin-bottom:5px;border-radius:3px;",
                    fluidRow(
                        column(3, checkboxInput(paste0("counts_geneset_",i), label=NULL,
                                               value=isTRUE(isolate(input[[paste0("counts_geneset_",i)]]) %||% TRUE))),
                        column(9, textInput(paste0("counts_geneset_name_",i), NULL,
                                           value=group_name, placeholder="Name"))
                    ),
                    fluidRow(
                        column(6, div(class="palette-dropdown",
                            selectInput(paste0("counts_geneset_palette_",i), NULL,
                                choices=c("Basic"="limited","Web-safe"="websafe","Full"="square"),
                                selected="square"))),
                        column(6, uiOutput(paste0("counts_geneset_color_picker_",i)))
                    )
                )
            })
            tagList(
                if (length(sys_groups) > 0) hr(style="margin:6px 0;"),
                tags$p(tags$strong("From config:"), style="margin:4px 0;font-size:0.9em;"),
                do.call(tagList, set_controls)
            )
        } else NULL

        custom_ui <- tagList(
            if (length(sys_groups) > 0 || length(gene_sets) > 0) hr(style="margin:6px 0;"),
            tags$a(class="section-header", `data-target`="counts_custom_gene_content",
                   "Custom Gene List", tags$span(class="arrow", "\u25bc")),
            tags$div(id="counts_custom_gene_content", class="section-content",
                     style="display:none;",
                textInput("counts_hl_name", "Group label:", value="Highlight"),
                fluidRow(
                    column(5, div(class="palette-dropdown",
                        selectInput("counts_hl_palette", NULL,
                            choices=c("Basic"="limited","Web-safe"="websafe","Full"="square"),
                            selected="limited", width="100%"))),
                    column(7, uiOutput("counts_hl_color_picker"))
                ),
                textAreaInput("counts_hl_genes", "Gene list (one per line or comma-separated):",
                              placeholder="GENE1\nGENE2, GENE3", rows=4, width="100%")
            )
        )

        tagList(
            checkboxInput("counts_show_labels", "Label points", value=TRUE),
            conditionalPanel(
                condition="input.counts_show_labels",
                if (PLOTLY_AVAILABLE)
                    tags$p(style="font-size:0.8em;color:#888;margin:2px 0 4px 0;",
                           "\u2139 Labels only appear in static plot mode."),
                sliderInput("counts_label_scale", "Label scale:", min=0.5, max=2.5, value=1, step=0.1)
            ),
            sliderInput("counts_hl_point_size", "Highlight point size:",
                        min=0.1, max=5, value=1.5, step=0.1),
            hr(style="margin:8px 0;"),
            do.call(tagList, sys_ui),
            pred_ui,
            custom_ui
        )
    })

    get_counts_highlight_config <- reactive({
        highlight_config <- list()

        # System groups (neg ctrl, essential, non-essential)
        sys_groups <- get_system_gene_groups(rv$config, rv$config_dir, rv$neg_ctrl_override)
        for (sg in sys_groups) {
            pid <- paste0("counts_sys_", sg$key)
            if (!isTRUE(input[[pid]])) next
            entry <- sg
            col <- input[[paste0(pid, "_color")]]; if (!is.null(col)) entry$color <- col
            highlight_config <- c(highlight_config, list(entry))
        }

        # Predefined gene sets
        gene_sets <- counts_predefined_gene_sets()
        for (i in seq_along(gene_sets)) {
            if (!isTRUE(input[[paste0("counts_geneset_", i)]])) next
            modified <- gene_sets[[i]]
            custom_name  <- input[[paste0("counts_geneset_name_",  i)]]
            custom_color <- input[[paste0("counts_geneset_color_", i)]]
            if (!is.null(custom_name)  && nzchar(custom_name))  modified$group_name <- custom_name
            if (!is.null(custom_color))                          modified$color      <- custom_color
            highlight_config <- c(highlight_config, list(modified))
        }

        # Custom gene list
        raw <- input$counts_hl_genes %||% ""
        if (nzchar(trimws(raw))) {
            tokens <- trimws(unlist(strsplit(gsub("[,;|\t]+", "\n", raw), "\n")))
            genes  <- unique(tokens[tokens != ""])
            if (length(genes) > 0) {
                group_name  <- if (nzchar(trimws(input$counts_hl_name %||% ""))) input$counts_hl_name else "Highlight"
                group_color <- rv$counts_hl_color %||% "#E41A1C"
                highlight_config <- c(highlight_config,
                    list(list(group_name=group_name, color=group_color, genes=genes, is_custom=TRUE)))
            }
        }

        if (length(highlight_config) > 0) highlight_config else NULL
    })

    output$counts_has_data <- reactive({
        !is.null(rv$count_matrix_file) || !is.null(rv$cpm_matrix_file)
    })
    outputOptions(output, "counts_has_data", suspendWhenHidden = FALSE)

    output$counts_no_config_msg <- renderUI({
        if (is.null(rv$output_dir))
            tags$div(style="padding:15px;background:#fff3cd;border:1px solid #ffc107;border-radius:4px;margin-bottom:10px;",
                     tags$strong("No experiment loaded."), " Load a config.yaml using the Configuration tab.")
    })

    output$counts_matrix_status <- renderUI({
        has_count <- !is.null(rv$count_matrix_file)
        has_cpm   <- !is.null(rv$cpm_matrix_file)
        if (!has_count && !has_cpm)
            tags$p(style="color:#856404;font-size:0.85em;", "\u26a0 No matrix files found.")
        else if (!has_count)
            tags$p(style="color:#856404;font-size:0.85em;", "\u26a0 Count matrix not found; CPM only.")
        else if (!has_cpm)
            tags$p(style="color:#856404;font-size:0.85em;", "\u26a0 CPM matrix not found; counts only.")
        else NULL
    })

    observe({
        has_count <- !is.null(rv$count_matrix_file)
        has_cpm   <- !is.null(rv$cpm_matrix_file)
        choices <- c()
        if (has_count) choices <- c(choices, c("Counts"="count"))
        if (has_cpm)   choices <- c(choices, c("CPM"="cpm"))
        updateSelectInput(session, "counts_matrix_type", choices=choices,
                          selected=if (has_count) "count" else if (has_cpm) "cpm" else NULL)
    })

    counts_matrix_data <- reactive({
        type <- input$counts_matrix_type %||% "count"
        fp <- if (type == "count") rv$count_matrix_file else rv$cpm_matrix_file
        req(!is.null(fp))
        read.table(fp, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    })

    counts_sample_names <- reactive({
        df <- counts_matrix_data()
        colnames(df)[sapply(df, is.numeric)]
    })

    # Per-plot config UI (cascading collapsible sections)
    output$counts_plot_configs <- renderUI({
        samples <- counts_sample_names()
        ids <- rv$counts_plot_ids
        lapply(seq_along(ids), function(i) {
            pid <- ids[i]
            x_id <- paste0("counts_x_", pid)
            y_id <- paste0("counts_y_", pid)
            # Y is always the first sample (reference/baseline).
            # X cycles through samples[2], samples[3], … for successive plots,
            # so each new plot defaults to a different comparison.
            n_samp <- length(samples)
            default_y <- if (n_samp >= 1) samples[1] else NULL
            default_x <- if (n_samp >= 2) {
                idx <- ((i - 1L) %% (n_samp - 1L)) + 2L
                samples[idx]
            } else if (n_samp >= 1) samples[1] else NULL
            cur_x <- isolate(input[[x_id]]) %||% default_x
            cur_y <- isolate(input[[y_id]]) %||% default_y
            if (!is.null(cur_x) && !(cur_x %in% samples)) cur_x <- default_x
            if (!is.null(cur_y) && !(cur_y %in% samples)) cur_y <- default_y
            content_id <- paste0("counts_content_", pid)
            n_ids <- length(ids)
            plot_name <- rv$counts_plot_names[[pid]] %||% paste0("Plot ", i)
            rename_btn <- tags$span(
                style = "cursor:pointer;color:#bbb;font-size:0.85em;line-height:1;padding:0 3px;",
                title = "Rename",
                onclick = paste0("event.stopPropagation();Shiny.setInputValue('counts_rename_plot','",
                                 pid, "',{priority:'event'})"),
                HTML("&#9998;"))
            up_btn <- if (i > 1)
                tags$span(class = "counts-move-btn",
                    title = "Move up",
                    onclick = paste0("event.stopPropagation();Shiny.setInputValue('counts_move_up','",
                                     pid, "',{priority:'event'})"),
                    HTML("&#8593;"))
            else NULL
            down_btn <- if (i < n_ids)
                tags$span(class = "counts-move-btn",
                    title = "Move down",
                    onclick = paste0("event.stopPropagation();Shiny.setInputValue('counts_move_down','",
                                     pid, "',{priority:'event'})"),
                    HTML("&#8595;"))
            else NULL
            remove_btn <- if (n_ids > 1)
                tags$span(
                    style = "cursor:pointer;color:#bbb;font-size:1.1em;line-height:1;padding-left:4px;",
                    onclick = paste0("event.stopPropagation();Shiny.setInputValue('counts_remove_plot','",
                                     pid, "',{priority:'event'})"),
                    HTML("&times;")
                )
            else NULL
            tags$div(`data-plot-id` = pid,
                tags$a(class = "section-header counts-plot-header", `data-target` = content_id,
                       style = "display:flex; align-items:center;",
                       tags$span(style = "flex:1;", plot_name),
                       rename_btn, up_btn, down_btn, remove_btn,
                       tags$span(class = "arrow", style = "margin-left:10px;", "\u25bc")),
                tags$div(id = content_id, class = "section-content counts-plot-content", style = "display:block;",
                    tags$div(style = "display:flex; align-items:center; gap:5px; margin-bottom:3px;",
                        tags$label("X", style = "white-space:nowrap; flex-shrink:0; margin:0; font-size:0.85em; font-weight:600;"),
                        tags$div(style = "flex:1; min-width:0;",
                            selectInput(x_id, NULL, choices = samples, selected = cur_x, width = "100%"))
                    ),
                    tags$div(style = "display:flex; align-items:center; gap:5px; margin-bottom:2px;",
                        tags$label("Y", style = "white-space:nowrap; flex-shrink:0; margin:0; font-size:0.85em; font-weight:600;"),
                        tags$div(style = "flex:1; min-width:0;",
                            selectInput(y_id, NULL, choices = samples, selected = cur_y, width = "100%"))
                    )
                )
            )
        })
    })

    observeEvent(input$counts_add_plot, {
        rv$counts_plot_id_counter <- rv$counts_plot_id_counter + 1L
        new_id <- paste0("p", rv$counts_plot_id_counter)
        # Name is based on counter so it's unique even after deletions
        rv$counts_plot_names[[new_id]] <- paste0("Plot ", rv$counts_plot_id_counter)
        rv$counts_plot_ids <- c(rv$counts_plot_ids, new_id)
    })

    observeEvent(input$counts_remove_plot, {
        pid <- input$counts_remove_plot
        ids <- rv$counts_plot_ids
        if (length(ids) > 1) {
            rv$counts_plot_ids <- ids[ids != pid]
            rv$counts_plot_names[[pid]] <- NULL
        }
    })

    # Rename: show modal pre-filled with current name
    observeEvent(input$counts_rename_plot, {
        pid <- input$counts_rename_plot
        rv$counts_renaming_pid <- pid
        current_name <- rv$counts_plot_names[[pid]] %||% pid
        showModal(modalDialog(
            title = "Rename plot",
            textInput("counts_rename_input", NULL, value = current_name,
                      placeholder = "Enter plot name"),
            footer = tagList(
                modalButton("Cancel"),
                actionButton("counts_rename_confirm", "Rename", class = "btn-primary btn-sm")
            ),
            size = "s", easyClose = TRUE
        ))
    })

    observeEvent(input$counts_rename_confirm, {
        pid <- rv$counts_renaming_pid
        new_name <- trimws(input$counts_rename_input %||% "")
        if (!is.null(pid) && nzchar(new_name))
            rv$counts_plot_names[[pid]] <- new_name
        removeModal()
        rv$counts_renaming_pid <- NULL
    })

    # Reset names: renumber all plots by current display order
    observeEvent(input$counts_reset_names, {
        ids <- rv$counts_plot_ids
        for (i in seq_along(ids))
            rv$counts_plot_names[[ids[i]]] <- paste0("Plot ", i)
    })

    # Drag-to-reorder: SortableJS reports new ID order after a drag
    observeEvent(input$counts_plot_order, {
        new_order <- input$counts_plot_order
        if (is.character(new_order) && setequal(new_order, rv$counts_plot_ids))
            rv$counts_plot_ids <- new_order
    })

    # Fallback: ↑/↓ buttons when SortableJS is unavailable
    observeEvent(input$counts_move_up, {
        ids <- rv$counts_plot_ids
        idx <- match(input$counts_move_up, ids)
        if (!is.na(idx) && idx > 1L) {
            ids[c(idx - 1L, idx)] <- ids[c(idx, idx - 1L)]
            rv$counts_plot_ids <- ids
        }
    })
    observeEvent(input$counts_move_down, {
        ids <- rv$counts_plot_ids
        idx <- match(input$counts_move_down, ids)
        if (!is.na(idx) && idx < length(ids)) {
            ids[c(idx, idx + 1L)] <- ids[c(idx + 1L, idx)]
            rv$counts_plot_ids <- ids
        }
    })

    # Notify user about SortableJS availability
    output$counts_sortable_status <- renderUI({
        ok <- input$counts_sortable_ok
        if (is.null(ok) || isTRUE(ok)) return(NULL)
        tags$div(
            style = "padding:3px 6px;background:#fff3cd;border:1px solid #ffc107;border-radius:3px;font-size:0.8em;margin-bottom:4px;",
            "\u26a0 Drag-to-reorder unavailable (SortableJS not loaded \u2014 offline?). Use the \u2191\u2193 buttons to reorder."
        )
    })

    counts_plot_pairs <- reactive({
        ids <- rv$counts_plot_ids
        samples <- counts_sample_names()
        pairs <- lapply(seq_along(ids), function(i) {
            pid <- ids[i]
            x_id <- paste0("counts_x_", pid)
            y_id <- paste0("counts_y_", pid)
            x <- input[[x_id]]
            y <- input[[y_id]]
            nm <- rv$counts_plot_names[[pid]] %||% paste0("Plot ", i)
            if (!is.null(x) && x %in% samples && !is.null(y) && y %in% samples)
                list(x=x, y=y, name=nm)
            else NULL
        })
        Filter(Negate(is.null), pairs)
    })

    output$counts_log_warning <- renderUI({
        transform   <- input$counts_transform %||% "none"
        pseudocount <- as.numeric(input$counts_pseudocount %||% 1)
        # Only warn when transform is active and no pseudocount (zeros will be dropped)
        if (transform == "none" || pseudocount > 0) return(NULL)
        df <- tryCatch(counts_matrix_data(), error=function(e) NULL)
        if (is.null(df)) return(NULL)
        pairs <- counts_plot_pairs()
        if (length(pairs) == 0) return(NULL)
        all_cols <- unique(unlist(lapply(pairs, function(p) c(p$x, p$y))))
        all_cols <- intersect(all_cols, colnames(df))
        n <- sum(apply(df[, all_cols, drop=FALSE], 1, function(r) any(r <= 0, na.rm=TRUE)))
        if (n > 0)
            tags$p(style="color:#856404;font-size:0.82em;margin-top:4px;",
                   paste0("\u26a0 ", n, " row(s) with \u22640 excluded. Set pseudocount\u00a0>\u00a00 to retain them."))
    })

    # Shared plot parameters helper
    counts_plot_params <- reactive({
        list(
            matrix_type         = input$counts_matrix_type %||% "count",
            level               = input$counts_level %||% "sgrna",
            gene_summary        = input$counts_gene_summary %||% "mean",
            gene_min            = as.numeric(input$counts_gene_min %||% 0),
            gene_filter         = input$counts_gene_filter %||% "none",
            transform           = input$counts_transform %||% "none",
            pseudocount         = as.numeric(input$counts_pseudocount %||% 1),
            show_identity  = isTRUE(input$counts_show_identity),
            identity_color = rv$counts_id_color,
            identity_lty   = input$counts_id_lty %||% "dashed",
            identity_lwd   = input$counts_id_lwd %||% 1,
            point_size     = input$counts_point_size %||% 0.6,
            point_alpha    = input$counts_point_alpha %||% 0.35,
            title_fmt      = input$counts_title_fmt %||% "full",
            highlight_config  = get_counts_highlight_config(),
            show_labels       = isTRUE(input$counts_show_labels),
            label_scale       = input$counts_label_scale %||% 1,
            hl_point_size     = input$counts_hl_point_size %||% 1.5,
            show_hist         = isTRUE(input$counts_show_hist) && GGEXTRA_AVAILABLE,
            aspect_ratio1     = isTRUE(input$counts_aspect_ratio1)
        )
    })

    # Processed (summarized + transformed) data for all current plots — used by data download
    counts_processed_data <- reactive({
        df     <- counts_matrix_data()
        pairs  <- counts_plot_pairs()
        req(length(pairs) > 0)
        params <- counts_plot_params()
        type_label <- switch(params$matrix_type, count = "count", norm = "normalized",
                             params$matrix_type)
        lapply(pairs, function(pair) {
            pd <- counts_prepare_pair_data(df, pair$x, pair$y,
                      level        = params$level,
                      gene_summary = params$gene_summary,
                      gene_min     = params$gene_min,
                      gene_filter  = params$gene_filter)
            pd$x <- counts_apply_transform(pd$x, params$transform, params$pseudocount)
            pd$y <- counts_apply_transform(pd$y, params$transform, params$pseudocount)
            x_label <- counts_axis_label(pair$x, type_label, params$transform, params$pseudocount)
            y_label <- counts_axis_label(pair$y, type_label, params$transform, params$pseudocount)
            id_col  <- if (params$level == "sgrna") "sgrna_id" else "gene"
            out <- data.frame(pd$lbl, pd$x, pd$y, stringsAsFactors = FALSE)
            colnames(out) <- c(id_col, x_label, y_label)
            list(name = pair$name, data = out)
        })
    })

    # Returns list of individual ggplot objects (used by both renderPlot and download)
    counts_plot_list <- reactive({
        df <- counts_matrix_data()
        pairs <- counts_plot_pairs()
        req(length(pairs) > 0)
        params <- counts_plot_params()
        lapply(pairs, function(pair) {
            p <- create_single_counts_scatter(df, pair$x, pair$y,
                matrix_type=params$matrix_type,
                level=params$level, gene_summary=params$gene_summary,
                gene_min=params$gene_min, gene_filter=params$gene_filter,
                transform=params$transform, pseudocount=params$pseudocount,
                show_identity=params$show_identity, identity_color=params$identity_color,
                identity_lty=params$identity_lty, identity_lwd=params$identity_lwd,
                plot_name=pair$name, title_fmt=params$title_fmt,
                point_size=params$point_size, point_alpha=params$point_alpha,
                highlight_config=params$highlight_config, config_dir=rv$config_dir,
                show_labels=params$show_labels, label_scale=params$label_scale,
                hl_point_size=params$hl_point_size, aspect_ratio1=params$aspect_ratio1)
            if (params$show_hist) {
                # Move legend to bottom so the right marginal histogram isn't blocked by it
                p <- p + theme(legend.position = "bottom")
                ggExtra::ggMarginal(p, type="histogram", margins="both", size=3,
                                    fill="gray70", color=NA)
            } else p
        })
    })

    # Interactive multi-plot: subplot of individual ggplotly objects
    create_counts_plotly <- reactive({
        df <- counts_matrix_data()
        pairs <- counts_plot_pairs()
        req(length(pairs) > 0)
        params <- counts_plot_params()
        ncols  <- max(1L, as.integer(input$counts_ncols %||% 2))
        nrows  <- ceiling(length(pairs) / ncols)

        plots <- lapply(pairs, function(pair) {
            p <- create_single_counts_scatter(df, pair$x, pair$y,
                     matrix_type=params$matrix_type,
                     level=params$level, gene_summary=params$gene_summary,
                     gene_min=params$gene_min, gene_filter=params$gene_filter,
                     transform=params$transform, pseudocount=params$pseudocount,
                     show_identity=params$show_identity, identity_color=params$identity_color,
                     identity_lty=params$identity_lty, identity_lwd=params$identity_lwd,
                     point_size=params$point_size, point_alpha=params$point_alpha,
                     plot_name=pair$name, title_fmt=params$title_fmt,
                     highlight_config=params$highlight_config, config_dir=rv$config_dir,
                     show_labels=FALSE, label_scale=params$label_scale,
                     hl_point_size=params$hl_point_size, aspect_ratio1=FALSE)
            ggplotly(p, tooltip="text")
        })

        if (length(plots) == 1) {
            fig <- plots[[1]]
            fig <- do.call(plotly::layout, c(list(fig), setNames(
                list(list(scaleanchor="yaxis", scaleratio=1, constrain="domain"),
                     list(constrain="domain")),
                c("xaxis", "yaxis"))))
        } else {
            fig <- plotly::subplot(plots, nrows=nrows, shareX=FALSE, shareY=FALSE,
                                   titleX=TRUE, titleY=TRUE, margin=0.06)
            for (i in seq_len(length(plots))) {
                xa <- if (i == 1) "xaxis" else paste0("xaxis", i)
                ya <- if (i == 1) "yaxis" else paste0("yaxis", i)
                axis_args <- setNames(
                    list(list(scaleanchor=ya, scaleratio=1, constrain="domain"),
                         list(constrain="domain")),
                    c(xa, ya)
                )
                fig <- do.call(plotly::layout, c(list(fig), axis_args))
            }
        }
        fig %>% plotly::layout(hovermode = "closest")
    })

    # Dynamic height: use number of plot slots, not data-dependent pairs,
    # so req() failures in the data chain never corrupt the graphics device.
    counts_static_height <- reactive({
        n_plots <- max(1L, length(rv$counts_plot_ids))
        ncols   <- max(1L, as.integer(input$counts_ncols %||% 2))
        nrows   <- ceiling(n_plots / ncols)
        px_per  <- max(150L, as.integer(input$counts_plot_px %||% 350))
        nrows * px_per
    })

    output$counts_plotly_ui <- renderUI({
        n_plots   <- max(1L, length(rv$counts_plot_ids))
        ncols     <- max(1L, as.integer(input$counts_ncols %||% 2))
        eff_ncols <- min(ncols, n_plots)
        px_h      <- max(150L, as.integer(input$counts_plot_px %||% 350))
        px_w      <- max(150L, as.integer(input$counts_plot_w  %||% 350))
        nrows     <- ceiling(n_plots / eff_ncols)
        max_w     <- paste0(eff_ncols * px_w, "px")
        tags$div(
            style = paste0("max-width:", max_w, "; width:100%;"),
            plotly::plotlyOutput("counts_scatter_plotly", height=paste0(nrows * px_h, "px"))
        )
    })
    output$counts_static_ui <- renderUI({
        n_plots   <- max(1L, length(rv$counts_plot_ids))
        ncols     <- max(1L, as.integer(input$counts_ncols %||% 2))
        eff_ncols <- min(ncols, n_plots)
        px_w      <- max(150L, as.integer(input$counts_plot_w %||% 350))
        max_w     <- paste0(eff_ncols * px_w, "px")
        h         <- paste0(counts_static_height(), "px")
        tags$div(
            style = paste0("max-width:", max_w, "; width:100%;"),
            plotOutput("counts_scatter_static", height=h, width="100%")
        )
    })

    if (PLOTLY_AVAILABLE) {
        output$counts_scatter_plotly <- plotly::renderPlotly({ create_counts_plotly() })
    }
    output$counts_scatter_static <- renderPlot({
        plots <- counts_plot_list()
        ncols <- max(1L, as.integer(input$counts_ncols %||% 2))
        draw_plot_grid(plots, ncols)
    }, height=function() counts_static_height())

    output$download_counts_scatter <- downloadHandler(
        filename = function() {
            ext <- if (tolower(input$counts_dl_fmt) == "pdf") "pdf" else "png"
            paste0(input$counts_matrix_type %||% "count",
                   "_scatter_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
        },
        content = function(file) {
            plots <- isolate(counts_plot_list())
            ncols <- max(1L, as.integer(isolate(input$counts_ncols) %||% 2))
            ext <- tolower(input$counts_dl_fmt)
            w <- input$counts_dl_w %||% 7; h <- input$counts_dl_h %||% 7
            dpi <- as.numeric(input$counts_dl_dpi %||% 300)
            if (ext == "pdf") {
                if (capabilities("cairo")) grDevices::cairo_pdf(file, width=w, height=h)
                else grDevices::pdf(file, width=w, height=h, useDingbats=FALSE)
            } else {
                grDevices::png(file, width=w*dpi, height=h*dpi, res=dpi,
                               type=if (capabilities("cairo")) "cairo" else "windows")
            }
            draw_plot_grid(plots, ncols)
            grDevices::dev.off()
        },
        contentType = "application/octet-stream"
    )

    output$download_counts_tsv <- downloadHandler(
        filename = function() {
            paste0(input$counts_matrix_type %||% "count",
                   "_counts_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
        },
        content = function(file) {
            processed <- isolate(counts_processed_data())
            combined  <- do.call(rbind, lapply(processed, function(p) {
                cbind(plot_name = p$name, p$data, stringsAsFactors = FALSE)
            }))
            write.table(combined, file, sep = "\t", row.names = FALSE, quote = FALSE)
        },
        contentType = "text/tab-separated-values"
    )

    if (OPENXLSX2_AVAILABLE || OPENXLSX_AVAILABLE) {
        output$download_counts_xlsx <- downloadHandler(
            filename = function() {
                paste0(input$counts_matrix_type %||% "count",
                       "_counts_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
            },
            content = function(file) {
                processed <- isolate(counts_processed_data())
                params    <- isolate(counts_plot_params())
                settings_df <- data.frame(
                    Setting = c("Matrix type", "Level", "Transform", "Pseudocount",
                                "Gene summary", "Min count filter", "Filter mode"),
                    Value   = c(params$matrix_type, params$level, params$transform,
                                as.character(params$pseudocount), params$gene_summary,
                                as.character(params$gene_min), params$gene_filter),
                    stringsAsFactors = FALSE
                )
                if (OPENXLSX2_AVAILABLE) {
                    wb <- openxlsx2::wb_workbook()
                    for (p in processed) {
                        sn <- substr(gsub("[^A-Za-z0-9 _-]", "", p$name), 1, 31)
                        if (nchar(sn) == 0) sn <- "Plot"
                        wb <- openxlsx2::wb_add_worksheet(wb, sn)
                        wb <- openxlsx2::wb_add_data(wb, sn, p$data)
                    }
                    wb <- openxlsx2::wb_add_worksheet(wb, "Settings")
                    wb <- openxlsx2::wb_add_data(wb, "Settings", settings_df)
                    openxlsx2::wb_save(wb, file)
                } else {
                    wb <- openxlsx::createWorkbook()
                    for (p in processed) {
                        sn <- substr(gsub("[^A-Za-z0-9 _-]", "", p$name), 1, 31)
                        if (nchar(sn) == 0) sn <- "Plot"
                        openxlsx::addWorksheet(wb, sn)
                        openxlsx::writeData(wb, sn, p$data)
                    }
                    openxlsx::addWorksheet(wb, "Settings")
                    openxlsx::writeData(wb, "Settings", settings_df)
                    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
                }
            },
            contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    }

    create_interactive_volcano <- reactive({
        req(rv$gene_data, input$selection_type)
        rank_method <- input$rank_method; reverse_rank <- input$reverse_rank_order
        result_info <- rv$result_files[[input$comparison]]
        fdr1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05)
        fdr2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
        create_volcano_plot_interactive(rv$gene_data, result_info$comparison_name, input$selection_type,
            get_highlight_config(), rv$config_dir, fdr1, fdr2, input$show_labels, input$label_scale,
            text_scale=input$text_scale, title_scale=input$title_scale,
            show_fdr_line_1=input$show_fdr_line_1, show_fdr_line_2=input$show_fdr_line_2,
            custom_title=input$volcano_title, point_size=input$point_size,
            fdr_color_1=get_fdr_color_1(), fdr_color_2=get_fdr_color_2(), plot_metric=input$plot_metric)
    })

    create_interactive_ranked <- reactive({
        req(rv$gene_data, input$selection_type)
        rank_method <- input$rank_method; reverse_rank <- input$reverse_rank_order
        result_info <- rv$result_files[[input$comparison]]
        fdr1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05)
        fdr2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
        create_ranked_plot_interactive(rv$gene_data, result_info$comparison_name, input$selection_type,
            get_highlight_config(), rv$config_dir, fdr1, fdr2, input$show_labels, input$label_scale,
            scale_point_size=input$scale_point_size, color_by_fdr=input$color_by_fdr,
            text_scale=input$text_scale, title_scale=input$title_scale,
            custom_title=input$ranked_title, point_size=input$point_size,
            fdr_color_1=get_fdr_color_1(), fdr_color_2=get_fdr_color_2(),
            rank_method=input$rank_method, reverse_rank=input$reverse_rank_order,
            sig_use_plot_thresholds=if(is.null(input$use_plot_thresholds)) TRUE else input$use_plot_thresholds,
            plot_metric=input$plot_metric)
    })

    output$volcano_plot <- {
        if (PLOTLY_AVAILABLE) plotly::renderPlotly({ render_plot_with_fallback(create_interactive_volcano(), use_plotly=if(is.null(input$use_interactive_plots)) TRUE else input$use_interactive_plots, tooltip="text") })
        else renderPlot({ create_interactive_volcano() })
    }
    output$volcano_plot_static <- renderPlot({ create_interactive_volcano() })
    output$ranked_plot <- {
        if (PLOTLY_AVAILABLE) plotly::renderPlotly({ render_plot_with_fallback(create_interactive_ranked(), use_plotly=if(is.null(input$use_interactive_plots)) TRUE else input$use_interactive_plots, tooltip="text") })
        else renderPlot({ create_interactive_ranked() })
    }
    output$ranked_plot_static <- renderPlot({ create_interactive_ranked() })

    save_plot <- function(file, plot, ext) {
        if (ext == "html") {
            if (PLOTLY_AVAILABLE) { htmlwidgets::saveWidget(ggplotly(plot, tooltip="text") %>% layout(hovermode="closest"), file=file); return() }
            ext <- "png"
        }
        if (ext == "pdf") {
            if (capabilities("cairo")) ggsave(filename=file, plot=plot, width=input$download_width, height=input$download_height, units="in", device=grDevices::cairo_pdf, bg="white")
            else ggsave(filename=file, plot=plot, width=input$download_width, height=input$download_height, units="in", device=function(...) grDevices::pdf(..., useDingbats=FALSE), bg="white")
        } else {
            ggsave(filename=file, plot=plot, width=input$download_width, height=input$download_height, units="in", device="png", dpi=as.numeric(input$download_dpi), bg="white")
        }
    }

    output$download_volcano <- downloadHandler(
        filename = function() { ext <- if (tolower(input$download_format) == "pdf") "pdf" else "png"; paste0(input$comparison %||% "comparison", "_", input$selection_type %||% "neg", "_volcano_", format(Sys.time(),"%Y%m%d_%H%M%S"), ".", ext) },
        content = function(file) {
            req(rv$gene_data, input$selection_type, input$comparison)
            result_info <- rv$result_files[[input$comparison]]
            fdr1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05); fdr2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
            plot <- create_volcano_plot_interactive(rv$gene_data, result_info$comparison_name, input$selection_type,
                get_highlight_config(), rv$config_dir, fdr1, fdr2, input$show_labels, input$label_scale,
                text_scale=input$text_scale, title_scale=input$title_scale,
                show_fdr_line_1=input$show_fdr_line_1, show_fdr_line_2=input$show_fdr_line_2,
                custom_title=input$volcano_title, point_size=input$point_size, fdr_color_1=get_fdr_color_1(), fdr_color_2=get_fdr_color_2(), plot_metric=input$plot_metric)
            save_plot(file, plot, tolower(input$download_format))
        },
        contentType = "application/octet-stream"
    )

    output$download_ranked <- downloadHandler(
        filename = function() { ext <- if (tolower(input$download_format) == "pdf") "pdf" else "png"; paste0(input$comparison %||% "comparison", "_", input$selection_type %||% "neg", "_ranked_", format(Sys.time(),"%Y%m%d_%H%M%S"), ".", ext) },
        content = function(file) {
            req(rv$gene_data, input$selection_type, input$comparison)
            result_info <- rv$result_files[[input$comparison]]
            fdr1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05); fdr2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
            plot <- create_ranked_plot_interactive(rv$gene_data, result_info$comparison_name, input$selection_type,
                get_highlight_config(), rv$config_dir, fdr1, fdr2, input$show_labels, input$label_scale,
                scale_point_size=input$scale_point_size, color_by_fdr=input$color_by_fdr,
                text_scale=input$text_scale, title_scale=input$title_scale,
                custom_title=input$ranked_title, point_size=input$point_size, fdr_color_1=get_fdr_color_1(), fdr_color_2=get_fdr_color_2(),
                rank_method=input$rank_method, reverse_rank=input$reverse_rank_order,
                sig_use_plot_thresholds=if(is.null(input$use_plot_thresholds)) TRUE else input$use_plot_thresholds, plot_metric=input$plot_metric)
            save_plot(file, plot, tolower(input$download_format))
        },
        contentType = "application/octet-stream"
    )
}

# ── Plot Functions (unchanged from v1.0) ─────────────────────────────────────

create_volcano_plot_interactive <- function(gene_data, comparison_name, result_type,
                                            highlight_config, config_dir, fdr_thresh1, fdr_thresh2,
                                            show_labels, label_scale, text_scale = 1, title_scale = 1,
                                            show_fdr_line_1 = TRUE, show_fdr_line_2 = TRUE,
                                            custom_title = "", point_size = 2,
                                            fdr_color_1 = "#FF0000", fdr_color_2 = "#FFA500",
                                            plot_metric = "lfc") {
    if (result_type == "neg") { lfc_col <- "neg_lfc"; fdr_col <- "neg_fdr"; title_suffix <- "Negative Selection" }
    else { lfc_col <- "pos_lfc"; fdr_col <- "pos_fdr"; title_suffix <- "Positive Selection" }
    score_col <- if (result_type == "neg") "neg_score" else "pos_score"
    colnames(gene_data) <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(gene_data)))
    lfc_col <- tolower(lfc_col); fdr_col <- tolower(fdr_col); score_col <- tolower(score_col)
    gene_id_col <- if ("id" %in% colnames(gene_data)) "id" else if ("gene" %in% colnames(gene_data)) "gene" else colnames(gene_data)[1]
    plot_data <- gene_data %>% filter(!is.na(.data[[lfc_col]]) & !is.na(.data[[fdr_col]]))
    min_nonzero_fdr <- min(plot_data[[fdr_col]][plot_data[[fdr_col]] > 0], na.rm = TRUE)
    replacement_fdr <- min(1e-200, min_nonzero_fdr / 10)
    plot_data <- plot_data %>% mutate(adjusted_fdr = ifelse(.data[[fdr_col]] == 0, replacement_fdr, .data[[fdr_col]]),
                                       neg_log10_fdr = -log10(adjusted_fdr), highlight_group = "none", highlight_color = "gray50")
    if (plot_metric == "score") {
        min_nonzero_score <- min(plot_data[[score_col]][plot_data[[score_col]] > 0], na.rm = TRUE)
        replacement_score <- min(1e-200, min_nonzero_score / 10)
        plot_data <- plot_data %>% mutate(score_value = -log10(ifelse(.data[[score_col]] == 0, replacement_score, .data[[score_col]])))
        x_col <- "score_value"; x_label <- "-log10(Score)"; metric_label <- "-log10(Score)"
    } else { x_col <- lfc_col; x_label <- "Log2 Fold Change"; metric_label <- "LFC" }
    genes_to_label <- NULL; legend_data <- NULL
    if (!is.null(highlight_config) && length(highlight_config) > 0) {
        for (i in seq_along(highlight_config)) {
            group_info <- highlight_config[[i]]
            group_color <- group_info$color %||% "black"; group_name <- group_info$group_name %||% paste0("Group ", i)
            if (!is.null(group_info$is_custom) && group_info$is_custom) highlight_genes <- group_info$genes
            else {
                gene_file <- group_info$file
                if (!is.null(config_dir) && !file.exists(gene_file)) gene_file <- file.path(config_dir, gene_file)
                highlight_genes <- load_gene_list(gene_file)
            }
            if (length(highlight_genes) > 0)
                plot_data <- plot_data %>% mutate(
                    highlight_group = ifelse(.data[[gene_id_col]] %in% highlight_genes & highlight_group == "none", group_name, highlight_group),
                    highlight_color = ifelse(.data[[gene_id_col]] %in% highlight_genes & highlight_color == "gray50", group_color, highlight_color))
        }
        genes_to_label <- plot_data %>% filter(highlight_group != "none")
        legend_data <- tibble(group = sapply(highlight_config, function(x) x$group_name %||% "Unknown"), color = sapply(highlight_config, function(x) x$color %||% "black"))
    }
    if (plot_metric == "score") x_limit <- ceiling(max(plot_data[[x_col]], na.rm=TRUE) * 10) / 10
    else x_limit <- ceiling(max(abs(plot_data[[x_col]]), na.rm=TRUE) * 10) / 10
    fdr_legend_data <- tibble()
    if (show_fdr_line_1) fdr_legend_data <- bind_rows(fdr_legend_data, tibble(label=paste0("FDR = ", fdr_thresh1), color=fdr_color_1))
    if (show_fdr_line_2) fdr_legend_data <- bind_rows(fdr_legend_data, tibble(label=paste0("FDR = ", fdr_thresh2), color=fdr_color_2))
    p <- ggplot(plot_data, aes(x = .data[[x_col]], y = neg_log10_fdr))
    if (plot_metric == "lfc") p <- p + geom_vline(xintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7)
    if (show_fdr_line_1) p <- p + geom_hline(yintercept = -log10(fdr_thresh1), linetype = "dashed", color = fdr_color_1, linewidth = 0.5)
    if (show_fdr_line_2) p <- p + geom_hline(yintercept = -log10(fdr_thresh2), linetype = "dashed", color = fdr_color_2, linewidth = 0.5)
    p <- p +
        geom_point(data=plot_data %>% filter(highlight_group=="none"),
                   aes(text=paste0(.data[[gene_id_col]],"<br>",metric_label,": ",round(.data[[x_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))),
                   color="gray50", size=point_size*0.75, alpha=0.5) +
        geom_point(data=plot_data %>% filter(highlight_group!="none"),
                   aes(color=highlight_group, text=paste0(.data[[gene_id_col]],"<br>",metric_label,": ",round(.data[[x_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))),
                   size=point_size, alpha=0.8)
    if (!is.null(legend_data) && !is.null(genes_to_label) && nrow(genes_to_label) > 0) {
        color_values <- setNames(legend_data$color, legend_data$group)
        p <- p + scale_color_manual(name="Gene Group", values=color_values) + guides(color=guide_legend(override.aes=list(size=4)))
        if (show_labels && !is.null(genes_to_label) && nrow(genes_to_label) > 0)
            p <- p + geom_text_repel(data=genes_to_label, aes(label=.data[[gene_id_col]], color=highlight_group),
                                     size=3*label_scale, max.overlaps=Inf, box.padding=0.8, point.padding=0.8, min.segment.length=0, force=3, force_pull=1, fontface="bold", show.legend=FALSE)
    }
    if (nrow(fdr_legend_data) > 0) {
        y_positions <- c(); if (show_fdr_line_1) y_positions <- c(y_positions, -log10(fdr_thresh1)); if (show_fdr_line_2) y_positions <- c(y_positions, -log10(fdr_thresh2))
        fdr_annotations <- tibble(x=rep(x_limit*0.95, length(y_positions)), y=y_positions, label=fdr_legend_data$label, color=fdr_legend_data$color)
        p <- p + geom_text(data=fdr_annotations, aes(x=x, y=y, label=label), color=fdr_annotations$color, hjust=1, vjust=-0.5, size=3*text_scale, fontface="bold")
    }
    plot_title <- if (!is.null(custom_title) && nzchar(custom_title)) custom_title else paste0(gsub(":", " vs. ", gsub("_", " ", comparison_name)), " (", title_suffix, ")")
    p + labs(title=plot_title, x=x_label, y="-log10(FDR)") +
        (if (plot_metric=="score") scale_x_continuous(limits=c(0,x_limit)) else scale_x_continuous(limits=c(-x_limit,x_limit))) +
        theme_minimal() + theme(plot.title=element_text(hjust=0.5,size=14*text_scale*title_scale), axis.title=element_text(size=11*text_scale),
                                axis.text=element_text(size=10*text_scale), legend.title=element_text(size=11*text_scale),
                                legend.text=element_text(size=10*text_scale), plot.background=element_rect(fill="white",color=NA),
                                panel.background=element_rect(fill="white",color=NA), legend.position="right")
}

create_ranked_plot_interactive <- function(gene_data, comparison_name, result_type,
                                           highlight_config, config_dir, fdr_thresh1, fdr_thresh2,
                                           show_labels, label_scale, scale_point_size = TRUE, color_by_fdr = TRUE,
                                           text_scale = 1, title_scale = 1, custom_title = "", point_size = 2,
                                           fdr_color_1 = "#FF0000", fdr_color_2 = "#FFA500",
                                           sig_use_plot_thresholds = TRUE, rank_method = "lfc", reverse_rank = FALSE,
                                           plot_metric = "lfc") {
    if (result_type == "neg") { lfc_col <- "neg_lfc"; fdr_col <- "neg_fdr"; title_suffix <- "Negative Selection" }
    else { lfc_col <- "pos_lfc"; fdr_col <- "pos_fdr"; title_suffix <- "Positive Selection" }
    score_col <- if (result_type == "neg") "neg_score" else "pos_score"
    colnames(gene_data) <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(gene_data)))
    lfc_col <- tolower(lfc_col); fdr_col <- tolower(fdr_col); score_col <- tolower(score_col)
    rank_col <- tolower(if (result_type == "neg") "neg_rank" else "pos_rank")
    gene_id_col <- if ("id" %in% colnames(gene_data)) "id" else if ("gene" %in% colnames(gene_data)) "gene" else colnames(gene_data)[1]
    plot_data <- gene_data %>% filter(!is.na(.data[[lfc_col]]) & !is.na(.data[[fdr_col]])) %>%
        mutate(significance = case_when(.data[[fdr_col]] < fdr_thresh1 ~ fdr_color_1,
                                        .data[[fdr_col]] >= fdr_thresh1 & .data[[fdr_col]] < fdr_thresh2 ~ fdr_color_2, TRUE ~ "black"),
               point_size = case_when(.data[[fdr_col]] < 0.0001 ~ 4, .data[[fdr_col]] < 0.001 ~ 3, .data[[fdr_col]] < 0.01 ~ 2.5,
                                      .data[[fdr_col]] < fdr_thresh1 ~ 2, .data[[fdr_col]] < fdr_thresh2 ~ 1.5, TRUE ~ 1),
               highlight_group = "none", highlight_color = NA_character_)
    if (plot_metric == "score") {
        min_nonzero_score <- min(plot_data[[score_col]][plot_data[[score_col]] > 0], na.rm = TRUE)
        replacement_score <- min(1e-200, min_nonzero_score / 10)
        plot_data <- plot_data %>% mutate(score_value = -log10(ifelse(.data[[score_col]] == 0, replacement_score, .data[[score_col]])))
        y_col <- "score_value"; y_label <- "-log10(Score)"; metric_label <- "-log10(Score)"
    } else { y_col <- lfc_col; y_label <- "Log2 Fold Change"; metric_label <- "LFC" }
    if (rank_method == "rra" && rank_col %in% colnames(plot_data)) plot_data <- plot_data %>% arrange(.data[[rank_col]])
    else plot_data <- plot_data %>% arrange(.data[[lfc_col]])
    plot_data <- plot_data %>% mutate(rank = row_number())
    genes_to_label <- NULL; legend_data <- NULL
    if (!is.null(highlight_config) && length(highlight_config) > 0) {
        for (i in seq_along(highlight_config)) {
            group_info <- highlight_config[[i]]; group_color <- group_info$color %||% "black"; group_name <- group_info$group_name %||% paste0("Group ", i)
            if (!is.null(group_info$is_custom) && group_info$is_custom) highlight_genes <- group_info$genes
            else { gene_file <- group_info$file; if (!is.null(config_dir) && !file.exists(gene_file)) gene_file <- file.path(config_dir, gene_file); highlight_genes <- load_gene_list(gene_file) }
            if (length(highlight_genes) > 0)
                plot_data <- plot_data %>% mutate(
                    highlight_group = ifelse(.data[[gene_id_col]] %in% highlight_genes & highlight_group == "none", group_name, highlight_group),
                    highlight_color = ifelse(.data[[gene_id_col]] %in% highlight_genes & is.na(highlight_color), group_color, highlight_color))
        }
        genes_to_label <- plot_data %>% filter(highlight_group != "none")
        legend_data <- tibble(group=sapply(highlight_config, function(x) x$group_name %||% "Unknown"), color=sapply(highlight_config, function(x) x$color %||% "black"))
    }
    plot_data <- plot_data %>% mutate(
        fdr_category = factor(case_when(significance==fdr_color_1~paste0("FDR < ",fdr_thresh1), significance==fdr_color_2~paste0(fdr_thresh1," \u2264 FDR < ",fdr_thresh2), TRUE~paste0("FDR \u2265 ",fdr_thresh2)),
                              levels=c(paste0("FDR < ",fdr_thresh1),paste0(fdr_thresh1," \u2264 FDR < ",fdr_thresh2),paste0("FDR \u2265 ",fdr_thresh2))),
        fdr_bin = factor(case_when(.data[[fdr_col]]<1e-5~"<1e-5",.data[[fdr_col]]<1e-4~"<1e-4",.data[[fdr_col]]<1e-3~"< 0.001",.data[[fdr_col]]<1e-2~"< 0.01",
                                   .data[[fdr_col]]<fdr_thresh1~paste0("< ",fdr_thresh1),.data[[fdr_col]]<fdr_thresh2~paste0("< ",fdr_thresh2),TRUE~paste0("\u2265 ",fdr_thresh2)),
                         levels=c("<1e-5","<1e-4","< 0.001","< 0.01",paste0("< ",fdr_thresh1),paste0("< ",fdr_thresh2),paste0("\u2265 ",fdr_thresh2))))
    size_values <- setNames(c(4.5,4,3.5,3,2.5,2,1.5), c("<1e-5","<1e-4","< 0.001","< 0.01",paste0("< ",fdr_thresh1),paste0("< ",fdr_thresh2),paste0("\u2265 ",fdr_thresh2)))
    highlight_colors <- if (!is.null(highlight_config)) sapply(highlight_config, function(x) x$color %||% "black") else character(0)
    highlight_names  <- if (!is.null(highlight_config)) sapply(highlight_config, function(x) x$group_name %||% "Unknown") else character(0)
    zero_line <- if (plot_metric == "lfc") geom_hline(yintercept=0, color="black", linewidth=0.8, alpha=0.7) else NULL
    if (color_by_fdr) {
        # Two layers: background points (fill only, no outline) + highlighted points (fill + colored outline).
        # Keeping layers separate avoids the col=NA rendering issue where shape=21 points with an NA
        # stroke color are suppressed by some graphics devices.
        fdr_cat_levels  <- levels(plot_data$fdr_category)
        fdr_fill_values <- setNames(c(fdr_color_1, fdr_color_2, "black"), fdr_cat_levels)
        scaled_size_values <- size_values * (point_size / 2)
        plot_data_bg <- plot_data %>% filter(highlight_group == "none")
        plot_data_hi <- plot_data %>% filter(highlight_group != "none")
        if (scale_point_size) {
            p <- ggplot(plot_data, aes(x=rank, y=.data[[y_col]])) + zero_line +
                geom_point(data=plot_data_bg,
                           aes(fill=fdr_category, size=fdr_bin,
                               text=paste0(.data[[gene_id_col]],"<br>Rank: ",rank,"<br>",metric_label,": ",round(.data[[y_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))),
                           shape=21, stroke=0, alpha=0.8)
            if (nrow(plot_data_hi) > 0)
                p <- p + geom_point(data=plot_data_hi,
                                    aes(fill=fdr_category, color=highlight_group, size=fdr_bin,
                                        text=paste0(.data[[gene_id_col]],"<br>Rank: ",rank,"<br>",metric_label,": ",round(.data[[y_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))),
                                    shape=21, stroke=1.0, alpha=0.8)
            p <- p + scale_size_manual(name="FDR (size)", values=scaled_size_values,
                                       guide=guide_legend(override.aes=list(fill="black", color="black")))
        } else {
            p <- ggplot(plot_data, aes(x=rank, y=.data[[y_col]])) + zero_line +
                geom_point(data=plot_data_bg,
                           aes(fill=fdr_category,
                               text=paste0(.data[[gene_id_col]],"<br>Rank: ",rank,"<br>",metric_label,": ",round(.data[[y_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))),
                           shape=21, stroke=0, size=point_size, alpha=0.8)
            if (nrow(plot_data_hi) > 0)
                p <- p + geom_point(data=plot_data_hi,
                                    aes(fill=fdr_category, color=highlight_group,
                                        text=paste0(.data[[gene_id_col]],"<br>Rank: ",rank,"<br>",metric_label,": ",round(.data[[y_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))),
                                    shape=21, stroke=1.0, size=point_size, alpha=0.8)
        }
        p <- p +
            scale_fill_manual(name="Significance", values=fdr_fill_values,
                               guide=guide_legend(override.aes=list(size=4, color=NA)))
        if (!is.null(genes_to_label) && nrow(genes_to_label) > 0 && length(highlight_names) > 0)
            p <- p + scale_color_manual(name="Gene Group",
                               values=setNames(highlight_colors, highlight_names),
                               breaks=highlight_names,
                               guide=guide_legend(override.aes=list(size=4, fill="white")))
    } else {
        if (scale_point_size) {
            p <- ggplot(plot_data, aes(x=rank,y=.data[[y_col]])) + zero_line +
                geom_point(data=plot_data %>% filter(highlight_group=="none"), aes(size=fdr_bin,text=paste0(.data[[gene_id_col]],"<br>Rank: ",rank,"<br>",metric_label,": ",round(.data[[y_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))), color="gray60", alpha=0.6) +
                geom_point(data=plot_data %>% filter(highlight_group!="none"), aes(color=highlight_group,size=fdr_bin,text=paste0(.data[[gene_id_col]],"<br>Rank: ",rank,"<br>",metric_label,": ",round(.data[[y_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))), alpha=0.8) +
                scale_size_manual(name="FDR (size)",values=size_values*(point_size/2))
        } else {
            p <- ggplot(plot_data, aes(x=rank,y=.data[[y_col]])) + zero_line +
                geom_point(data=plot_data %>% filter(highlight_group=="none"), aes(text=paste0(.data[[gene_id_col]],"<br>Rank: ",rank,"<br>",metric_label,": ",round(.data[[y_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))), color="gray60", size=point_size, alpha=0.6) +
                geom_point(data=plot_data %>% filter(highlight_group!="none"), aes(color=highlight_group,text=paste0(.data[[gene_id_col]],"<br>Rank: ",rank,"<br>",metric_label,": ",round(.data[[y_col]],3),"<br>FDR: ",round(.data[[fdr_col]],5))), size=point_size, alpha=0.8)
        }
        p <- p + scale_color_manual(name="Gene Group",values=c(setNames(highlight_colors,highlight_names),"none"="gray60"),breaks=if(length(highlight_names)>0) highlight_names else NULL) +
            guides(color=guide_legend(override.aes=list(size=4)))
    }
    if (show_labels && !is.null(genes_to_label) && nrow(genes_to_label) > 0)
        p <- p + geom_text_repel(data=genes_to_label, aes(label=.data[[gene_id_col]]),
                                  color=ifelse(!is.na(genes_to_label$highlight_color), genes_to_label$highlight_color, "black"),
                                  size=3*label_scale, max.overlaps=Inf, box.padding=0.8, point.padding=0.8, min.segment.length=0, force=3, force_pull=1, fontface="bold")
    plot_title <- if (!is.null(custom_title) && nzchar(custom_title)) custom_title else paste0(gsub(":", " vs. ", gsub("_", " ", comparison_name)), " (", title_suffix, ")")
    p <- p + labs(title=plot_title, x=if(rank_method=="rra") "Rank (RRA)" else "Rank (LFC)", y=y_label) +
        theme_minimal() + theme(plot.title=element_text(hjust=0.5,size=14*text_scale*title_scale), axis.title=element_text(size=11*text_scale),
                                axis.text=element_text(size=10*text_scale), legend.title=element_text(size=11*text_scale), legend.text=element_text(size=10*text_scale),
                                panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
                                plot.background=element_rect(fill="white",color=NA), panel.background=element_rect(fill="white",color=NA),
                                legend.position="right", legend.box="vertical")
    if (reverse_rank) p <- p + scale_x_reverse()
    if (is.null(legend_data)) p <- p + guides(color="none")
    return(p)
}

# ── Counts Scatterplot ────────────────────────────────────────────────────────

# Draw a list of ggplots in a grid layout to the CURRENT open device.
# Uses only base grid + ggplot2::ggplotGrob — no cowplot or extra packages needed.
# cowplot::plot_grid internally opens a null PDF device to measure plots, which
# fails on some Windows systems ("failed to load default encoding / WinAnsi.enc").
draw_plot_grid <- function(plots, ncols) {
    if (length(plots) == 0) return(invisible(NULL))
    if (length(plots) == 1) { print(plots[[1]]); return(invisible(NULL)) }
    grobs <- lapply(plots, function(p) if (inherits(p, "gg")) ggplot2::ggplotGrob(p) else p)
    nc <- min(max(1L, as.integer(ncols)), length(grobs))
    nr <- ceiling(length(grobs) / nc)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr, nc)))
    for (k in seq_along(grobs)) {
        r <- ceiling(k / nc)
        c <- ((k - 1L) %% nc) + 1L
        grid::pushViewport(grid::viewport(layout.pos.row = r, layout.pos.col = c))
        grid::grid.draw(grobs[[k]])
        grid::popViewport()
    }
    grid::popViewport()
    invisible(NULL)
}

# Prepare a (x, y, lbl) data frame for one scatter plot at the requested level.
# For gene level: optionally filter low-count sgRNAs then summarize by target.
counts_prepare_pair_data <- function(matrix_data, x_col, y_col,
                                     level="sgrna",
                                     gene_summary="mean", gene_min=0, gene_filter="none") {
    id_col     <- if ("sgrna_id"     %in% colnames(matrix_data)) "sgrna_id"
                  else colnames(matrix_data)[1]
    target_col <- if ("sgrna_target" %in% colnames(matrix_data)) "sgrna_target"
                  else colnames(matrix_data)[2]

    if (level == "sgrna") {
        return(data.frame(x=matrix_data[[x_col]], y=matrix_data[[y_col]],
                          lbl=matrix_data[[id_col]], gene=matrix_data[[target_col]],
                          stringsAsFactors=FALSE))
    }

    # Gene level: build working frame
    sub <- data.frame(gene=matrix_data[[target_col]],
                      x=matrix_data[[x_col]], y=matrix_data[[y_col]],
                      stringsAsFactors=FALSE)
    sub <- sub[!is.na(sub$x) & !is.na(sub$y), , drop=FALSE]

    # Filter low-count sgRNAs before summarising
    if (gene_min > 0 && gene_filter != "none") {
        keep <- switch(gene_filter,
            either = sub$x >= gene_min | sub$y >= gene_min,
            both   = sub$x >= gene_min & sub$y >= gene_min,
            rep(TRUE, nrow(sub)))
        sub <- sub[keep, , drop=FALSE]
    }

    if (nrow(sub) == 0)
        return(data.frame(x=numeric(0), y=numeric(0), lbl=character(0)))

    # Summarise by gene using tapply
    agg_fn <- switch(gene_summary, mean=mean, median=median, max=max, mean)
    x_agg  <- tapply(sub$x, sub$gene, agg_fn, na.rm=TRUE)
    y_agg  <- tapply(sub$y, sub$gene, agg_fn, na.rm=TRUE)
    genes  <- names(x_agg)
    data.frame(x=as.numeric(x_agg), y=as.numeric(y_agg), lbl=genes,
               stringsAsFactors=FALSE)
}

# Apply a named transform to a numeric vector.
# pseudocount is added before the log; if pseudocount==0, values <=0 become NA.
counts_apply_transform <- function(v, transform, pseudocount) {
    if (transform == "none") return(v)
    vp <- v + pseudocount
    if (pseudocount == 0) vp[v <= 0] <- NA
    switch(transform, log10 = log10(vp), log2 = log2(vp), vp)
}

# Build an axis label that describes the transform applied.
counts_axis_label <- function(col_name, type_label, transform, pseudocount) {
    if (transform == "none") return(paste0(col_name, " (", type_label, ")"))
    base_str <- if (transform == "log10") "log\u2081\u2080" else "log\u2082"
    if (pseudocount == 0)
        paste0(col_name, " (", base_str, " ", type_label, ")")
    else
        paste0(col_name, " (", base_str, "(", type_label, "+", pseudocount, "))")
}

# Compute shared axis limits from already-transformed values.
counts_axis_lim <- function(vals, transform) {
    vals <- vals[!is.na(vals) & is.finite(vals)]
    if (length(vals) == 0) return(c(0, 1))
    r   <- range(vals)
    # Apply consistent padding to both bounds
    pad <- diff(r) * 0.04
    if (pad == 0) pad <- 0.5
    # Ensure lower bound has meaningful padding even when minimum is 0
    lower <- r[1] - pad
    upper <- r[2] + pad
    c(lower, upper)
}

# Single scatter (used for both interactive and static multi-plot)
create_single_counts_scatter <- function(matrix_data, x_col, y_col,
                                          matrix_type="count",
                                          level="sgrna",
                                          gene_summary="mean", gene_min=0, gene_filter="none",
                                          transform="none", pseudocount=1,
                                          show_identity=TRUE, identity_color="#FF0000",
                                          identity_lty="dashed", identity_lwd=1,
                                          point_size=0.6, point_alpha=0.35,
                                          show_stats=TRUE, point_color=NULL,
                                          highlight_config=NULL, config_dir=NULL,
                                          show_labels=TRUE, label_scale=1.0, hl_point_size=NULL,
                                          plot_name=NULL, title_fmt="full",
                                          aspect_ratio1=TRUE) {
    base_df <- counts_prepare_pair_data(matrix_data, x_col, y_col,
                                        level=level, gene_summary=gene_summary,
                                        gene_min=gene_min, gene_filter=gene_filter)
    x_t <- counts_apply_transform(base_df$x, transform, pseudocount)
    y_t <- counts_apply_transform(base_df$y, transform, pseudocount)
    df  <- data.frame(x=x_t, y=y_t, lbl=base_df$lbl, stringsAsFactors=FALSE)
    if ("gene" %in% colnames(base_df)) df$gene <- base_df$gene  # needed for sgRNA-level highlight matching
    df  <- df[!is.na(df$x) & !is.na(df$y) & is.finite(df$x) & is.finite(df$y), , drop=FALSE]
    if (nrow(df) == 0)
        return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data.", size=5) + theme_void())

    # Determine gene ID column
    gene_id_col <- if ("id" %in% colnames(df)) "id" else if ("gene" %in% colnames(df)) "gene" else "lbl"

    # Apply gene highlighting if config provided
    highlight_legend_data <- NULL
    if (!is.null(highlight_config) && length(highlight_config) > 0) {
        df$highlight_group <- "none"
        df$highlight_color <- "gray60"
        for (i in seq_along(highlight_config)) {
            group_info <- highlight_config[[i]]
            group_color <- group_info$color %||% "black"
            group_name <- group_info$group_name %||% paste0("Group ", i)
            if (!is.null(group_info$is_custom) && group_info$is_custom)
                highlight_genes <- group_info$genes
            else {
                gene_file <- group_info$file
                if (!is.null(config_dir) && !file.exists(gene_file))
                    gene_file <- file.path(config_dir, gene_file)
                highlight_genes <- load_gene_list(gene_file)
            }
            if (length(highlight_genes) > 0)
                df <- df %>% mutate(
                    highlight_group = ifelse(.data[[gene_id_col]] %in% highlight_genes & highlight_group == "none", group_name, highlight_group),
                    highlight_color = ifelse(.data[[gene_id_col]] %in% highlight_genes & highlight_color == "gray60", group_color, highlight_color))
        }
        highlight_legend_data <- tibble(group = sapply(highlight_config, function(x) x$group_name %||% "Unknown"),
                                       color = sapply(highlight_config, function(x) x$color %||% "black"))
    } else {
        df$highlight_group <- "none"
        df$highlight_color <- "gray60"
    }

    # Set default point color based on level if not provided (only used when no highlighting)
    if (is.null(point_color))
        point_color <- if (level == "sgrna") "steelblue" else "purple"
    if (is.null(hl_point_size))
        hl_point_size <- max(point_size * 2, 1.2)

    # Calculate correlation coefficients if stats will be shown
    pearson_r <- NA
    spearman_rho <- NA
    if (show_stats && nrow(df) > 1) {
        pearson_r <- cor(df$x, df$y, method="pearson", use="complete.obs")
        spearman_rho <- cor(df$x, df$y, method="spearman", use="complete.obs")
    }

    lim        <- counts_axis_lim(c(df$x, df$y), transform)
    type_label <- if (matrix_type == "cpm") "CPM" else "Counts"
    x_lab      <- counts_axis_label(x_col, type_label, transform, pseudocount)
    y_lab      <- counts_axis_label(y_col, type_label, transform, pseudocount)

    p <- ggplot(df, aes(x=x, y=y,
                        text=paste0(lbl, "<br>", x_col, ": ", round(x, 2),
                                    "<br>", y_col, ": ", round(y, 2))))

    # Render points: two layers if highlighting is active, single layer otherwise
    genes_to_label <- NULL
    if (!is.null(highlight_config) && length(highlight_config) > 0) {
        df_bg <- df %>% filter(highlight_group == "none")
        df_hi <- df %>% filter(highlight_group != "none")
        p <- p + geom_point(data=df_bg, size=point_size, alpha=point_alpha*0.6, color="gray60")
        if (nrow(df_hi) > 0) {
            p <- p + geom_point(data=df_hi, aes(color=highlight_group), size=hl_point_size, alpha=point_alpha)
            genes_to_label <- df_hi
        }
    } else {
        p <- p + geom_point(size=point_size, alpha=point_alpha, color=point_color)
    }

    if (show_identity)
        p <- p + geom_abline(slope=1, intercept=0, color=identity_color,
                             linetype=identity_lty, linewidth=identity_lwd)
    p <- p + scale_x_continuous(limits=lim) + scale_y_continuous(limits=lim)

    # Add legend for highlight groups if present
    if (!is.null(highlight_legend_data) && nrow(highlight_legend_data) > 0 &&
            !is.null(genes_to_label) && nrow(genes_to_label) > 0) {
        color_values <- setNames(highlight_legend_data$color, highlight_legend_data$group)
        p <- p + scale_color_manual(name="Gene Group", values=color_values) +
                 guides(color=guide_legend(override.aes=list(size=3)))
    }

    # Gene labels via ggrepel (static mode only — caller passes show_labels=FALSE for plotly)
    if (show_labels && !is.null(genes_to_label) && nrow(genes_to_label) > 0)
        p <- p + ggrepel::geom_text_repel(
            data=genes_to_label, aes(label=lbl, color=highlight_group),
            size=3*label_scale, max.overlaps=Inf,
            box.padding=0.8, point.padding=0.8,
            min.segment.length=0, force=3, force_pull=1,
            fontface="bold", show.legend=FALSE)

    has_name  <- !is.null(plot_name) && nzchar(plot_name)
    x_disp    <- gsub("_", " ", x_col)
    y_disp    <- gsub("_", " ", y_col)
    title_str <- switch(title_fmt,
        full = paste0(x_disp, " vs ", y_disp, if (has_name) paste0(" (", plot_name, ")") else ""),
        xy   = paste0(x_disp, " vs ", y_disp),
        name = if (has_name) plot_name else paste0(x_disp, " vs ", y_disp),
               paste0(x_disp, " vs ", y_disp)
    )

    # Build subtitle with optional correlation statistics
    subtitle_parts <- c()
    if (show_stats) {
        subtitle_parts <- c(paste0("n = ", nrow(df)))
        if (!is.na(pearson_r))
            subtitle_parts <- c(subtitle_parts, paste0("r = ", round(pearson_r, 3)))
        if (!is.na(spearman_rho))
            subtitle_parts <- c(subtitle_parts, paste0("rho = ", round(spearman_rho, 3)))
    }
    subtitle_str <- paste(subtitle_parts, collapse=" | ")

    p + labs(title=title_str, subtitle=subtitle_str, x=x_lab, y=y_lab) +
        theme_minimal() +
        theme(plot.title=element_text(hjust=0.5, size=11),
              plot.subtitle=element_text(hjust=0.5, size=9),
              axis.title=element_text(size=9), axis.text=element_text(size=8),
              plot.background=element_rect(fill="white", color=NA),
              panel.background=element_rect(fill="white", color=NA),
              legend.position=if (!is.null(highlight_legend_data) && nrow(highlight_legend_data) > 0) "right" else "none",
              aspect.ratio=if (aspect_ratio1) 1 else NULL)
}

shinyApp(ui = ui, server = server)

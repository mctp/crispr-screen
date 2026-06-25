#!/usr/bin/env Rscript

# A shiny app to explore CRISPR screen results.
# Brian Magnuson (bmagnuso@umich.edu) -- Michigan Center for Translational Pathology

library(shiny)
library(tidyverse)
library(yaml)
library(ggrepel)
library(colourpicker)
library(shinyFiles)

# Use Cairo for bitmap rendering in headless server environments
# without this, ggplotly() fails with "unable to start device PNG"
if (!capabilities("X11")) {
    options(bitmapType = "cairo")
}

# Step 1: Conditional plotly import with availability flag
PLOTLY_AVAILABLE <- requireNamespace("plotly", quietly = TRUE)
if (PLOTLY_AVAILABLE) {
    library(plotly)
}

# Global UI constants
PICKER_WIDTH <- "60px"
PALETTE_WIDTH <- "60px"
DEFAULT_FDR_COLOR_1 <- "#FF0000"  # Red
DEFAULT_FDR_COLOR_2 <- "#FF7F00"  # Orange (actual colourpicker limited palette color)

# Note: We use colourpicker's built-in "limited" palette.
# DEFAULT_FDR_COLOR_2 is set to #FF7F00, which is the actual orange in that palette.
# (We previously tried #FFA500 but that color doesn't exist in the limited palette!)

# Inlined helper functions from ranked_gene_plots.R (previously sourced)

# define null-coalescing operator
`%||%` <- function(lhs, rhs) {
    if (!is.null(lhs) && length(lhs) > 0) lhs else rhs
}

is_docker <- function() {
    return(file.exists("/.dockerenv"))
}

parse_config <- function(file) {
    if (grepl("\\.yaml$", file) || grepl("\\.yml$", file)) {
        config_content <- readLines(file)
        config_content <- gsub(":\\s*\\.$", ": \"./\"", config_content)
        config <- yaml::yaml.load(paste(config_content, collapse = "\n"))
        names(config) <- tolower(names(config))
    } else {
        lines <- readLines(file, warn = FALSE)
        for (line in lines) {
            if (grepl("=", line) && !grepl("^#", line)) {
                parts <- strsplit(line, "=")[[1]]
                key <- trimws(parts[1])
                value <- trimws(parts[2])
                assign(key, value, envir = .GlobalEnv)
            }
        }
        config_values <- ls(.GlobalEnv)
        config_values <- config_values[!config_values %in% c("lines", "key_value", "key", "value")]
        config <- lapply(config_values, function(x) get(x, envir = .GlobalEnv))
        names(config) <- config_values
    }

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

    return(config)
}

# Extract comparison string from config entry (handles both old and new format)
extract_comparison_string <- function(comparison_entry) {
    if (is.character(comparison_entry)) {
        return(comparison_entry)
    } else if (is.list(comparison_entry) && !is.null(comparison_entry$comparison)) {
        return(comparison_entry$comparison)
    } else {
        stop("Invalid comparison format in config")
    }
}

# Extract plot config from comparison entry (returns NULL if not present)
extract_plot_config <- function(comparison_entry) {
    if (is.list(comparison_entry) && !is.null(comparison_entry$plot_config)) {
        return(comparison_entry$plot_config)
    }
    return(NULL)
}

parse_comparison_with_replicates <- function(comparison_str) {
    parts <- strsplit(comparison_str, ":")[[1]]
    if (length(parts) != 2) {
        stop(paste("Invalid comparison format:", comparison_str))
    }

    treatment_part <- parts[1]
    control_part <- parts[2]

    treatment_samples <- trimws(strsplit(treatment_part, ",")[[1]])
    control_samples <- trimws(strsplit(control_part, ",")[[1]])

    treatment_name <- generate_group_name(treatment_samples)
    control_name <- generate_group_name(control_samples)

    comparison_name <- paste0(treatment_name, "_vs_", control_name)

    return(list(
        treatment_samples = treatment_samples,
        control_samples = control_samples,
        comparison_name = comparison_name
    ))
}

generate_group_name <- function(sample_list) {
    if (length(sample_list) == 1) {
        return(sample_list[1])
    }
    base_names <- character(length(sample_list))
    for (i in seq_along(sample_list)) {
        sample <- sample_list[i]
        if (grepl("_\\d+$", sample)) {
            base_name <- sub("_\\d+$", "", sample)
            base_names[i] <- base_name
        } else {
            base_names[i] <- NA_character_
        }
    }
    if (all(!is.na(base_names)) && length(unique(base_names)) == 1) {
        return(paste0(base_names[1], "_reps"))
    }
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
                comparison_name = comparison_name,
                gene_file = gene_file,
                output_dir = rra_dir,
                original_comparison = comparison_str,
                plot_config = plot_config
            )
            cat("Found RRA results for:", comparison_str, "->", comparison_name, "\n")
            if (!is.null(plot_config)) {
                cat("  Plot config available with suffix:", plot_config$suffix %||% "(none)", "\n")
            }
        } else {
            cat("Warning: No RRA results found for:", comparison_str, "at", gene_file, "\n")
        }
    }
    return(result_files)
}
        
# Load gene list from file
load_gene_list <- function(file_path) {
    if (!file.exists(file_path)) {
        cat("Warning: Gene list file not found:", file_path, "\n")
        return(character(0))
    }
    genes <- readLines(file_path, warn = FALSE)
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    return(genes)
}

# takes a ggplot object and optionally converts it to plotly
render_plot_with_fallback <- function(ggplot_obj, use_plotly = TRUE, tooltip = "text") {
    if (PLOTLY_AVAILABLE && use_plotly) {
        return(ggplotly(ggplot_obj, tooltip = tooltip) %>% layout(hovermode = "closest"))
    } else {
        return(ggplot_obj)
    }
}

# UI Definition
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
                cursor: pointer;
                padding: 8px 10px;
                margin: 5px 0;
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 4px;
                user-select: none;
                display: block;
                text-decoration: none;
                color: #212529;
                font-weight: 500;
            }
            .section-header:hover {
                background-color: #e9ecef;
                text-decoration: none;
                color: #212529;
            }
            .section-header .arrow {
                float: right;
                transition: transform 0.2s;
            }
            .section-content {
                padding: 10px 5px;
                border-left: 2px solid #dee2e6;
                margin-left: 5px;
                margin-bottom: 10px;
            }
            /* Normalize selectize dropdown height and alignment for PALETTE dropdowns only */
            .palette-dropdown .selectize-control { margin-bottom: 0; }
            .palette-dropdown .selectize-control.single .selectize-input {
                min-height: 32px;
                height: 32px;
                padding-top: 4px;
                padding-bottom: 4px;
            }
            .palette-dropdown .selectize-control.single .selectize-input.input-active {
                min-height: 32px;
                height: 32px;
            }
            /* Ensure baseline alignment with neighboring inputs */
            .palette-dropdown .selectize-control.single { display: inline-flex; align-items: center; }
            /* Reduce internal input shim to avoid height jumps */
            .palette-dropdown .selectize-input > input { display: none; }
            /* Enforce fixed width on palette selectize to prevent content-based resizing */
            .palette-dropdown .selectize-control.single .selectize-input { width: ", PALETTE_WIDTH, " !important; max-width: ", PALETTE_WIDTH, " !important; min-width: ", PALETTE_WIDTH, " !important; }
        "))),
        tags$script(HTML("
            $(document).on('shiny:connected', function() {
                $('.section-header').on('click', function() {
                    var targetId = $(this).data('target');
                    var content = $('#' + targetId);
                    var arrow = $(this).find('.arrow');
                    
                    content.slideToggle(200);
                    if (arrow.text() === '▼') {
                        arrow.text('▶');
                    } else {
                        arrow.text('▼');
                    }
                });
                // Press Enter in config path input triggers Load Config
                $('#config_path_input').on('keypress', function(e) {
                    if (e.which === 13) {
                        e.preventDefault();
                        $('#load_config_path').click();
                    }
                });
            });
            // Copy config path to clipboard
            $(document).on('click', '#copy_config_path', function() {
                var val = $('#config_path_text').val();
                if (navigator.clipboard && window.isSecureContext) {
                    navigator.clipboard.writeText(val);
                } else {
                    var $temp = $('<textarea>');
                    $('body').append($temp);
                    $temp.val(val).select();
                    document.execCommand('copy');
                    $temp.remove();
                }
            });
            
            // Fix fdr_color_2 picker display when it shows black instead of orange
            // This is a workaround for colourpicker bug with orange in limited palette
            $(document).ready(function() {
                // Wait for pickers to initialize, then fix fdr_color_2 if it's black
                setTimeout(function() {
                    var picker2 = $('#fdr_color_2');
                    var expectedColor = '", DEFAULT_FDR_COLOR_2, "';
                    if (picker2.length && picker2.val() === '#000000') {
                        picker2.val(expectedColor);
                        picker2.css('background-color', expectedColor);
                        picker2.css('color', '#000');
                    }
                }, 500);
            });
        "))
    ),
    titlePanel("CRISPR Screen Explorer"),
    
    sidebarLayout(
        sidebarPanel(
            width = 3,
            class = "sidebar",
            
            # Config file selector (local browse, no upload)
            fluidRow(
                column(6, shinyFiles::shinyFilesButton("config_select", "Browse...", 
                                         title = "Select a config.yaml", multiple = FALSE)),
                column(6, actionButton("reload_default", "Reload Default", class = "btn-default btn-xs"))
            ),
            textInput("config_path_input", "Config Path:", value = "", placeholder = "Paste or browse for config.yaml"),
            actionButton("load_config_path", "Load Config", class = "btn-primary btn-sm", style = "width: 100%; margin-bottom: 8px;"),
            checkboxInput("show_debug", "Show debug messages", value = TRUE),
            uiOutput("config_error_message"),
            
            # Comparison selector
            selectInput("comparison", "Select Comparison:",
                       choices = NULL),
            uiOutput("global_sig_warning"),
            
            # Selection type
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
            
            # checkbox for interactive plots
            if (PLOTLY_AVAILABLE) {
                checkboxInput("use_interactive_plots", "Interactive Plots", value = TRUE)
            } else {
                tags$p(style = "font-size: 0.9em; color: #999; margin: 8px 0;",
                       "Interactive plots not available.")
            },
            
            # Collapsible sections
            tags$a(class = "section-header", `data-target` = "highlight_content",
                   "Highlight Gene Sets", tags$span(class = "arrow", "▼")),
            tags$div(id = "highlight_content", class = "section-content", style = "display: block;",
                uiOutput("gene_set_checkboxes")
            ),
            
            tags$a(class = "section-header", `data-target` = "significant_content",
                   "Significant Genes", tags$span(class = "arrow", "▶")),
            tags$div(id = "significant_content", class = "section-content", style = "display: none;",
                checkboxInput("highlight_significant", "Highlight significant genes", value = FALSE),
                numericInput("sig_top_n", "Max genes per threshold:", value = 10, min = 1, max = 1000, step = 1),
                checkboxInput("use_plot_thresholds", "Plot opt thresh", value = TRUE),
                uiOutput("sig_warning"),
                
                conditionalPanel(
                    condition = "!input.use_plot_thresholds",
                    fluidRow(
                        column(6,
                            checkboxInput("sig_use_thresh1", "Thresh 1", value = TRUE),
                            div(class = "palette-dropdown",
                                selectInput("sig_palette_1", "",
                                           choices = c("L" = "limited", "W" = "websafe", "F" = "square"),
                                           selected = "limited", width = PALETTE_WIDTH)
                            ),
                            uiOutput("sig_color_picker_1")
                        ),
                        column(6,
                            checkboxInput("sig_use_thresh2", "Thresh 2", value = TRUE),
                            div(class = "palette-dropdown",
                                selectInput("sig_palette_2", "",
                                           choices = c("L" = "limited", "W" = "websafe", "F" = "square"),
                                           selected = "limited", width = PALETTE_WIDTH)
                            ),
                            uiOutput("sig_color_picker_2")
                        )
                    )
                )
            ),
            
            tags$a(class = "section-header", `data-target` = "custom_content",
                   "Custom Gene List", tags$span(class = "arrow", "▶")),
            tags$div(id = "custom_content", class = "section-content", style = "display: none;",
                fileInput("custom_gene_file", "Upload Gene List File",
                         accept = c(".txt", ".csv", ".tsv")),
                fluidRow(
                    column(6, textInput("custom_list_name", "List Name:", value = "Custom")),
                    column(6, uiOutput("custom_color_picker"))
                ),
                div(class = "palette-dropdown",
                    selectInput("custom_palette_type", "",
                               choices = c("L" = "limited", "W" = "websafe", "F" = "square"),
                               selected = "limited", width = PALETTE_WIDTH)
                ),
                conditionalPanel(
                    condition = "input.custom_palette_type == 'websafe'",
                    selectInput("websafe_order", "Web-Safe Ordering:",
                               choices = c("RGB grid (R fast)" = "rgbgrid", "HSV (by hue)" = "hsv"),
                               selected = "rgbgrid")
                ),
                textAreaInput("custom_genes", "Gene List:",
                             placeholder = "GENE1 GENE2\nGENE3, GENE4",
                             rows = 3),
                fluidRow(
                    column(4, actionButton("apply_custom", "Apply", class = "btn-primary btn-sm", style = "width: 100%;")),
                    column(4, actionButton("clear_custom", "Clear Custom", class = "btn-default btn-sm", style = "width: 100%;")),
                    column(4, checkboxInput("show_custom", "Show on plot", value = TRUE))
                ),
                uiOutput("custom_warning")
            ),
            
                 tags$a(class = "section-header", `data-target` = "plot_content",
                     "Plot Options", tags$span(class = "arrow", "▶")),
                 tags$div(id = "plot_content", class = "section-content", style = "display: none;",
                    tags$style(HTML(
                        paste0(
                            ".compact-row .shiny-input-container{margin-bottom:0;}\n",
                            ".compact-row input:not([type=\"checkbox\"]), .compact-row select {height: 32px;}\n",
                            ".compact-row .checkbox {display: inline-block; width: auto; margin: 0;}\n",
                            ".compact-row .checkbox label {display: inline-flex; align-items: center; gap: 6px; margin: 0; white-space: nowrap;}\n",
                            ".compact-row .checkbox input[type=\"checkbox\"] {position: static; margin: 0;}\n"
                        )
                    )),
                textInput("volcano_title", "Volcano Title:", value = "", placeholder = "Auto"),
                textInput("ranked_title", "Ranked Title:", value = "", placeholder = "Auto"),
                sliderInput("point_size", "Base Point Size:", 
                           min = 0.5, max = 4, value = 2, step = 0.25),
                fluidRow(
                    column(12,
                        tags$div(tags$strong("FDR Threshold 1:")),
                            tags$div(
                                class = "compact-row",
                                style = "display: flex; gap: 8px; align-items: center; flex-wrap: wrap;",
                                textInput("fdr_threshold_1", label = NULL, value = "0.05", placeholder = "0.05", width = "60px"),
                            uiOutput("fdr_color_picker_1"),
                                div(class = "palette-dropdown",
                                    selectInput("fdr_palette_1", label = NULL,
                                            choices = c("L" = "limited", "W" = "websafe", "F" = "square"),
                                        selected = "limited", width = PALETTE_WIDTH)
                                ),
                            checkboxInput("show_fdr_line_1", label = "Show line", value = TRUE)
                        )
                    )
                ),
                fluidRow(
                    column(12,
                        tags$div(tags$strong("FDR Threshold 2:")),
                            tags$div(
                                class = "compact-row",
                                style = "display: flex; gap: 8px; align-items: center; flex-wrap: wrap;",
                                textInput("fdr_threshold_2", label = NULL, value = "0.1", placeholder = "0.1", width = "60px"),
                            uiOutput("fdr_color_picker_2"),
                                div(class = "palette-dropdown",
                                    selectInput("fdr_palette_2", label = NULL,
                                            choices = c("L" = "limited", "W" = "websafe", "F" = "square"),
                                        selected = "limited", width = PALETTE_WIDTH)
                                ),
                            checkboxInput("show_fdr_line_2", label = "Show line", value = TRUE)
                        )
                    )
                ),
                checkboxInput("show_labels", "Show Gene Labels", value = TRUE),
                sliderInput("label_scale", "Label Scale:", 
                           min = 0.5, max = 2.5, value = 1, step = 0.1),
                checkboxInput("scale_point_size", "Scale point size by FDR", value = TRUE),
                checkboxInput("color_by_fdr", "Color points by FDR", value = TRUE),
                sliderInput("text_scale", "Text Scale:", 
                           min = 0.5, max = 2.5, value = 1, step = 0.1)
            ),
            
            tags$a(class = "section-header", `data-target` = "download_content",
                   "Download Options", tags$span(class = "arrow", "▶")),
            tags$div(id = "download_content", class = "section-content", style = "display: none;",
                fluidRow(
                    column(6, selectInput("download_format", "Format:", 
                               choices = c("PNG" = "png", "PDF" = "pdf"),
                               selected = "png")),
                    column(6, selectInput("download_dpi", "DPI:", 
                               choices = c("72" = 72, "150" = 150, "300" = 300, 
                                          "600" = 600, "1200" = 1200),
                               selected = "300"))
                ),
                fluidRow(
                    column(6, sliderInput("download_width", "Width (in):", 
                               min = 4, max = 20, value = 10, step = 0.5)),
                    column(6, sliderInput("download_height", "Height (in):", 
                               min = 4, max = 20, value = 8, step = 0.5))
                ),
                downloadButton("download_volcano", "Download Volcano", style = "width: 100%; margin-bottom: 5px;"),
                downloadButton("download_ranked", "Download Ranked", style = "width: 100%;")
            )
        ),
        
        mainPanel(
            width = 9,
            fluidRow(
                column(12,
                    h3("Ranked Gene Plot"),
                    if (PLOTLY_AVAILABLE) {
                        conditionalPanel(
                            condition = "input.use_interactive_plots",
                            plotly::plotlyOutput("ranked_plot", height = "500px")
                        )
                    } else {
                        plotOutput("ranked_plot", height = "500px")
                    },
                    if (PLOTLY_AVAILABLE) {
                        conditionalPanel(
                            condition = "!input.use_interactive_plots",
                            plotOutput("ranked_plot_static", height = "500px")
                        )
                    }
                )
            ),
            hr(),
            fluidRow(
                column(12,
                    h3("Volcano Plot"),
                    if (PLOTLY_AVAILABLE) {
                        conditionalPanel(
                            condition = "input.use_interactive_plots",
                            plotly::plotlyOutput("volcano_plot", height = "500px")
                        )
                    } else {
                        plotOutput("volcano_plot", height = "500px")
                    },
                    if (PLOTLY_AVAILABLE) {
                        conditionalPanel(
                            condition = "!input.use_interactive_plots",
                            plotOutput("volcano_plot_static", height = "500px")
                        )
                    }
                )
            )
        )
    )
)

# Server Logic
server <- function(input, output, session) {
    cat("\n========== SERVER FUNCTION STARTING ==========\n")
    cat("CHECKPOINT 1: Function body entered\n")
    
    # Generate 216 web-safe colors (6x6x6 RGB cube)
    # Generate web-safe colours (6x6x6 RGB cube)
    cat("CHECKPOINT 2: About to create websafe colors\n")
    steps <- c("00", "33", "66", "99", "CC", "FF")
    websafe_colors_all <- as.vector(outer(steps, steps, FUN = function(r, g) paste0(r, g)))
    websafe_colors_all <- unlist(lapply(steps, function(b) paste0("#", paste0(websafe_colors_all, b))))
    
    # Two ordering strategies:
    # 1) RGB grid (R fast within G blocks, B groups) — matches many legacy tables
    websafe_colors_rgbgrid <- unlist(lapply(steps, function(b)
        unlist(lapply(steps, function(g)
            unlist(lapply(steps, function(r) paste0("#", r, g, b)))))))
    
    # 2) HSV ordering (Hue -> Saturation -> Value) — intuitive by hue families
    hsv_vals <- grDevices::rgb2hsv(col2rgb(websafe_colors_rgbgrid), maxColorValue = 255)
    ord_hsv <- order(hsv_vals["h", ], hsv_vals["s", ], hsv_vals["v", ])
    websafe_colors_hsv <- websafe_colors_rgbgrid[ord_hsv]
    
    cat("CHECKPOINT 3: Created websafe colors\n")
    
    # Default selection used by pickers; will be switched based on input
    websafe_colors <- websafe_colors_rgbgrid
    
    # Simple debug logging function - console only
    log_debug <- function(msg) {
        cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n"))
    }
    
    cat("CHECKPOINT 4: About to create reactive values\n")
    # Reactive values
    rv <- reactiveValues(
        config = NULL,
            config_dir = NULL,
        gene_data = NULL,
        comparisons = NULL,
        result_files = NULL,
        custom_gene_list = NULL,
        gene_sets_enabled = list(),
        warning_message = NULL,
        sig_warning = NULL,
        config_error = NULL,
        has_highlight_sets = FALSE,
        sig_color_1 = DEFAULT_FDR_COLOR_1,
        sig_color_2 = DEFAULT_FDR_COLOR_2,
        fdr_color_1 = DEFAULT_FDR_COLOR_1,
        fdr_color_2 = DEFAULT_FDR_COLOR_2
    )
    
    cat("CHECKPOINT 5: Reactive values created\n")
    
    # download format selector (HTML available only with plotly)
    output$download_format_selector <- renderUI({
        if (PLOTLY_AVAILABLE) {
            selectInput("download_format", "Format:", 
                       choices = c("PNG" = "png", "PDF" = "pdf", "HTML (Interactive)" = "html"),
                       selected = "png")
        } else {
            selectInput("download_format", "Format:", 
                       choices = c("PNG" = "png", "PDF" = "pdf"),
                       selected = "png")
        }
    })
    
    # Helper to clear all state and UI
    clear_state <- function() {
        rv$config <- NULL
        rv$gene_data <- NULL
        rv$result_files <- NULL
        rv$comparisons <- NULL
        updateSelectInput(session, "comparison", choices = NULL)
    }

    # Helper to normalize comparisons (handle both formats)
    normalize_comparisons <- function(comparisons) {
        if (is.null(comparisons) || length(comparisons) == 0) return(NULL)
        normalized <- list()
        for (item in comparisons) {
            if (is.character(item)) {
                # Format 1: just a value like "treatment:control"
                normalized <- c(normalized, list(item))
            } else if (is.list(item) && "comparison" %in% names(item)) {
                # Format 2: key:value like comparison: "treatment:control" with optional plot_config
                # Preserve the full structure to retain plot_config
                normalized <- c(normalized, list(item))
            } else if (is.list(item) && length(item) == 1 && is.character(item[[1]])) {
                # Other single-value list formats
                normalized <- c(normalized, list(item[[1]]))
            }
        }
        if (length(normalized) == 0) return(NULL)
        return(normalized)
    }

    # Helper to load config by path (sets base dir, results, comparison choices)
    load_config_from_path <- function(config_path) {
        cat("LOAD_CONFIG: Starting with path:", config_path, "\n")
        
        # Clear previous error
        rv$config_error <- NULL
        
        # Validate path
        if (is.null(config_path) || !nzchar(config_path)) {
            cat("LOAD_CONFIG: ERROR - path is null or empty\n")
            rv$config_error <- "No config file selected."
            clear_state()
            return(FALSE)
        }
        if (!file.exists(config_path)) {
            cat("LOAD_CONFIG: ERROR - file does not exist:", config_path, "\n")
            rv$config_error <- paste0("Config file not found: ", config_path)
            clear_state()
            return(FALSE)
        }
        
        cat("LOAD_CONFIG: File exists, parsing...\n")
        # Parse config with error handling
        rv$config_dir <- dirname(config_path)
        tryCatch({
            rv$config <- parse_config(config_path)
            cat("LOAD_CONFIG: Config parsed successfully\n")
        }, error = function(e) {
            cat("LOAD_CONFIG: ERROR parsing config:", e$message, "\n")
            rv$config_error <- paste0("Error parsing config: ", e$message)
            rv$config <- NULL
        })
        
        if (is.null(rv$config)) {
            if (is.null(rv$config_error)) {
                rv$config_error <- "Failed to parse config file (returned NULL)."
            }
            cat("LOAD_CONFIG: ERROR - config is NULL\n")
            clear_state()
            return(FALSE)
        }
        
        cat("LOAD_CONFIG: Config loaded, checking comparisons field\n")
        # Validate comparisons field
        comparisons_raw <- rv$config$comparisons
        if (is.null(comparisons_raw)) {
            cat("LOAD_CONFIG: ERROR - no comparisons field in config\n")
            rv$config_error <- "Config missing 'comparisons' field."
            clear_state()
            return(FALSE)
        }
        
        cat("LOAD_CONFIG: Normalizing comparisons\n")
        # Normalize comparisons
        comparisons <- normalize_comparisons(comparisons_raw)
        if (is.null(comparisons) || length(comparisons) == 0) {
            cat("LOAD_CONFIG: ERROR - comparisons is null or empty after normalization\n")
            rv$config_error <- "Config 'comparisons' field is empty or invalid."
            clear_state()
            return(FALSE)
        }
        
        cat("LOAD_CONFIG: Found", length(comparisons), "comparisons\n")
        # Load output dir and experiment name
        output_dir <- rv$config$output_dir %||% "output"
        experiment_name <- rv$config$experiment_name %||% "experiment"
        cat("LOAD_CONFIG: output_dir =", output_dir, ", experiment_name =", experiment_name, "\n")
        
        # Resolve output_dir: try relative to config first, then absolute
        output_dir_relative <- file.path(rv$config_dir, output_dir)
        cat("LOAD_CONFIG: Checking output_dir_relative:", output_dir_relative, "\n")
        if (file.exists(output_dir_relative)) {
            output_dir <- output_dir_relative
            cat("LOAD_CONFIG: Using relative path:", output_dir, "\n")
        } else if (!file.exists(output_dir)) {
            cat("LOAD_CONFIG: ERROR - neither relative nor absolute path exists\n")
            # Neither relative nor absolute path exists
            rv$config_error <- paste0("Output directory not found. Tried: ", output_dir_relative, " and ", output_dir)
            clear_state()
            return(FALSE)
        } else {
            cat("LOAD_CONFIG: Using absolute path:", output_dir, "\n")
        }
        
        cat("LOAD_CONFIG: Finding RRA results in:", output_dir, "\n")
        # Find result files
        tryCatch({
            rv$result_files <- find_rra_results(output_dir, experiment_name, comparisons)
            cat("LOAD_CONFIG: Found", length(rv$result_files), "result files\n")
            if (is.null(rv$result_files) || length(rv$result_files) == 0) {
                cat("LOAD_CONFIG: ERROR - no result files found\n")
                rv$config_error <- paste0("No result files found in: ", output_dir)
                clear_state()
                return(FALSE)
            }
            comparison_names <- names(rv$result_files)
            cat("LOAD_CONFIG: Comparison names:", paste(comparison_names, collapse = ", "), "\n")
            updateSelectInput(session, "comparison", choices = comparison_names)
            cat("LOAD_CONFIG: Updated comparison dropdown\n")
        }, error = function(e) {
            cat("LOAD_CONFIG: ERROR finding results:", e$message, "\n")
            rv$config_error <- paste0("Error finding result files: ", e$message)
            clear_state()
            return(FALSE)
        })
        
        # Update path input field
        updateTextInput(session, "config_path_input", value = config_path)
        cat("LOAD_CONFIG: SUCCESS\n")
        return(TRUE)
    }

    # Configure local file browser for config.yaml selection
    vols <- shinyFiles::getVolumes()()
    roots <- vols
    default_root_path <- Sys.getenv("APP_ROOT_PATH", unset = "/data")
    if (file.exists(default_root_path)) {
        roots <- c(roots, CRISPR = default_root_path)
    }
    roots <- c(roots, WD = getwd(), Home = normalizePath("~", winslash = "/", mustWork = FALSE))

    # Note: shinyFileChoose doesn't support allowDirCreate parameter (only shinyDirChoose does)
    # Sorting is available via the "Sort content" dropdown button in the file browser
    if ("CRISPR" %in% names(roots)) {
        shinyFiles::shinyFileChoose(input, id = "config_select", roots = roots,
                                    filetypes = c("yaml", "yml"),
                                    defaultRoot = "CRISPR", defaultPath = ".")
    } else if ("WD" %in% names(roots)) {
        shinyFiles::shinyFileChoose(input, id = "config_select", roots = roots,
                                    filetypes = c("yaml", "yml"),
                                    defaultRoot = "WD", defaultPath = ".")
    } else {
        shinyFiles::shinyFileChoose(input, id = "config_select", roots = roots,
                                    filetypes = c("yaml", "yml"))
    }

    # Render config error messages
    output$config_error_message <- renderUI({
        if (!is.null(rv$config_error)) {
            tags$div(style = "margin-top: 8px; padding: 8px; background-color: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; border-radius: 4px; font-size: 0.9em;",
                tags$strong("Config Error: "), rv$config_error)
        }
    })

    # Debug output removed - using console logging only

    observeEvent(input$config_select, {
        paths <- shinyFiles::parseFilePaths(roots, input$config_select)
        if (nrow(paths) > 0) {
            cfg <- as.character(paths$datapath[1])
            updateTextInput(session, "config_path_input", value = cfg)
            load_config_from_path(cfg)
        }
    }, ignoreInit = TRUE)
    
    # Debug checkbox observer - test that debugging works
    observeEvent(input$show_debug, {
        if (input$show_debug) {
            cat("[DEBUG] DEBUG PANEL OPENED - Debugging is working!\n")
        }
    })
    
    # Load config from text input button
    observeEvent(input$load_config_path, {
        if (!is.null(input$config_path_input) && nzchar(input$config_path_input)) {
            load_config_from_path(input$config_path_input)
        } else {
            rv$config_error <- "Please enter a config path."
        }
    })
    
    # Reload default config button
    observeEvent(input$reload_default, {
        default_config <- file.path(default_root_path, "config.yaml")
        if (file.exists(default_config)) {
            load_config_from_path(default_config)
        } else if (file.exists("config.yaml")) {
            load_config_from_path("config.yaml")
        } else {
            rv$config_error <- "No default config found."
            clear_state()
        }
    })

    # On startup, try default paths without uploading
    cat("CHECKPOINT 6: About to enter isolate block\n")
    isolate({
        cat("ISOLATE BLOCK: Starting isolate block\n")
        default_config <- file.path(default_root_path, "config.yaml")
        cat("ISOLATE BLOCK: default_root_path =", default_root_path, "\n")
        cat("ISOLATE BLOCK: default_config =", default_config, "\n")
        cat("ISOLATE BLOCK: file.exists(default_config) =", file.exists(default_config), "\n")
        
        if (file.exists(default_config)) {
            cat("ISOLATE BLOCK: Loading config from:", default_config, "\n")
            load_config_from_path(default_config)
            cat("ISOLATE BLOCK: Config loaded successfully\n")
        } else if (file.exists("config.yaml")) {
            cat("ISOLATE BLOCK: Loading config from current dir: config.yaml\n")
            load_config_from_path("config.yaml")
            cat("ISOLATE BLOCK: Config loaded successfully\n")
        } else {
            cat("ISOLATE BLOCK: No config.yaml found\n")
        }
        cat("ISOLATE BLOCK: Finished\n")
    })
    
    # Load gene data when comparison changes
    observeEvent(input$comparison, {
        req(input$comparison, rv$result_files)
        
        result_info <- rv$result_files[[input$comparison]]
        if (!is.null(result_info) && file.exists(result_info$gene_file)) {
            rv$gene_data <- read.delim(result_info$gene_file, sep = "\t", header = TRUE)
            
            # Check if highlight sets exist
            has_sets <- !is.null(result_info$plot_config) && 
                       !is.null(result_info$plot_config$highlight_genes) &&
                       length(result_info$plot_config$highlight_genes) > 0
            rv$has_highlight_sets <- has_sets
            
            # Auto-enable significant genes if no predefined sets AND no custom gene list
            if (!has_sets && (is.null(rv$custom_gene_list) || length(rv$custom_gene_list$genes) == 0)) {
                updateCheckboxInput(session, "highlight_significant", value = TRUE)
            }
            
            # Initialize gene set checkboxes state
            if (has_sets) {
                for (i in seq_along(result_info$plot_config$highlight_genes)) {
                    group_info <- result_info$plot_config$highlight_genes[[i]]
                    group_name <- group_info$group_name %||% paste0("Group ", i)
                    if (is.null(rv$gene_sets_enabled[[group_name]])) {
                        rv$gene_sets_enabled[[group_name]] <- TRUE
                    }
                }
            }
        }
    })

    # Clear custom gene list handler
    observeEvent(input$clear_custom, {
        rv$custom_gene_list <- NULL
        updateTextAreaInput(session, "custom_genes", value = "")
        rv$warning_message <- NULL
        # If no predefined highlight sets, default to significant genes
        if (!is.null(input$comparison) && !is.null(rv$result_files)) {
            result_info <- rv$result_files[[input$comparison]]
            has_sets <- !is.null(result_info$plot_config) &&
                       !is.null(result_info$plot_config$highlight_genes) &&
                       length(result_info$plot_config$highlight_genes) > 0
            if (!has_sets) {
                updateCheckboxInput(session, "highlight_significant", value = TRUE)
            }
        }
    })
    
    # Render gene set checkboxes with editable names and colors
    output$gene_set_checkboxes <- renderUI({
        req(input$comparison, rv$result_files)
        
        result_info <- rv$result_files[[input$comparison]]
        if (is.null(result_info$plot_config) || 
            is.null(result_info$plot_config$highlight_genes)) {
            return(p("No predefined gene sets available"))
        }
        
        gene_sets <- result_info$plot_config$highlight_genes
        controls <- lapply(seq_along(gene_sets), function(i) {
            group_info <- gene_sets[[i]]
            group_name <- group_info$group_name %||% paste0("Group ", i)
            group_color <- group_info$color %||% "black"
            
            div(
                style = "border: 1px solid #ddd; padding: 5px; margin-bottom: 5px; border-radius: 3px;",
                fluidRow(
                    column(3, checkboxInput(
                        paste0("geneset_", i),
                        label = NULL,
                        value = rv$gene_sets_enabled[[group_name]] %||% TRUE
                    )),
                    column(9, textInput(
                        paste0("geneset_name_", i),
                        NULL,
                        value = group_name,
                        placeholder = "Name"
                    ))
                ),
                fluidRow(
                    column(6, div(class = "palette-dropdown",
                        selectInput(
                            paste0("geneset_palette_", i),
                            NULL,
                            choices = c("F" = "square", "W" = "websafe", "L" = "limited"),
                            selected = "square"
                        )
                    )),
                    column(6, uiOutput(paste0("geneset_color_picker_", i)))
                )
            )
        })
        
        do.call(tagList, controls)
    })
    
    # FDR threshold color pickers (Plot Options) - respond to palette selection
    output$fdr_color_picker_1 <- renderUI({
        palette_type <- if (is.null(input$fdr_palette_1)) "limited" else input$fdr_palette_1
        # Use isolate to prevent reactive dependency on input$fdr_color_1
        current_color <- isolate({
            if (!is.null(input$fdr_color_1)) input$fdr_color_1 else DEFAULT_FDR_COLOR_1
        })
        
        # isolate({
        #     rv$debug_info <- paste0(rv$debug_info, "\n[RENDER_1] Creating picker with palette=", palette_type, ", value=", current_color)
        # })
        
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput(
                "fdr_color_1", NULL,
                value = current_color,
                palette = "limited",
                allowedCols = colors_vec,
                returnName = FALSE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        } else {
            # Use built-in palette (limited or square)
            colourpicker::colourInput(
                "fdr_color_1", NULL,
                value = current_color,
                palette = palette_type,
                returnName = FALSE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        }
    })
    
    output$fdr_color_picker_2 <- renderUI({
        palette_type <- if (is.null(input$fdr_palette_2)) "limited" else input$fdr_palette_2
        # Use isolate to prevent reactive dependency on input$fdr_color_2
        current_color <- isolate({
            if (!is.null(input$fdr_color_2)) input$fdr_color_2 else DEFAULT_FDR_COLOR_2
        })
        
        # isolate({
        #     rv$debug_info <- paste0(rv$debug_info, "\n[RENDER_2] Creating picker with palette=", palette_type, ", value=", current_color)
        # })
        
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput(
                "fdr_color_2", NULL,
                value = current_color,
                palette = "limited",
                allowedCols = colors_vec,
                returnName = FALSE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        } else {
            # Use built-in palette (limited or square)
            colourpicker::colourInput(
                "fdr_color_2", NULL,
                value = current_color,
                palette = palette_type,
                returnName = FALSE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        }
    })

    # Ensure FDR color pickers render even when hidden (collapsed section)
    outputOptions(output, "fdr_color_picker_1", suspendWhenHidden = FALSE)
    outputOptions(output, "fdr_color_picker_2", suspendWhenHidden = FALSE)
    
    # Force correct colors after pickers initialize
    # Helper function to get the correct FDR color, working around colourpicker bug
    # If picker shows black (#000000), it means it failed to initialize with orange, so use default
    get_fdr_color_1 <- reactive({
        color <- input$fdr_color_1
        if (is.null(color)) return(DEFAULT_FDR_COLOR_1)
        # isolate({
        #     rv$debug_info <- paste0(rv$debug_info, "\n[1] get_fdr_color_1: input=", color)
        # })
        return(color)
    })
    
    get_fdr_color_2 <- reactive({
        color <- input$fdr_color_2
        if (is.null(color)) {
            # isolate({
            #     rv$debug_info <- paste0(rv$debug_info, "\n[2] get_fdr_color_2: NULL, using ", DEFAULT_FDR_COLOR_2)
            # })
            return(DEFAULT_FDR_COLOR_2)
        }
        # If picker initialized to black instead of orange, override with orange
        if (color == "#000000") {
            # isolate({
            #     rv$debug_info <- paste0(rv$debug_info, "\n[2] get_fdr_color_2: got BLACK, overriding with ", DEFAULT_FDR_COLOR_2)
            # })
            return(DEFAULT_FDR_COLOR_2)
        }
        # isolate({
        #     rv$debug_info <- paste0(rv$debug_info, "\n[2] get_fdr_color_2: using ", color)
        # })
        return(color)
    })
    
    # Monitor and try to fix the visual display of fdr_color_2 when it shows black
    observe({
        req(!is.null(input$fdr_color_2))
        # isolate({
        #     rv$debug_info <- paste0(rv$debug_info, "\n[OBSERVE] Checking fdr_color_2, value=", input$fdr_color_2)
        # })
        if (input$fdr_color_2 == "#000000") {
            # isolate({
            #     rv$debug_info <- paste0(rv$debug_info, " -> ATTEMPTING JS FIX")
            # })
            # Use JavaScript to manually update the input's visual appearance
            js_code <- sprintf("
                console.log('JS: Attempting to fix fdr_color_2');
                var picker = $('#fdr_color_2');
                console.log('JS: picker found:', picker.length > 0);
                if (picker.length) {
                    console.log('JS: Current value:', picker.val());
                    picker.val('%s');
                    picker.css('background-color', '%s');
                    console.log('JS: New value set to:', picker.val());
                    picker.trigger('change');
                }
            ", DEFAULT_FDR_COLOR_2, DEFAULT_FDR_COLOR_2)
            shinyjs::runjs(js_code)
        } else {
            # isolate({
            #     rv$debug_info <- paste0(rv$debug_info, " -> OK (not black)")
            # })
        }
    })
    
    # Monitor what renderUI is actually creating
    observe({
        req(!is.null(input$fdr_palette_2))
        # isolate({
        #     current_input_val <- if (!is.null(input$fdr_color_2)) input$fdr_color_2 else "NULL"
        #     rv$debug_info <- paste0(rv$debug_info, 
        #         "\n[RENDER_MONITOR] fdr_color_picker_2 may have rendered, palette=", 
        #         input$fdr_palette_2, ", current input$fdr_color_2=", current_input_val)
        # })
    })
    
    # Significant genes color pickers
    output$sig_color_picker_1 <- renderUI({
        palette_type <- if (is.null(input$sig_palette_1)) "limited" else input$sig_palette_1
        current_color <- if (!is.null(input$sig_color_1)) input$sig_color_1 else rv$sig_color_1
        
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput(
                "sig_color_1", NULL,
                value = current_color,
                palette = "limited",
                allowedCols = colors_vec,
                returnName = TRUE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        } else {
            colourpicker::colourInput(
                "sig_color_1", NULL,
                value = current_color,
                palette = palette_type,
                returnName = TRUE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        }
    })
    
    output$sig_color_picker_2 <- renderUI({
        palette_type <- if (is.null(input$sig_palette_2)) "limited" else input$sig_palette_2
        current_color <- if (!is.null(input$sig_color_2)) input$sig_color_2 else rv$sig_color_2
        
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput(
                "sig_color_2", NULL,
                value = current_color,
                palette = "limited",
                allowedCols = colors_vec,
                returnName = TRUE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        } else {
            colourpicker::colourInput(
                "sig_color_2", NULL,
                value = current_color,
                palette = palette_type,
                returnName = TRUE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        }
    })
    
    output$sig_color_picker_1_plot <- renderUI({
        palette_type <- if (is.null(input$sig_palette_1_plot)) "limited" else input$sig_palette_1_plot
        current_color <- if (!is.null(input$sig_color_1_plot)) input$sig_color_1_plot else rv$sig_color_1
        
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput(
                "sig_color_1_plot", NULL,
                value = current_color,
                palette = "limited",
                allowedCols = colors_vec,
                returnName = TRUE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        } else {
            colourpicker::colourInput(
                "sig_color_1_plot", NULL,
                value = current_color,
                palette = palette_type,
                returnName = TRUE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        }
    })
    
    output$sig_color_picker_2_plot <- renderUI({
        palette_type <- if (is.null(input$sig_palette_2_plot)) "limited" else input$sig_palette_2_plot
        current_color <- if (!is.null(input$sig_color_2_plot)) input$sig_color_2_plot else rv$sig_color_2
        
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourpicker::colourInput(
                "sig_color_2_plot", NULL,
                value = current_color,
                palette = "limited",
                allowedCols = colors_vec,
                returnName = TRUE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        } else {
            colourpicker::colourInput(
                "sig_color_2_plot", NULL,
                value = current_color,
                palette = palette_type,
                returnName = TRUE,
                showColour = "both",
                allowTransparent = FALSE,
                width = PICKER_WIDTH
            )
        }
    })
    
    # Update FDR threshold colors when changed in Plot Options
    # NOTE: These observeEvent blocks are REMOVED because they were causing initialization issues.
    # The problem: updateColourInput() triggers observeEvent even with ignoreInit=TRUE,
    # which then overwrites rv with the wrong initial value (black instead of orange).
    # Instead, rv$fdr_color_1 and rv$fdr_color_2 remain as DEFAULT constants,
    # and the renderUI blocks read directly from input$fdr_color_X for user changes.
    # 
    # observeEvent(input$fdr_color_1, {
    #     if (!is.null(input$fdr_color_1)) {
    #         rv$fdr_color_1 <- input$fdr_color_1
    #         if (!is.null(input$use_plot_thresholds) && input$use_plot_thresholds) {
    #             rv$sig_color_1 <- input$fdr_color_1
    #         }
    #     }
    # }, ignoreInit = TRUE)
    # 
    # observeEvent(input$fdr_color_2, {
    #     if (!is.null(input$fdr_color_2)) {
    #         rv$fdr_color_2 <- input$fdr_color_2
    #         if (!is.null(input$use_plot_thresholds) && input$use_plot_thresholds) {
    #             rv$sig_color_2 <- input$fdr_color_2
    #         }
    #     }
    # }, ignoreInit = TRUE)
    
    # Update stored colors when changed
    observeEvent(input$sig_color_1, {
        if (!is.null(input$sig_color_1)) rv$sig_color_1 <- input$sig_color_1
    })
    observeEvent(input$sig_color_2, {
        if (!is.null(input$sig_color_2)) rv$sig_color_2 <- input$sig_color_2
    })
    observeEvent(input$sig_color_1_plot, {
        if (!is.null(input$sig_color_1_plot)) rv$sig_color_1 <- input$sig_color_1_plot
    })
    observeEvent(input$sig_color_2_plot, {
        if (!is.null(input$sig_color_2_plot)) rv$sig_color_2 <- input$sig_color_2_plot
    })
    
    # When switching to "Plot opt thresh", sync colors from plot options
    observeEvent(input$use_plot_thresholds, {
        if (!is.null(input$use_plot_thresholds) && input$use_plot_thresholds) {
            # Set to current FDR threshold colors from Plot Options
            rv$sig_color_1 <- get_fdr_color_1()
            rv$sig_color_2 <- get_fdr_color_2()
        }
    })
    
    # Render custom color picker based on palette selection
    output$custom_color_picker <- renderUI({
        palette_type <- if (is.null(input$custom_palette_type)) "square" else input$custom_palette_type
        current_color <- if (is.null(input$custom_color)) "#00FF00" else input$custom_color
        
        if (palette_type == "websafe") {
            colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
            colourInput(
                "custom_color", NULL,
                value = current_color,
                palette = "limited",
                allowedCols = colors_vec,
                returnName = TRUE,
                showColour = "both",
                width = PICKER_WIDTH
            )
        } else if (palette_type == "limited") {
            colourInput(
                "custom_color", NULL,
                value = current_color,
                palette = "limited",
                returnName = TRUE,
                showColour = "both",
                width = PICKER_WIDTH
            )
        } else {
            colourInput(
                "custom_color", NULL,
                value = current_color,
                palette = "square",
                returnName = TRUE,
                showColour = "both",
                width = PICKER_WIDTH
            )
        }
    })

    # Ensure custom color picker renders even when section is collapsed
    outputOptions(output, "custom_color_picker", suspendWhenHidden = FALSE)
    
    # Dynamically render color pickers for each gene set
    observe({
        req(input$comparison, rv$result_files)
        
        result_info <- rv$result_files[[input$comparison]]
        if (!is.null(result_info$plot_config) && 
            !is.null(result_info$plot_config$highlight_genes)) {
            
            gene_sets <- result_info$plot_config$highlight_genes
            lapply(seq_along(gene_sets), function(i) {
                output[[paste0("geneset_color_picker_", i)]] <- renderUI({
                    palette_type <- input[[paste0("geneset_palette_", i)]]
                    if (is.null(palette_type)) palette_type <- "square"
                    
                    current_color <- input[[paste0("geneset_color_", i)]]
                    if (is.null(current_color)) {
                        current_color <- gene_sets[[i]]$color %||% "black"
                    }
                    
                    if (palette_type == "websafe") {
                        colors_vec <- if (!is.null(input$websafe_order) && input$websafe_order == "hsv") websafe_colors_hsv else websafe_colors
                        colourInput(
                            paste0("geneset_color_", i),
                            NULL,
                            value = current_color,
                            palette = "limited",
                            allowedCols = colors_vec,
                            returnName = TRUE,
                            showColour = "both",
                            width = PICKER_WIDTH
                        )
                    } else if (palette_type == "limited") {
                        colourInput(
                            paste0("geneset_color_", i),
                            NULL,
                            value = current_color,
                            palette = "limited",
                            returnName = TRUE,
                            showColour = "both",
                            width = PICKER_WIDTH
                        )
                    } else {
                        colourInput(
                            paste0("geneset_color_", i),
                            NULL,
                            value = current_color,
                            palette = "square",
                            returnName = TRUE,
                            showColour = "both",
                            width = PICKER_WIDTH
                        )
                    }
                })
            })
        }
    })
    
    # Load custom gene list from file
    observeEvent(input$custom_gene_file, {
        req(input$custom_gene_file)
        
        # Read file handling different line endings (Unix, Mac, Windows)
        file_content <- readLines(input$custom_gene_file$datapath, warn = FALSE)

        # Preserve one-gene-per-line from uploaded file. Trim and drop empty lines.
        genes <- trimws(file_content)
        genes <- genes[genes != ""]

        if (length(genes) > 0) {
            # Update the text area with loaded genes using newlines (preserve original delimiter)
            updateTextAreaInput(session, "custom_genes", value = paste(genes, collapse = "\n"))
        }
    })
    
    # Apply custom gene list with validation
    observeEvent(input$apply_custom, {
        req(input$custom_genes)
        cat("[APPLY_CUSTOM] Apply clicked. custom_genes length:", nchar(as.character(input$custom_genes)), " rv$gene_data present:", !is.null(rv$gene_data), "\n")
            raw <- input$custom_genes
            # Normalize newlines
            raw <- gsub("\r\n", "\n", raw)
            raw <- gsub("\r", "\n", raw)
            lines <- unlist(strsplit(raw, "\n"))
            lines <- lines[trimws(lines) != ""]

            # Tokenize entire text for fallback stats
            # Try delimiters in priority order: newline, comma, semicolon, colon, whitespace.
            raw_norm <- raw
            raw_norm <- gsub("\r\n", "\n", raw_norm)
            raw_norm <- gsub("\r", "\n", raw_norm)

            try_split <- function(text, delim_pattern, split_func) {
                toks <- split_func(text)
                toks <- trimws(toks)
                toks <- toks[toks != ""]
                # if any token still contains forbidden separators, it's likely wrong split
                bad <- any(grepl("[[[:space:]],;:]", toks))
                return(list(tokens = toks, ok = !bad))
            }

            # Define split functions for each candidate delimiter
            split_newline <- function(t) unlist(strsplit(t, "\n"))
            split_comma <- function(t) unlist(strsplit(t, ","))
            split_semicolon <- function(t) unlist(strsplit(t, ";"))
            split_colon <- function(t) unlist(strsplit(t, ":"))
            split_whitespace <- function(t) unlist(strsplit(t, "\\s+"))

            candidates_try <- list(
                newline = try_split(raw_norm, "\n", split_newline),
                comma = try_split(raw_norm, ",", split_comma),
                semicolon = try_split(raw_norm, ";", split_semicolon),
                colon = try_split(raw_norm, ":", split_colon),
                whitespace = try_split(raw_norm, "\\s+", split_whitespace)
            )

            # pick first 'ok' candidate in priority order; else use whitespace fallback
            delim_order <- c("newline", "comma", "semicolon", "colon", "whitespace")
            chosen <- NULL
            for (d in delim_order) {
                if (length(candidates_try[[d]]$tokens) > 0 && candidates_try[[d]]$ok) {
                    chosen <- d
                    tokens <- candidates_try[[d]]$tokens
                    break
                }
            }
            if (is.null(chosen)) {
                chosen <- "whitespace"
                tokens <- candidates_try$whitespace$tokens
            }

            # strip surrounding punctuation from tokens (commas, semicolons, colons, quotes)
            tokens <- gsub("^[\"'`,;:[[:space:]]]+|[\"'`,;:[[:space:]]]+$", "", tokens)
            tokens <- tokens[tokens != ""]

            all_tokens <- tokens

            # Build line_last fallback: treat each line as entry unless it contains explicit separators
            lines_norm <- unlist(strsplit(raw_norm, "\n"))
            lines_norm <- trimws(lines_norm)
            lines_norm <- lines_norm[lines_norm != ""]
            if (length(lines_norm) > 0) {
                line_last <- sapply(lines_norm, function(l) {
                    if (grepl("[,;:]", l)) {
                        parts <- unlist(strsplit(l, "[,;:]+"))
                        parts <- trimws(parts)
                        parts <- parts[parts != ""]
                        if (length(parts) == 0) return(NA_character_)
                        return(tail(parts, 1))
                    } else {
                        return(l)
                    }
                })
                line_last <- trimws(line_last)
                line_last <- line_last[!is.na(line_last) & line_last != ""]
            } else {
                line_last <- character(0)
            }

            cat("[APPLY_CUSTOM] Tokenization chosen delimiter:", chosen, "-> tokens:", length(all_tokens), "line_last:", length(line_last), "\n")

            if (length(all_tokens) == 0 && length(line_last) == 0) {
                rv$warning_message <- "No tokens parsed from input."
                rv$custom_gene_list <- NULL
                return()
            }

            # Candidate genes: prefer last-column tokens; include all tokens as backup
            candidates <- unique(c(line_last, all_tokens))
        
            # Filter out non-gene patterns
            candidates <- candidates[!grepl("^[0-9]+$", candidates)]               # pure numbers
            candidates <- candidates[!grepl("^[ACGTNacgtn]{15,}$", candidates)]     # long nucleotide strings
            candidates <- candidates[!grepl("^.*_[0-9]+$", candidates)]             # sgRNA IDs with underscore + digits
        
            if (length(candidates) == 0) {
                rv$warning_message <- paste0(
                    "Parsed ", length(all_tokens), " tokens (", length(unique(all_tokens)), " unique). After filtering, no gene-like tokens left.")
                rv$custom_gene_list <- NULL
                return()
            }
            genes <- candidates

        # Normalize tokens: strip hidden/non-printable characters and trim
        genes <- gsub("[^[:print:]]", "", genes)
        genes <- trimws(genes)
        genes <- genes[genes != ""]

        # If no gene dataset is loaded yet, accept the list but skip validation
        if (is.null(rv$gene_data)) {
            rv$custom_gene_list <- list(
                name = input$custom_list_name %||% "Custom",
                genes = unique(genes),
                color = input$custom_color
            )
            rv$warning_message <- paste0("Loaded ", length(unique(genes)), " candidate genes; dataset not loaded so no validation performed.")
            cat("[APPLY_CUSTOM] No gene data available; stored ", length(unique(genes)), " genes for later validation.\n")
            return()
        }
        
        # Get gene identifier column from current data (auto-select by best overlap)
        colnames_clean <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(rv$gene_data)))
        gene_data_temp <- rv$gene_data
        colnames(gene_data_temp) <- colnames_clean

        candidate_cols <- unique(c(
            intersect(c("id", "gene"), colnames_clean),
            grep("gene|symbol|name", colnames_clean, value = TRUE)
        ))
        if (length(candidate_cols) == 0) {
            candidate_cols <- colnames_clean
        }

        best_col <- NULL
        best_hits <- -1
        best_avail <- NULL

        for (col in candidate_cols) {
            avail <- gene_data_temp[[col]]
            hits <- sum(genes %in% avail, na.rm = TRUE)
            if (hits > best_hits) {
                best_hits <- hits
                best_col <- col
                best_avail <- avail
            }
        }

        gene_id_col <- if (!is.null(best_col)) best_col else colnames_clean[1]
        available_genes <- if (!is.null(best_avail)) best_avail else gene_data_temp[[gene_id_col]]
        cat("[APPLY_CUSTOM] gene_id_col chosen:", gene_id_col, "(hits:", best_hits, ") first values:", paste(head(available_genes, 5), collapse=", "), "\n")
        
        # Check which genes match (case-sensitive)
        matched_mask <- genes %in% available_genes
        # Map matched tokens back to the dataset's exact casing so later %in% checks succeed
        matched_tokens <- genes[matched_mask]
        if (length(matched_tokens) > 0) {
            matched_genes <- unique(available_genes[match(matched_tokens, available_genes)])
        } else {
            matched_genes <- character(0)
        }
        unmatched_genes <- genes[!matched_mask]
        
        if (length(matched_genes) == 0) {
            rv$warning_message <- paste0(
                "Parsed ", length(all_tokens), " tokens (", length(unique(all_tokens)), " unique; ", length(genes), " candidate gene tokens). ",
                "None matched genes in dataset. Candidate tokens: ",
                paste(head(unmatched_genes, 25), collapse = " "),
                if (length(unmatched_genes) > 25) paste0(" (and ", length(unmatched_genes) - 25, " more)") else ""
            )
            rv$custom_gene_list <- NULL
        } else {
            rv$custom_gene_list <- list(
                name = input$custom_list_name,
                genes = unique(matched_genes),  # deduplicate matched genes
                color = input$custom_color
            )
            
            if (length(unmatched_genes) > 0) {
                rv$warning_message <- paste0(
                    "Parsed ", length(all_tokens), " tokens (", length(unique(all_tokens)), " unique; ", length(genes), " candidate). ",
                    "Matched ", length(matched_genes), " genes. ",
                    length(unmatched_genes), " non-matching candidate tokens: ",
                    paste(head(unmatched_genes, 20), collapse = " "),
                    if (length(unmatched_genes) > 20) paste0(" (and ", length(unmatched_genes) - 20, " more)") else ""
                )
            } else {
                rv$warning_message <- paste0(
                    "Successfully loaded ", length(matched_genes), " unique genes from ",
                    length(all_tokens), " tokens (", length(unique(all_tokens)), " unique; ", length(genes), " candidate)."
                )
            }
        }
    })
    
    # Render custom list warnings
    output$custom_warning <- renderUI({
        if (!is.null(rv$warning_message)) {
            warning_color <- if (grepl("Successfully", rv$warning_message)) {
                "green"
            } else if (grepl("None matched|No tokens parsed", rv$warning_message)) {
                "red"
            } else {
                "orange"
            }
            
            div(
                style = paste0("margin-top: 10px; padding: 8px; border-radius: 4px; ",
                              "background-color: ", warning_color, "22; ",
                              "border: 1px solid ", warning_color, ";"),
                p(style = paste0("margin: 0; color: ", warning_color, "; font-size: 12px;"),
                  rv$warning_message)
            )
        }
    })
    
    # Validate FDR thresholds
    get_valid_fdr_threshold <- function(value, default) {
        if (is.null(value) || value == "") return(default)
        num_val <- suppressWarnings(as.numeric(value))
        if (is.na(num_val) || num_val < 0 || num_val > 1) return(default)
        return(num_val)
    }
    
    # Helper to get significant genes based on FDR thresholds
    get_significant_genes <- function(gene_data, result_type, thresh1, thresh2, use_thresh1, use_thresh2, top_n = 10, rank_by_rra = FALSE, reverse_rank = FALSE) {
        # DEBUG: Trace function inputs - commented out to avoid console spam
        # cat("DEBUG get_significant_genes: rank_by_rra=", rank_by_rra, ", reverse_rank=", reverse_rank, "\n", file=stderr())
        
        if (result_type == "neg") {
            fdr_col <- "neg_fdr"
            lfc_col <- "neg_lfc"
            rank_col <- "neg_rank"
        } else {
            fdr_col <- "pos_fdr"
            lfc_col <- "pos_lfc"
            rank_col <- "pos_rank"
        }
        
        colnames(gene_data) <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(gene_data)))
        fdr_col <- tolower(fdr_col)
        lfc_col <- tolower(lfc_col)
        rank_col <- tolower(rank_col)
        
        gene_id_col <- if ("id" %in% colnames(gene_data)) "id" 
                      else if ("gene" %in% colnames(gene_data)) "gene" 
                      else colnames(gene_data)[1]
        
        sig_groups <- list()
        
        if (use_thresh1) {
            ordered_t1 <- gene_data %>%
                filter(!is.na(.data[[fdr_col]]) & .data[[fdr_col]] < thresh1)
            # Primary sort: smallest FDR; Tie-break: RRA rank or LFC direction
            if (rank_by_rra && rank_col %in% colnames(ordered_t1)) {
                if (reverse_rank) {
                    ordered_t1 <- ordered_t1 %>% arrange(.data[[fdr_col]], desc(.data[[rank_col]]))
                } else {
                    ordered_t1 <- ordered_t1 %>% arrange(.data[[fdr_col]], .data[[rank_col]])
                }
            } else {
                if (result_type == "neg") {
                    if (reverse_rank) {
                        ordered_t1 <- ordered_t1 %>% arrange(.data[[fdr_col]], desc(.data[[lfc_col]]))
                    } else {
                        ordered_t1 <- ordered_t1 %>% arrange(.data[[fdr_col]], .data[[lfc_col]])
                    }
                } else {
                    if (reverse_rank) {
                        ordered_t1 <- ordered_t1 %>% arrange(.data[[fdr_col]], .data[[lfc_col]])
                    } else {
                        ordered_t1 <- ordered_t1 %>% arrange(.data[[fdr_col]], desc(.data[[lfc_col]]))
                    }
                }
            }
            genes_t1 <- ordered_t1 %>%
                head(top_n) %>%
                pull(.data[[gene_id_col]])
            if (length(genes_t1) > 0) {
                sig_groups <- c(sig_groups, list(list(genes = genes_t1, threshold = thresh1)))
            }
        }
        
        if (use_thresh2) {
            ordered_t2 <- gene_data %>%
                filter(!is.na(.data[[fdr_col]]) & .data[[fdr_col]] >= thresh1 & .data[[fdr_col]] < thresh2)
            if (rank_by_rra && rank_col %in% colnames(ordered_t2)) {
                if (reverse_rank) {
                    ordered_t2 <- ordered_t2 %>% arrange(.data[[fdr_col]], desc(.data[[rank_col]]))
                } else {
                    ordered_t2 <- ordered_t2 %>% arrange(.data[[fdr_col]], .data[[rank_col]])
                }
            } else {
                if (result_type == "neg") {
                    if (reverse_rank) {
                        ordered_t2 <- ordered_t2 %>% arrange(.data[[fdr_col]], desc(.data[[lfc_col]]))
                    } else {
                        ordered_t2 <- ordered_t2 %>% arrange(.data[[fdr_col]], .data[[lfc_col]])
                    }
                } else {
                    if (reverse_rank) {
                        ordered_t2 <- ordered_t2 %>% arrange(.data[[fdr_col]], .data[[lfc_col]])
                    } else {
                        ordered_t2 <- ordered_t2 %>% arrange(.data[[fdr_col]], desc(.data[[lfc_col]]))
                    }
                }
            }
            genes_t2 <- ordered_t2 %>%
                head(top_n) %>%
                pull(.data[[gene_id_col]])
            if (length(genes_t2) > 0) {
                sig_groups <- c(sig_groups, list(list(genes = genes_t2, threshold = thresh2)))
            }
        }
        
        return(sig_groups)
    }
    
    # Build highlight config based on selections
    # Calculate significant genes separately - this runs independently
    significant_genes_data <- reactive({
        req(input$highlight_significant, input$selection_type, rv$gene_data)
        
        use_plot_thresh <- if (is.null(input$use_plot_thresholds)) TRUE else input$use_plot_thresholds
        
        if (use_plot_thresh) {
            # Use thresholds from Plot Options
            thresh1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05)
            thresh2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
            use_t1 <- if (is.null(input$sig_use_thresh1)) TRUE else input$sig_use_thresh1
            use_t2 <- if (is.null(input$sig_use_thresh2)) TRUE else input$sig_use_thresh2
            color1 <- if (!is.null(input$sig_color_1_plot)) input$sig_color_1_plot else rv$sig_color_1
            color2 <- if (!is.null(input$sig_color_2_plot)) input$sig_color_2_plot else rv$sig_color_2
        } else {
            # Use custom thresholds
            thresh1 <- get_valid_fdr_threshold(input$sig_fdr_threshold_1, 0.05)
            thresh2 <- get_valid_fdr_threshold(input$sig_fdr_threshold_2, 0.1)
            use_t1 <- TRUE
            use_t2 <- TRUE
            color1 <- if (!is.null(input$sig_color_1)) input$sig_color_1 else rv$sig_color_1
            color2 <- if (!is.null(input$sig_color_2)) input$sig_color_2 else rv$sig_color_2
        }
        
        top_n <- if (is.null(input$sig_top_n) || is.na(input$sig_top_n)) 10 else as.integer(input$sig_top_n)
        rank_by_rra <- isTRUE(input$rank_method == "rra")
        reverse_rank <- isTRUE(input$reverse_rank_order)
        
        sig_groups <- get_significant_genes(rv$gene_data, input$selection_type, thresh1, thresh2, use_t1, use_t2, top_n, rank_by_rra, reverse_rank)
        
        list(sig_groups = sig_groups, thresh2 = thresh2, color1 = color1, color2 = color2)
    })
    
    # Update warning message whenever significant genes change
    observe({
        req(input$highlight_significant)
        if (input$highlight_significant) {
            sig_data <- significant_genes_data()
            highlight_config <- get_highlight_config()
            
            if ((is.null(highlight_config) || length(highlight_config) == 0) && length(sig_data$sig_groups) == 0) {
                rv$sig_warning <- paste0("No genes with FDR < ", sig_data$thresh2, " for ", if (input$selection_type == "neg") "negative" else "positive", " selection.")
            } else {
                rv$sig_warning <- NULL
            }
        } else {
            rv$sig_warning <- NULL
        }
    })
    
    get_highlight_config <- reactive({
        req(input$comparison, rv$result_files)
        # Access rank inputs to establish dependency, but don't block on them
        rank_method <- input$rank_method
        reverse_rank <- input$reverse_rank_order
        
        result_info <- rv$result_files[[input$comparison]]
        highlight_config <- list()
        
        # Add enabled predefined gene sets with user-modified names and colors
        if (!is.null(result_info$plot_config) && 
            !is.null(result_info$plot_config$highlight_genes)) {
            
            gene_sets <- result_info$plot_config$highlight_genes
            for (i in seq_along(gene_sets)) {
                checkbox_id <- paste0("geneset_", i)
                if (!is.null(input[[checkbox_id]]) && input[[checkbox_id]]) {
                    # Get user-modified name and color
                    custom_name <- input[[paste0("geneset_name_", i)]]
                    custom_color <- input[[paste0("geneset_color_", i)]]
                    
                    modified_set <- gene_sets[[i]]
                    if (!is.null(custom_name) && custom_name != "") {
                        modified_set$group_name <- custom_name
                    }
                    if (!is.null(custom_color)) {
                        modified_set$color <- custom_color
                    }
                    
                    highlight_config <- c(highlight_config, list(modified_set))
                }
            }
        }
        
        # Add custom gene list if present and show checkbox is enabled
        if (!is.null(rv$custom_gene_list) && !is.null(input$show_custom) && input$show_custom) {
            custom_config <- list(
                group_name = rv$custom_gene_list$name,
                color = rv$custom_gene_list$color,
                genes = rv$custom_gene_list$genes,
                is_custom = TRUE
            )
            highlight_config <- c(highlight_config, list(custom_config))
        }
        
        # Add significant genes if enabled
        if (!is.null(input$highlight_significant) && input$highlight_significant && !is.null(rv$gene_data)) {
            sig_data <- significant_genes_data()
            sig_groups <- sig_data$sig_groups
            
            for (i in seq_along(sig_groups)) {
                group <- sig_groups[[i]]
                sig_config <- list(
                    group_name = paste0("FDR < ", group$threshold),
                    color = if (i == 1) sig_data$color1 else sig_data$color2,
                    genes = group$genes,
                    is_custom = TRUE
                )
                highlight_config <- c(highlight_config, list(sig_config))
            }
        }
        
        return(if (length(highlight_config) > 0) highlight_config else NULL)
    })

    # Render significant warning UI
    output$global_sig_warning <- renderUI({
        if (isTRUE(input$highlight_significant) && !is.null(rv$sig_warning) && nzchar(rv$sig_warning)) {
            tags$div(style = "margin-top: 6px; padding: 6px; background-color: #fff3cd; color: #856404; border: 1px solid #ffeeba; border-radius: 4px; font-size: 0.9em;",
                rv$sig_warning
            )
        }
    })
    outputOptions(output, "global_sig_warning", suspendWhenHidden = FALSE)
    
    # Create volcano plot
    create_interactive_volcano <- reactive({
        req(rv$gene_data, input$selection_type)
        # Access rank inputs to establish dependency, but don't block on them
        rank_method <- input$rank_method
        reverse_rank <- input$reverse_rank_order
        
        result_info <- rv$result_files[[input$comparison]]
        comparison_name <- result_info$comparison_name
        
        fdr1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05)
        fdr2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
        
        plot <- create_volcano_plot_interactive(
            rv$gene_data,
            comparison_name,
            input$selection_type,
            get_highlight_config(),
            rv$config_dir,
            fdr1,
            fdr2,
            input$show_labels,
            input$label_scale,
            text_scale = input$text_scale,
            show_fdr_line_1 = input$show_fdr_line_1,
            show_fdr_line_2 = input$show_fdr_line_2,
            custom_title = input$volcano_title,
            point_size = input$point_size,
            fdr_color_1 = get_fdr_color_1(),
            fdr_color_2 = get_fdr_color_2(),
            plot_metric = input$plot_metric
        )
        
        return(plot)
    })
    
    # Create ranked plot
    create_interactive_ranked <- reactive({
        req(rv$gene_data, input$selection_type)
        # Access rank inputs to establish dependency, but don't block on them
        rank_method <- input$rank_method
        reverse_rank <- input$reverse_rank_order
        
        result_info <- rv$result_files[[input$comparison]]
        comparison_name <- result_info$comparison_name
        
        fdr1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05)
        fdr2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
        
        plot <- create_ranked_plot_interactive(
            rv$gene_data,
            comparison_name,
            input$selection_type,
            get_highlight_config(),
            rv$config_dir,
            fdr1,
            fdr2,
            input$show_labels,
            input$label_scale,
            scale_point_size = input$scale_point_size,
            color_by_fdr = input$color_by_fdr,
            text_scale = input$text_scale,
            custom_title = input$ranked_title,
            point_size = input$point_size,
            fdr_color_1 = get_fdr_color_1(),
            fdr_color_2 = get_fdr_color_2(),
            rank_method = input$rank_method,
            reverse_rank = input$reverse_rank_order,
            sig_use_plot_thresholds = if (is.null(input$use_plot_thresholds)) TRUE else input$use_plot_thresholds,
            plot_metric = input$plot_metric
        )
        
        return(plot)
    })
    
    # Render plots
    output$volcano_plot <- {
        if (PLOTLY_AVAILABLE) {
            plotly::renderPlotly({
                plot_obj <- create_interactive_volcano()
                use_plotly <- if (is.null(input$use_interactive_plots)) TRUE else input$use_interactive_plots
                render_plot_with_fallback(plot_obj, use_plotly = use_plotly, tooltip = "text")
            })
        } else {
            renderPlot({
                create_interactive_volcano()
            })
        }
    }
    
    output$volcano_plot_static <- renderPlot({
        create_interactive_volcano()
    })
    
    output$ranked_plot <- {
        if (PLOTLY_AVAILABLE) {
            plotly::renderPlotly({
                plot_obj <- create_interactive_ranked()
                use_plotly <- if (is.null(input$use_interactive_plots)) TRUE else input$use_interactive_plots
                render_plot_with_fallback(plot_obj, use_plotly = use_plotly, tooltip = "text")
            })
        } else {
            renderPlot({
                create_interactive_ranked()
            })
        }
    }
    
    output$ranked_plot_static <- renderPlot({
        create_interactive_ranked()
    })
    
    # Download handlers
    output$download_volcano <- downloadHandler(
        filename = function() {
            tryCatch({
                cat("[VOLCANO_DOWNLOAD] Filename function called\n")
                ext <- tolower(input$download_format)
                ext <- if (ext %in% c("png", "pdf", "html")) ext else "png"
                comparison <- if (!is.null(input$comparison)) input$comparison else "comparison"
                selection <- if (!is.null(input$selection_type)) input$selection_type else "neg"
                fname <- paste0(comparison, "_", selection, "_volcano_", 
                       format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
                cat(paste0("[VOLCANO_DOWNLOAD] Filename: ", fname, "\n"))
                return(fname)
            }, error = function(e) {
                cat(paste0("[VOLCANO_DOWNLOAD] ERROR in filename function: ", e$message, "\n"))
                return(paste0("volcano_plot.", tolower(input$download_format)))
            })
        },
        content = function(file) {
            cat(paste0("[VOLCANO_DOWNLOAD] Content function called, file=", file, "\n"))
            req(rv$gene_data, input$selection_type, input$comparison)
            cat("[VOLCANO_DOWNLOAD] req() passed\n")
            
            result_info <- rv$result_files[[input$comparison]]
            comparison_name <- result_info$comparison_name
            
            fdr1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05)
            fdr2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
            
            # Create plot directly instead of using reactive
            plot <- create_volcano_plot_interactive(
                rv$gene_data,
                comparison_name,
                input$selection_type,
                get_highlight_config(),
                rv$config_dir,
                fdr1,
                fdr2,
                input$show_labels,
                input$label_scale,
                text_scale = input$text_scale,
                show_fdr_line_1 = input$show_fdr_line_1,
                show_fdr_line_2 = input$show_fdr_line_2,
                custom_title = input$volcano_title,
                point_size = input$point_size,
                fdr_color_1 = get_fdr_color_1(),
                fdr_color_2 = get_fdr_color_2(),
                plot_metric = input$plot_metric
            )
            
            cat("[VOLCANO_DOWNLOAD] Plot created successfully\n")
            ext <- tolower(input$download_format)
            cat(paste0("[VOLCANO_DOWNLOAD] Format: ", ext, ", width=", input$download_width, ", height=", input$download_height, "\n"))
            
            if (ext == "html") {
                if (PLOTLY_AVAILABLE) {
                    cat("[VOLCANO_DOWNLOAD] Saving as HTML (Interactive)\n")
                    plotly_obj <- ggplotly(plot, tooltip = "text") %>% layout(hovermode = "closest")
                    htmlwidgets::saveWidget(plotly_obj, file = file)
                    cat("[VOLCANO_DOWNLOAD] HTML saved successfully\n")
                } else {
                    cat("[VOLCANO_DOWNLOAD] ERROR: HTML export requested but plotly not available, falling back to PNG\n")
                    ext <- "png"
                }
            }
            
            if (ext == "pdf") {
                cat("[VOLCANO_DOWNLOAD] Saving as PDF\n")
                if (capabilities("cairo")) {
                    cat("[VOLCANO_DOWNLOAD] Using cairo_pdf\n")
                    ggsave(filename = file, plot = plot,
                           width = input$download_width,
                           height = input$download_height,
                           units = "in", device = grDevices::cairo_pdf, bg = "white")
                    cat("[VOLCANO_DOWNLOAD] PDF saved successfully\n")
                } else {
                    cat("[VOLCANO_DOWNLOAD] Using pdf() with useDingbats=FALSE\n")
                    ggsave(filename = file, plot = plot,
                           width = input$download_width,
                           height = input$download_height,
                           units = "in",
                           device = function(...) grDevices::pdf(..., useDingbats = FALSE),
                           bg = "white")
                    cat("[VOLCANO_DOWNLOAD] PDF saved successfully\n")
                }
            } else {
                cat("[VOLCANO_DOWNLOAD] Saving as PNG\n")
                ggsave(filename = file, plot = plot,
                       width = input$download_width,
                       height = input$download_height,
                       units = "in", device = "png",
                       dpi = as.numeric(input$download_dpi), bg = "white")
                cat("[VOLCANO_DOWNLOAD] PNG saved successfully\n")
            }
        }
    )
    
    output$download_ranked <- downloadHandler(
        filename = function() {
            tryCatch({
                cat("[RANKED_DOWNLOAD] Filename function called\n")
                ext <- tolower(input$download_format)
                ext <- if (ext %in% c("png", "pdf", "html")) ext else "png"
                comparison <- if (!is.null(input$comparison)) input$comparison else "comparison"
                selection <- if (!is.null(input$selection_type)) input$selection_type else "neg"
                fname <- paste0(comparison, "_", selection, "_ranked_", 
                       format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
                cat(paste0("[RANKED_DOWNLOAD] Filename: ", fname, "\n"))
                return(fname)
            }, error = function(e) {
                cat(paste0("[RANKED_DOWNLOAD] ERROR in filename function: ", e$message, "\n"))
                return(paste0("ranked_plot.", tolower(input$download_format)))
            })
        },
        content = function(file) {
            cat(paste0("[RANKED_DOWNLOAD] Content function called, file=", file, "\n"))
            req(rv$gene_data, input$selection_type, input$comparison)
            cat("[RANKED_DOWNLOAD] req() passed\n")
            
            result_info <- rv$result_files[[input$comparison]]
            comparison_name <- result_info$comparison_name
            
            fdr1 <- get_valid_fdr_threshold(input$fdr_threshold_1, 0.05)
            fdr2 <- get_valid_fdr_threshold(input$fdr_threshold_2, 0.1)
            
            # Create plot directly instead of using reactive
            plot <- create_ranked_plot_interactive(
                rv$gene_data,
                comparison_name,
                input$selection_type,
                get_highlight_config(),
                rv$config_dir,
                fdr1,
                fdr2,
                input$show_labels,
                input$label_scale,
                scale_point_size = input$scale_point_size,
                color_by_fdr = input$color_by_fdr,
                text_scale = input$text_scale,
                custom_title = input$ranked_title,
                point_size = input$point_size,
                fdr_color_1 = get_fdr_color_1(),
                fdr_color_2 = get_fdr_color_2(),
                rank_method = input$rank_method,
                reverse_rank = input$reverse_rank_order,
                sig_use_plot_thresholds = if (is.null(input$use_plot_thresholds)) TRUE else input$use_plot_thresholds,
                plot_metric = input$plot_metric
            )
            
            cat("[RANKED_DOWNLOAD] Ranked plot created successfully\n")
            ext <- tolower(input$download_format)
            cat(paste0("[RANKED_DOWNLOAD] Format: ", ext, ", width=", input$download_width, ", height=", input$download_height, "\n"))
            
            if (ext == "html") {
                if (PLOTLY_AVAILABLE) {
                    cat("[RANKED_DOWNLOAD] Saving as HTML (Interactive)\n")
                    plotly_obj <- ggplotly(plot, tooltip = "text") %>% layout(hovermode = "closest")
                    htmlwidgets::saveWidget(plotly_obj, file = file)
                    cat("[RANKED_DOWNLOAD] HTML saved successfully\n")
                } else {
                    cat("[RANKED_DOWNLOAD] ERROR: HTML export requested but plotly not available, falling back to PNG\n")
                    ext <- "png"
                }
            }
            
            if (ext == "pdf") {
                cat("[RANKED_DOWNLOAD] Saving as PDF\n")
                if (capabilities("cairo")) {
                    cat("[RANKED_DOWNLOAD] Using cairo_pdf\n")
                    ggsave(filename = file, plot = plot,
                           width = input$download_width,
                           height = input$download_height,
                           units = "in", device = grDevices::cairo_pdf, bg = "white")
                    cat("[RANKED_DOWNLOAD] PDF saved successfully\n")
                } else {
                    cat("[RANKED_DOWNLOAD] Using pdf() with useDingbats=FALSE\n")
                    ggsave(filename = file, plot = plot,
                           width = input$download_width,
                           height = input$download_height,
                           units = "in",
                           device = function(...) grDevices::pdf(..., useDingbats = FALSE),
                           bg = "white")
                    cat("[RANKED_DOWNLOAD] PDF saved successfully\n")
                }
            } else {
                cat("[RANKED_DOWNLOAD] Saving as PNG\n")
                ggsave(filename = file, plot = plot,
                       width = input$download_width,
                       height = input$download_height,
                       units = "in", device = "png",
                       dpi = as.numeric(input$download_dpi), bg = "white")
                cat("[RANKED_DOWNLOAD] PNG saved successfully\n")
            }
        }
    )
}

# Interactive plot functions with custom parameters
create_volcano_plot_interactive <- function(gene_data, comparison_name, result_type, 
                                           highlight_config, config_dir, fdr_thresh1, fdr_thresh2,
                                           show_labels, label_scale, text_scale = 1,
                                           show_fdr_line_1 = TRUE, show_fdr_line_2 = TRUE,
                                           custom_title = "", point_size = 2,
                                           fdr_color_1 = "#FF0000", fdr_color_2 = "#FFA500",
                                           plot_metric = "lfc") {
    
    if (result_type == "neg") {
        lfc_col <- "neg_lfc"
        fdr_col <- "neg_fdr"
        title_suffix <- "Negative Selection"
    } else {
        lfc_col <- "pos_lfc"
        fdr_col <- "pos_fdr"
        title_suffix <- "Positive Selection"
    }
    score_col <- if (result_type == "neg") "neg_score" else "pos_score"
    
    colnames(gene_data) <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(gene_data)))
    lfc_col <- tolower(lfc_col)
    fdr_col <- tolower(fdr_col)
    score_col <- tolower(score_col)
    
    gene_id_col <- if ("id" %in% colnames(gene_data)) "id" 
                   else if ("gene" %in% colnames(gene_data)) "gene" 
                   else colnames(gene_data)[1]
    
    plot_data <- gene_data %>%
        filter(!is.na(.data[[lfc_col]]) & !is.na(.data[[fdr_col]]))
    
    min_nonzero_fdr <- min(plot_data[[fdr_col]][plot_data[[fdr_col]] > 0], na.rm = TRUE)
    replacement_fdr <- min(1e-200, min_nonzero_fdr / 10)
    
    plot_data <- plot_data %>%
        mutate(
            adjusted_fdr = ifelse(.data[[fdr_col]] == 0, replacement_fdr, .data[[fdr_col]]),
            neg_log10_fdr = -log10(adjusted_fdr),
            highlight_group = "none",
            highlight_color = "gray50"
        )

    # Compute metric for x-axis
    if (plot_metric == "score") {
        if (!score_col %in% colnames(plot_data)) {
            stop(paste0("Missing required score column: ", score_col))
        }
        min_nonzero_score <- min(plot_data[[score_col]][plot_data[[score_col]] > 0], na.rm = TRUE)
        replacement_score <- min(1e-200, min_nonzero_score / 10)
        plot_data <- plot_data %>%
            mutate(score_value = -log10(ifelse(.data[[score_col]] == 0, replacement_score, .data[[score_col]])))
        x_col <- "score_value"
        x_label <- "-log10(Score)"
        metric_label <- "-log10(Score)"
    } else {
        x_col <- lfc_col
        x_label <- "Log2 Fold Change"
        metric_label <- "LFC"
    }
    
    genes_to_label <- NULL
    legend_data <- NULL
    
    if (!is.null(highlight_config) && length(highlight_config) > 0) {
        for (i in seq_along(highlight_config)) {
            group_info <- highlight_config[[i]]
            group_color <- group_info$color %||% "black"
            group_name <- group_info$group_name %||% paste0("Group ", i)
            
            # Load from file or use provided list
            if (!is.null(group_info$is_custom) && group_info$is_custom) {
                highlight_genes <- group_info$genes
            } else {
                gene_file <- group_info$file
                                # Resolve gene file path relative to config directory
                                if (!is.null(config_dir) && !file.exists(gene_file)) {
                                    gene_file <- file.path(config_dir, gene_file)
                                }
                highlight_genes <- load_gene_list(gene_file)
            }
            
            if (length(highlight_genes) > 0) {
                plot_data <- plot_data %>%
                    mutate(
                        highlight_group = ifelse(
                            .data[[gene_id_col]] %in% highlight_genes & highlight_group == "none",
                            group_name,
                            highlight_group
                        ),
                        highlight_color = ifelse(
                            .data[[gene_id_col]] %in% highlight_genes & highlight_color == "gray50",
                            group_color,
                            highlight_color
                        )
                    )
            }
        }
        
        genes_to_label <- plot_data %>% filter(highlight_group != "none")
        legend_data <- tibble(
            group = sapply(highlight_config, function(x) x$group_name %||% "Unknown"),
            color = sapply(highlight_config, function(x) x$color %||% "black")
        )
    }
    
    if (plot_metric == "score") {
        max_val <- max(plot_data[[x_col]], na.rm = TRUE)
        x_limit <- ceiling(max_val * 10) / 10
    } else {
        max_abs_val <- max(abs(plot_data[[x_col]]), na.rm = TRUE)
        x_limit <- ceiling(max_abs_val * 10) / 10
    }
    
    # Build FDR threshold legend data
    fdr_legend_data <- tibble()
    if (show_fdr_line_1) {
        fdr_legend_data <- bind_rows(fdr_legend_data, 
            tibble(label = paste0("FDR = ", fdr_thresh1), color = fdr_color_1))
    }
    if (show_fdr_line_2) {
        fdr_legend_data <- bind_rows(fdr_legend_data, 
            tibble(label = paste0("FDR = ", fdr_thresh2), color = fdr_color_2))
    }
    
    p <- ggplot(plot_data, aes(x = .data[[x_col]], y = neg_log10_fdr))
    
    # Add thicker vertical line at LFC=0 to show transition from negative to positive
    if (plot_metric == "lfc") {
        p <- p + geom_vline(xintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7)
    }
    
    if (show_fdr_line_1) {
        p <- p + geom_hline(yintercept = -log10(fdr_thresh1), linetype = "dashed", 
                           color = fdr_color_1, linewidth = 0.5)
    }
    if (show_fdr_line_2) {
        p <- p + geom_hline(yintercept = -log10(fdr_thresh2), linetype = "dashed", 
                           color = fdr_color_2, linewidth = 0.5)
    }
    
    p <- p + geom_point(
            data = plot_data %>% filter(highlight_group == "none"),
            aes(text = paste0(.data[[gene_id_col]], "<br>", metric_label, ": ", round(.data[[x_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))),
            color = "gray50", size = point_size * 0.75, alpha = 0.5
        ) +
        geom_point(
            data = plot_data %>% filter(highlight_group != "none"),
            aes(color = highlight_group, text = paste0(.data[[gene_id_col]], "<br>", metric_label, ": ", round(.data[[x_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))),
            size = point_size, alpha = 0.8
        )
    
    if (!is.null(legend_data)) {
        color_values <- setNames(legend_data$color, legend_data$group)
        p <- p + scale_color_manual(name = "Gene Group", values = color_values) +
            guides(color = guide_legend(override.aes = list(size = 4)))
        
        if (show_labels && !is.null(genes_to_label) && nrow(genes_to_label) > 0) {
            p <- p + geom_text_repel(
                data = genes_to_label,
                aes(label = .data[[gene_id_col]], color = highlight_group),
                size = 3 * label_scale,
                max.overlaps = Inf,
                box.padding = 0.8,
                point.padding = 0.8,
                min.segment.length = 0,
                force = 3,
                force_pull = 1,
                fontface = "bold",
                show.legend = FALSE
            )
        }
    }
    
    # Add FDR threshold annotations as text labels on the plot
    if (nrow(fdr_legend_data) > 0) {
        # Add text annotations for FDR thresholds on the right side of plot
        y_positions <- c()
        if (show_fdr_line_1) y_positions <- c(y_positions, -log10(fdr_thresh1))
        if (show_fdr_line_2) y_positions <- c(y_positions, -log10(fdr_thresh2))
        
        fdr_annotations <- tibble(
            x = rep(x_limit * 0.95, length(y_positions)),
            y = y_positions,
            label = fdr_legend_data$label,
            color = fdr_legend_data$color
        )
        
        p <- p + geom_text(
            data = fdr_annotations,
            aes(x = x, y = y, label = label),
            color = fdr_annotations$color,
            hjust = 1, vjust = -0.5,
            size = 3 * text_scale,
            fontface = "bold"
        )
    }
    
    # Format title
    plot_title <- if (custom_title != "" && !is.null(custom_title)) {
        custom_title
    } else {
        clean_name <- gsub("_", " ", comparison_name)
        paste0(clean_name, " (", title_suffix, ")")
    }
    
    p <- p +
        labs(
            title = plot_title,
            x = x_label,
            y = "-log10(FDR)"
        ) +
        (
            if (plot_metric == "score") {
                scale_x_continuous(limits = c(0, x_limit))
            } else {
                scale_x_continuous(limits = c(-x_limit, x_limit))
            }
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14 * text_scale),
            axis.title = element_text(size = 11 * text_scale),
            axis.text = element_text(size = 10 * text_scale),
            legend.title = element_text(size = 11 * text_scale),
            legend.text = element_text(size = 10 * text_scale),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            legend.position = "right"
        )
    
    return(p)
}

create_ranked_plot_interactive <- function(gene_data, comparison_name, result_type, 
                                          highlight_config, config_dir, fdr_thresh1, fdr_thresh2,
                                          show_labels, label_scale, scale_point_size = TRUE, color_by_fdr = TRUE, text_scale = 1,
                                          custom_title = "", point_size = 2,
                                          fdr_color_1 = "#FF0000", fdr_color_2 = "#FFA500",
                                          sig_use_plot_thresholds = TRUE, rank_method = "lfc", reverse_rank = FALSE,
                                          plot_metric = "lfc") {
    
    if (result_type == "neg") {
        lfc_col <- "neg_lfc"
        fdr_col <- "neg_fdr"
        title_suffix <- "Negative Selection"
    } else {
        lfc_col <- "pos_lfc"
        fdr_col <- "pos_fdr"
        title_suffix <- "Positive Selection"
    }
    score_col <- if (result_type == "neg") "neg_score" else "pos_score"
    
    colnames(gene_data) <- tolower(gsub("[^a-zA-Z0-9]", "_", colnames(gene_data)))
    lfc_col <- tolower(lfc_col)
    fdr_col <- tolower(fdr_col)
    score_col <- tolower(score_col)
    
    # Look for RRA rank column (could be neg_rank or pos_rank)
    rank_col <- if (result_type == "neg") "neg_rank" else "pos_rank"
    rank_col <- tolower(rank_col)
    
    gene_id_col <- if ("id" %in% colnames(gene_data)) "id" 
                   else if ("gene" %in% colnames(gene_data)) "gene" 
                   else colnames(gene_data)[1]
    
    plot_data <- gene_data %>%
        filter(!is.na(.data[[lfc_col]]) & !is.na(.data[[fdr_col]])) %>%
        mutate(
            significance = case_when(
                .data[[fdr_col]] < fdr_thresh1 ~ fdr_color_1,
                .data[[fdr_col]] >= fdr_thresh1 & .data[[fdr_col]] < fdr_thresh2 ~ fdr_color_2,
                TRUE ~ "black"
            ),
            point_size = case_when(
                .data[[fdr_col]] < 0.0001 ~ 4,
                .data[[fdr_col]] < 0.001 ~ 3,
                .data[[fdr_col]] < 0.01 ~ 2.5,
                .data[[fdr_col]] < fdr_thresh1 ~ 2,
                .data[[fdr_col]] < fdr_thresh2 ~ 1.5,
                TRUE ~ 1
            ),
            highlight_group = "none",
            highlight_color = NA_character_
        )

    # Compute metric for y-axis
    if (plot_metric == "score") {
        if (!score_col %in% colnames(plot_data)) {
            stop(paste0("Missing required score column: ", score_col))
        }
        min_nonzero_score <- min(plot_data[[score_col]][plot_data[[score_col]] > 0], na.rm = TRUE)
        replacement_score <- min(1e-200, min_nonzero_score / 10)
        plot_data <- plot_data %>%
            mutate(score_value = -log10(ifelse(.data[[score_col]] == 0, replacement_score, .data[[score_col]])))
        y_col <- "score_value"
        y_label <- "-log10(Score)"
        metric_label <- "-log10(Score)"
    } else {
        y_col <- lfc_col
        y_label <- "Log2 Fold Change"
        metric_label <- "LFC"
    }
    
    # Apply ranking method (always sort in same direction, scale reversal will flip display)
    rank_col <- tolower(rank_col)  # Ensure normalized
    if (rank_method == "rra" && rank_col %in% colnames(plot_data)) {
        # Sort by RRA rank
        plot_data <- plot_data %>% arrange(.data[[rank_col]])
    } else {
        # Sort by LFC (default)
        plot_data <- plot_data %>% arrange(.data[[lfc_col]])
    }
    
    plot_data <- plot_data %>% mutate(rank = row_number())
    
    genes_to_label <- NULL
    legend_data <- NULL
    
    if (!is.null(highlight_config) && length(highlight_config) > 0) {
        for (i in seq_along(highlight_config)) {
            group_info <- highlight_config[[i]]
            group_color <- group_info$color %||% "black"
            group_name <- group_info$group_name %||% paste0("Group ", i)
            
            if (!is.null(group_info$is_custom) && group_info$is_custom) {
                highlight_genes <- group_info$genes
            } else {
                gene_file <- group_info$file
                if (!is.null(config_dir) && !file.exists(gene_file)) {
                    gene_file <- file.path(config_dir, gene_file)
                }
                highlight_genes <- load_gene_list(gene_file)
            }
            
            if (length(highlight_genes) > 0) {
                plot_data <- plot_data %>%
                    mutate(
                        highlight_group = ifelse(
                            .data[[gene_id_col]] %in% highlight_genes & highlight_group == "none",
                            group_name,
                            highlight_group
                        ),
                        highlight_color = ifelse(
                            .data[[gene_id_col]] %in% highlight_genes & is.na(highlight_color),
                            group_color,
                            highlight_color
                        )
                    )
            }
        }
        
        genes_to_label <- plot_data %>% filter(highlight_group != "none")
        legend_data <- tibble(
            group = sapply(highlight_config, function(x) x$group_name %||% "Unknown"),
            color = sapply(highlight_config, function(x) x$color %||% "black")
        )
    }
    
    plot_data <- plot_data %>%
        mutate(
            fdr_category = factor(
                case_when(
                    significance == fdr_color_1 ~ paste0("FDR < ", fdr_thresh1),
                    significance == fdr_color_2 ~ paste0(fdr_thresh1, " ≤ FDR < ", fdr_thresh2),
                    TRUE ~ paste0("FDR ≥ ", fdr_thresh2)
                ),
                levels = c(paste0("FDR < ", fdr_thresh1), 
                          paste0(fdr_thresh1, " ≤ FDR < ", fdr_thresh2),
                          paste0("FDR ≥ ", fdr_thresh2))
            ),
            outline_group = ifelse(highlight_group == "none", "None", highlight_group)
        )
    
    # Build FDR bin (for size legend)
    plot_data <- plot_data %>% mutate(
        fdr_bin = factor(case_when(
            .data[[fdr_col]] < 1e-5 ~ "<1e-5",
            .data[[fdr_col]] < 1e-4 ~ "<1e-4",
            .data[[fdr_col]] < 1e-3 ~ "< 0.001",
            .data[[fdr_col]] < 1e-2 ~ "< 0.01",
            .data[[fdr_col]] < fdr_thresh1 ~ paste0("< ", fdr_thresh1),
            .data[[fdr_col]] < fdr_thresh2 ~ paste0("< ", fdr_thresh2),
            TRUE ~ paste0("\u2265 ", fdr_thresh2)
        ), levels = c("<1e-5","<1e-4","< 0.001","< 0.01", paste0("< ", fdr_thresh1), paste0("< ", fdr_thresh2), paste0("\u2265 ", fdr_thresh2)))
    )

    # Build size values map
    size_values <- setNames(
        c(4.5, 4, 3.5, 3, 2.5, 2, 1.5),
        c("<1e-5", "<1e-4", "< 0.001", "< 0.01", paste0("< ", fdr_thresh1), paste0("< ", fdr_thresh2), paste0("\u2265 ", fdr_thresh2))
    )

    if (color_by_fdr) {
        zero_line <- if (plot_metric == "lfc") geom_hline(yintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7) else NULL
        if (scale_point_size) {
            # Scale size_values by point_size
            scaled_size_values <- size_values * (point_size / 2)
            p <- ggplot(plot_data, aes(x = rank, y = .data[[y_col]])) +
                zero_line +
                geom_point(
                    data = plot_data %>% filter(highlight_group == "none"),
                    aes(color = fdr_category, size = fdr_bin, text = paste0(.data[[gene_id_col]], "<br>Rank: ", rank, "<br>", metric_label, ": ", round(.data[[y_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))),
                    alpha = 0.7
                ) +
                geom_point(
                    data = plot_data %>% filter(highlight_group != "none"),
                    aes(
                        fill = if (sig_use_plot_thresholds) fdr_category else highlight_group,
                        color = outline_group, size = fdr_bin, text = paste0(.data[[gene_id_col]], "<br>Rank: ", rank, "<br>", metric_label, ": ", round(.data[[y_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))
                    ),
                    shape = 21, stroke = 1.2, alpha = 0.8
                ) +
                scale_size_manual(name = "FDR (size)", values = scaled_size_values)
        } else {
            p <- ggplot(plot_data, aes(x = rank, y = .data[[y_col]])) +
                geom_point(
                    data = plot_data %>% filter(highlight_group == "none"),
                    aes(color = fdr_category, text = paste0(.data[[gene_id_col]], "<br>Rank: ", rank, "<br>", metric_label, ": ", round(.data[[y_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))),
                    size = point_size, alpha = 0.7
                ) +
                geom_point(
                    data = plot_data %>% filter(highlight_group != "none"),
                    aes(
                        fill = if (sig_use_plot_thresholds) fdr_category else highlight_group,
                        color = outline_group, text = paste0(.data[[gene_id_col]], "<br>Rank: ", rank, "<br>", metric_label, ": ", round(.data[[y_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))
                    ),
                    shape = 21, stroke = 1.2, size = point_size, alpha = 0.8
                )
        }
        # Combined color scale includes categories (not shown) + highlight outlines (shown)
        highlight_colors <- if (!is.null(highlight_config)) sapply(highlight_config, function(x) x$color %||% "black") else character(0)
        highlight_names <- if (!is.null(highlight_config)) sapply(highlight_config, function(x) x$group_name %||% "Unknown") else character(0)
        fdr_color_map <- setNames(c(fdr_color_1, fdr_color_2, "black"), levels(plot_data$fdr_category))
        outline_color_map <- setNames(highlight_colors, highlight_names)
        highlight_fill_map <- setNames(highlight_colors, highlight_names)
        p <- p +
            scale_color_manual(
                name = if (sig_use_plot_thresholds) "Gene Group" else "Outline",
                values = c(fdr_color_map, outline_color_map, "None" = NA),
                breaks = c(highlight_names, "None"),
                na.value = NA,
                guide = if (sig_use_plot_thresholds) guide_legend(override.aes = list(size = 4)) else "none"
            ) +
            (
                if (sig_use_plot_thresholds) {
                    scale_fill_manual(
                        name = "Significance",
                        values = c(fdr_color_1, fdr_color_2, "black"),
                        labels = levels(plot_data$fdr_category),
                        guide = guide_legend(override.aes = list(size = 4))
                    )
                } else {
                    scale_fill_manual(
                        name = "Gene Group",
                        values = highlight_fill_map,
                        breaks = if (length(highlight_names) > 0) highlight_names else NULL,
                        guide = guide_legend(override.aes = list(size = 4))
                    )
                }
            )
    } else {
        # When not coloring by FDR, all points colored by highlight group
        highlight_colors <- if (!is.null(highlight_config)) sapply(highlight_config, function(x) x$color %||% "black") else character(0)
        highlight_names <- if (!is.null(highlight_config)) sapply(highlight_config, function(x) x$group_name %||% "Unknown") else character(0)
        
        zero_line <- if (plot_metric == "lfc") geom_hline(yintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7) else NULL
        if (scale_point_size) {
            # Scale size_values by point_size
            scaled_size_values <- size_values * (point_size / 2)
            p <- ggplot(plot_data, aes(x = rank, y = .data[[y_col]])) +
                zero_line +
                geom_point(
                    data = plot_data %>% filter(highlight_group == "none"),
                    aes(size = fdr_bin, text = paste0(.data[[gene_id_col]], "<br>Rank: ", rank, "<br>", metric_label, ": ", round(.data[[y_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))),
                    color = "gray60", alpha = 0.6
                ) +
                geom_point(
                    data = plot_data %>% filter(highlight_group != "none"),
                    aes(color = highlight_group, size = fdr_bin, text = paste0(.data[[gene_id_col]], "<br>Rank: ", rank, "<br>", metric_label, ": ", round(.data[[y_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))),
                    alpha = 0.8
                ) +
                scale_size_manual(name = "FDR (size)", values = scaled_size_values)
        } else {
            p <- ggplot(plot_data, aes(x = rank, y = .data[[y_col]])) +
                geom_point(
                    data = plot_data %>% filter(highlight_group == "none"),
                    aes(text = paste0(.data[[gene_id_col]], "<br>Rank: ", rank, "<br>", metric_label, ": ", round(.data[[y_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))),
                    color = "gray60", size = point_size, alpha = 0.6
                ) +
                geom_point(
                    data = plot_data %>% filter(highlight_group != "none"),
                    aes(color = highlight_group, text = paste0(.data[[gene_id_col]], "<br>Rank: ", rank, "<br>", metric_label, ": ", round(.data[[y_col]], 3), "<br>FDR: ", round(.data[[fdr_col]], 5))),
                    size = point_size, alpha = 0.8
                )
        }
        
        # Map all groups including "none" (gray for non-highlighted)
        all_colors <- c(setNames(highlight_colors, highlight_names), "none" = "gray60")
        p <- p + scale_color_manual(
            name = "Gene Group",
            values = all_colors,
            breaks = if (length(highlight_names) > 0) highlight_names else NULL
        ) + guides(color = guide_legend(override.aes = list(size = 4)))
    }

    
    if (show_labels && !is.null(genes_to_label) && nrow(genes_to_label) > 0) {
        p <- p + geom_text_repel(
            data = genes_to_label,
            aes(label = .data[[gene_id_col]]),
            color = ifelse(!is.na(genes_to_label$highlight_color),
                          genes_to_label$highlight_color, "black"),
            size = 3 * label_scale,
            max.overlaps = Inf,
            box.padding = 0.8,
            point.padding = 0.8,
            min.segment.length = 0,
            force = 3,
            force_pull = 1,
            fontface = "bold"
        )
    }
    
    # Format title
    plot_title <- if (custom_title != "" && !is.null(custom_title)) {
        custom_title
    } else {
        clean_name <- gsub("_", " ", comparison_name)
        paste0(clean_name, " (", title_suffix, ")")
    }
    
    p <- p +
        labs(
            title = plot_title,
            x = if (rank_method == "rra") "Rank (RRA)" else "Rank (LFC)",
            y = y_label
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14 * text_scale),
            axis.title = element_text(size = 11 * text_scale),
            axis.title.x = element_text(size = 11 * text_scale),
            axis.text.y = element_text(size = 10 * text_scale),
            axis.text.x = element_text(size = 10 * text_scale),
            axis.ticks.x = element_line(),
            legend.title = element_text(size = 11 * text_scale),
            legend.text = element_text(size = 10 * text_scale),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            legend.position = "right",
            legend.box = "vertical"
        )
    
    # Reverse x-axis if requested
    if (reverse_rank) {
        p <- p + scale_x_reverse()
    }
    
    if (is.null(legend_data)) {
        p <- p + guides(color = "none")
    }
    
    return(p)
}

# Run the app
shinyApp(ui = ui, server = server)

# CRISPR Screen Explorer

Interactive Shiny app for exploring CRISPR screen results.

## Versions

- `app.R` — active app, currently v1.1
- `crispr-screen-app-v1.1.R` — versioned source for v1.1
- `crispr-screen-app-v1.0.R` — prior version, kept for archival reference

## R Packages

**Required:**
`shiny`, `tidyverse`, `yaml`, `ggrepel`, `colourpicker`, `shinyFiles`

**Optional (enhances functionality if installed):**
- `plotly` — interactive plots
- `DT` — interactive data tables
- `openxlsx2` or `openxlsx` — Excel export
- `ggExtra` — marginal histograms on count scatter plots

Install required packages via `install_crispr-screen.sh` (see `Dockerfile.R.4.4.3.crispr-screen`).

## Setup

Before running, create the `www/` directory and fetch the required JavaScript dependency:

```bash
mkdir -p www
curl -fsSL https://cdn.jsdelivr.net/npm/sortablejs@1.15.2/Sortable.min.js -o www/Sortable.min.js
```

## Running

```bash
Rscript app.R
# or
R -e "shiny::runApp('.')"
```

Set `APP_ROOT_PATH` to point the file browser at your experiments directory:

```bash
APP_ROOT_PATH=/path/to/experiments Rscript app.R
```

## Configuration

Copy `config.example.yaml` to your experiment directory as `config.yaml` and edit:

- `experiment_name` — used to locate output files
- `comparisons` — treatment:control pairs for RRA analysis
- `output_dir` — where pipeline output is written (default: `output/`)

The app loads a `config.yaml` to find all output files; no paths need to be set in the app itself.

## Expected Experiment Directory Structure

```
my_experiment/
└── crispr-screen/
    ├── config.yaml
    ├── sgrnas.txt
    ├── sample_metadata.yaml
    └── output/
        ├── <experiment_name>_sgrna_count_matrix.txt
        ├── <experiment_name>_sgrna_cpm_matrix.txt
        ├── <experiment_name>_barplot.png
        ├── count_distributions/
        │   └── <sample>_fastq_histogram.png  (and variants)
        ├── <sample>_fastq/
        │   └── <sample>_fastq.countsummary.txt
        ├── PCA_variance/
        │   └── pca_*.png
        └── <experiment_name>_rra_<comparison>/
            └── (MAGeCK RRA results)
```

The app auto-discovers output files from `output_dir` — run the pipeline first, then point the app at the `config.yaml`.

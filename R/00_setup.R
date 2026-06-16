# 00_setup.R
# Project setup and configuration

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

source(path("R", "utils.R"))

required_packages <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "readr",
  "stringr",
  "tibble",
  "janitor",
  "ggplot2"
)

missing_packages <- required_packages[
  !purrr::map_lgl(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running the pipeline:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

DIRS <- list(
  raw = path("data", "raw"),
  interim = path("data", "interim"),
  processed = path("data", "processed"),
  metadata = path("data", "metadata"),
  tables = path("outputs", "tables"),
  figures = path("outputs", "figures"),
  models = path("outputs", "models"),
  logs = path("outputs", "logs")
)

purrr::walk(DIRS, ensure_dir)

CONFIG <- list(
  project_name = "ethnic-density-network",

  census_table_codes = c(
    "uk021", "uk004", "uk015", "uk029", "uk030",
    "uk007a", "uk003", "uk017", "uk023", "uk044",
    "uk054", "uk045", "uk062", "uk063", "uk066",
    "uk067", "uk002"
  ),

  oa_lookup_pattern = "OA21_RGN22_LU",

  processed_dataset = path(
    "data",
    "processed",
    "oa_ethnic_density_network_analysis.rds"
  ),

  processed_dataset_csv = path(
    "data",
    "processed",
    "oa_ethnic_density_network_analysis.csv"
  ),

  variable_map_path = path("data", "metadata", "variable_map.csv"),

  variable_map_template_path = path(
    "data",
    "metadata",
    "variable_map_template.csv"
  ),

  ethnic_density_vars = c(
    "share_white",
    "share_asian",
    "share_indian",
    "share_pakistani",
    "share_bangladeshi",
    "share_black",
    "share_mixed",
    "share_other_ethnic"
  ),

  categorical_vars = c(
    "region"
  ),

  model = list(
    min_complete_cases = 500,
    transform_proportions = TRUE,
    mgm_lambda_selection = "EBIC",
    mgm_lambda_gamma = 0.25
  )
)

log_message("Project root: ", PROJECT_ROOT)
log_message("Setup complete.")
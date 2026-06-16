# 02_clean_variables.R
# Clean imported tables and standardise key identifiers

source(file.path("R", "00_setup.R"))

log_message("Starting cleaning step.")

imported_files <- list.files(
  DIRS$interim,
  pattern = "^imported_.*\\.rds$",
  full.names = TRUE
)

if (length(imported_files) == 0) {
  stop("No imported .rds files found. Run R/01_import.R first.")
}

clean_one_table <- function(file) {
  object_name <- tools::file_path_sans_ext(basename(file))
  log_message("Cleaning: ", object_name)

  df <- readRDS(file) |>
    janitor::clean_names()

  df <- df |>
    dplyr::mutate(
      dplyr::across(where(is.character), normalise_chr)
    )

  oa_col <- detect_oa_column(df, required = FALSE)

  if (!is.null(oa_col) && oa_col != "oa21cd") {
    df <- df |>
      dplyr::rename(oa21cd = dplyr::all_of(oa_col))
  }

  if ("oa21cd" %in% names(df)) {
    df <- df |>
      dplyr::mutate(oa21cd = normalise_chr(.data$oa21cd))
  }

  cleaned_path <- path(
    "data",
    "interim",
    stringr::str_replace(basename(file), "^imported_", "cleaned_")
  )

  safe_save_rds(df, cleaned_path)

  tibble::tibble(
    table = object_name,
    n_rows = nrow(df),
    n_cols = ncol(df),
    has_oa21cd = "oa21cd" %in% names(df),
    columns = paste(names(df), collapse = "; ")
  )
}

diagnostics <- purrr::map_dfr(imported_files, clean_one_table)

safe_write_csv(
  diagnostics,
  path("outputs", "tables", "cleaned_table_diagnostics.csv")
)

# Additional cleaning for OA-region lookup if present.
lookup_path <- path("data", "interim", "cleaned_oa_region_lookup.rds")

if (file.exists(lookup_path)) {
  lookup <- readRDS(lookup_path)

  oa_col <- find_column(
    lookup,
    candidates = c("oa21cd", "oa21_code", "oa_code"),
    regex = "oa21",
    required = TRUE
  )

  region_code_col <- find_column(
    lookup,
    candidates = c("rgn22cd", "region_code", "rgn_code"),
    regex = "rgn.*cd|region.*code",
    required = FALSE
  )

  region_name_col <- find_column(
    lookup,
    candidates = c("rgn22nm", "region_name", "rgn_name"),
    regex = "rgn.*nm|region.*name",
    required = FALSE
  )

  lookup2 <- lookup |>
    dplyr::rename(oa21cd = dplyr::all_of(oa_col)) |>
    dplyr::mutate(oa21cd = normalise_chr(.data$oa21cd))

  if (!is.null(region_code_col) && region_code_col != "rgn22cd") {
    lookup2 <- lookup2 |>
      dplyr::rename(rgn22cd = dplyr::all_of(region_code_col))
  }

  if (!is.null(region_name_col) && region_name_col != "region") {
    lookup2 <- lookup2 |>
      dplyr::rename(region = dplyr::all_of(region_name_col))
  }

  lookup2 <- lookup2 |>
    dplyr::select(dplyr::any_of(c("oa21cd", "rgn22cd", "region"))) |>
    dplyr::distinct()

  check_unique_key(lookup2, "oa21cd", "oa_region_lookup")

  safe_save_rds(
    lookup2,
    path("data", "interim", "oa_region_lookup_cleaned.rds")
  )
}

log_message("Cleaning step complete.")
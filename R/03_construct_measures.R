# 03_construct_measures.R
# Construct OA-level ethnic-density and contextual measures

source(file.path("R", "00_setup.R"))

log_message("Starting measure construction step.")

cleaned_files <- list.files(
  DIRS$interim,
  pattern = "^cleaned_.*\\.rds$",
  full.names = TRUE
)

cleaned_files <- cleaned_files[
  !stringr::str_detect(basename(cleaned_files), "oa_region_lookup")
]

if (length(cleaned_files) == 0) {
  stop("No cleaned table files found. Run R/02_clean_variables.R first.")
}

# Create variable-map template if not present.
# The user/Codex should fill this using Census metadata before final construction.
if (!file.exists(CONFIG$variable_map_path)) {
  template <- tibble::tibble(
    table_code = character(),
    source_col = character(),
    source_cols = character(),
    denominator_col = character(),
    analytic_var = character(),
    measure_type = character(),
    notes = character()
  )

  safe_write_csv(template, CONFIG$variable_map_template_path)

  stop(
    "No variable map found at: ", CONFIG$variable_map_path, "\n",
    "A template has been created at: ", CONFIG$variable_map_template_path, "\n",
    "Fill it with source columns and analytic variable names, then save as variable_map.csv."
  )
}

variable_map <- readr::read_csv(
  CONFIG$variable_map_path,
  show_col_types = FALSE
) |>
  janitor::clean_names()

assert_columns(
  variable_map,
  c("table_code", "analytic_var"),
  "variable_map"
)

if (!"source_col" %in% names(variable_map)) {
  variable_map$source_col <- NA_character_
}

if (!"source_cols" %in% names(variable_map)) {
  variable_map$source_cols <- NA_character_
}

if (!"denominator_col" %in% names(variable_map)) {
  variable_map$denominator_col <- NA_character_
}

parse_source_cols <- function(source_col, source_cols = NA_character_) {
  cols <- character()

  if (!is.na(source_cols) && source_cols != "") {
    cols <- stringr::str_split(source_cols, ";", simplify = FALSE)[[1]]
  } else if (!is.na(source_col) && source_col != "") {
    cols <- source_col
  }

  cols |>
    stringr::str_trim() |>
    purrr::map_chr(janitor::make_clean_names) |>
    purrr::discard(~ is.na(.x) || .x == "")
}

clean_mapped_col <- function(x) {
  purrr::map_chr(x, function(value) {
    if (is.na(value) || value == "") {
      return(NA_character_)
    }

    janitor::make_clean_names(value)
  })
}

variable_map <- variable_map |>
  dplyr::mutate(
    table_code = stringr::str_to_lower(.data$table_code),
    source_col = clean_mapped_col(.data$source_col),
    source_cols = purrr::map_chr(
      .data$source_cols,
      ~ paste(parse_source_cols(NA_character_, .x), collapse = ";")
    ),
    denominator_col = clean_mapped_col(.data$denominator_col),
    analytic_var = clean_mapped_col(.data$analytic_var),
    measure_type = .data$measure_type %||% "value"
  )

read_cleaned_table <- function(table_code) {
  f <- path("data", "interim", paste0("cleaned_", table_code, ".rds"))

  if (!file.exists(f)) {
    stop("Could not find cleaned table for table_code: ", table_code)
  }

  readRDS(f)
}

construct_from_table <- function(table_code, map_rows) {
  df <- read_cleaned_table(table_code)

  if (!"oa21cd" %in% names(df)) {
    stop("Table ", table_code, " does not contain oa21cd.")
  }

  mapped_source_cols <- purrr::map2(
    map_rows$source_col,
    map_rows$source_cols,
    parse_source_cols
  )

  needed <- unique(c(unlist(mapped_source_cols), map_rows$denominator_col))
  needed <- needed[!is.na(needed) & needed != ""]

  missing <- setdiff(needed, names(df))

  if (length(missing) > 0) {
    stop(
      "Table ", table_code, " is missing mapped columns: ",
      paste(missing, collapse = ", ")
    )
  }

  out <- df |>
    dplyr::select(oa21cd)

  for (i in seq_len(nrow(map_rows))) {
    src <- map_rows$source_col[[i]]
    srcs <- mapped_source_cols[[i]]
    denom <- map_rows$denominator_col[[i]]
    analytic <- map_rows$analytic_var[[i]]
    measure_type <- map_rows$measure_type[[i]]

    if (length(srcs) == 0) {
      stop("No source column mapped for analytic variable: ", analytic)
    }

    numerator <- if (length(srcs) == 1) {
      df[[srcs]]
    } else {
      rowSums(df[, srcs, drop = FALSE], na.rm = FALSE)
    }

    if (!is.na(denom) && denom != "") {
      out[[analytic]] <- safe_prop(numerator, df[[denom]])
    } else if (!is.na(measure_type) && measure_type == "proportion_already") {
      out[[analytic]] <- as.numeric(numerator)
    } else if (!is.na(measure_type) && measure_type == "count_sum") {
      out[[analytic]] <- numerator
    } else {
      out[[analytic]] <- numerator
    }
  }

  out |>
    dplyr::distinct()
}

constructed_tables <- variable_map |>
  dplyr::group_split(.data$table_code) |>
  purrr::map(function(map_rows) {
    table_code <- unique(map_rows$table_code)
    construct_from_table(table_code, map_rows)
  })

analysis_df <- purrr::reduce(
  constructed_tables,
  function(x, y) {
    left_join_checked(x, y, by = "oa21cd", x_name = "analysis_df", y_name = "constructed_table")
  }
)

# Join region lookup if available.
lookup_path <- path("data", "interim", "oa_region_lookup_cleaned.rds")

if (file.exists(lookup_path)) {
  lookup <- readRDS(lookup_path)

  analysis_df <- left_join_checked(
    analysis_df,
    lookup,
    by = "oa21cd",
    x_name = "analysis_df",
    y_name = "oa_region_lookup"
  )
} else {
  warning("No cleaned OA-region lookup found. Region will not be included.")
}

# Diagnostics.
missingness <- missingness_table(analysis_df)

safe_write_csv(
  missingness,
  path("outputs", "tables", "analysis_dataset_missingness.csv")
)

prop_vars <- names(analysis_df)[
  stringr::str_detect(names(analysis_df), "^share_|_share$|prop|proportion|density")
]

prop_checks <- range_check_prop(analysis_df, prop_vars)

safe_write_csv(
  prop_checks,
  path("outputs", "tables", "proportion_range_checks.csv")
)

safe_save_rds(analysis_df, CONFIG$processed_dataset)
safe_write_csv(analysis_df, CONFIG$processed_dataset_csv)

log_message("Analysis dataset rows: ", nrow(analysis_df))
log_message("Analysis dataset columns: ", ncol(analysis_df))
log_message("Measure construction step complete.")

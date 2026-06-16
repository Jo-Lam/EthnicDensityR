# 01_import.R
# Import raw source files without modifying them

source(file.path("R", "00_setup.R"))

log_message("Starting import step.")

raw_files <- list.files(
  DIRS$raw,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

metadata_files <- list.files(
  DIRS$metadata,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

all_source_files <- c(raw_files, metadata_files)

inventory <- make_file_inventory(all_source_files)

safe_write_csv(
  inventory,
  path("outputs", "tables", "raw_file_inventory.csv")
)

# Extract zip files to data/interim/_unzipped without modifying raw files.
zip_files <- raw_files[tolower(tools::file_ext(raw_files)) == "zip"]

unzipped_dir <- path("data", "interim", "_unzipped")
ensure_dir(unzipped_dir)

if (length(zip_files) > 0) {
  purrr::walk(zip_files, function(z) {
    out_dir <- file.path(
      unzipped_dir,
      tools::file_path_sans_ext(basename(z))
    )

    ensure_dir(out_dir)

    log_message("Extracting zip: ", z)
    utils::unzip(z, exdir = out_dir)
  })
}

candidate_data_files <- c(
  list.files(
    DIRS$raw,
    pattern = "\\.(csv|tsv|txt)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  ),
  list.files(
    unzipped_dir,
    pattern = "\\.(csv|tsv|txt)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
)

candidate_inventory <- make_file_inventory(candidate_data_files)

safe_write_csv(
  candidate_inventory,
  path("outputs", "tables", "candidate_data_file_inventory.csv")
)

find_table_file <- function(table_code) {
  hits <- candidate_data_files[
    stringr::str_detect(
      stringr::str_to_lower(basename(candidate_data_files)),
      stringr::str_to_lower(table_code)
    )
  ]

  if (length(hits) == 0) {
    warning("No candidate file found for table code: ", table_code)
    return(NULL)
  }

  if (length(hits) > 1) {
    log_message(
      "Multiple files found for ", table_code,
      ". Using first after sorting. Check inventory if incorrect."
    )
    hits <- sort(hits)
  }

  hits[[1]]
}

imported_tables <- list()

for (table_code in CONFIG$census_table_codes) {
  f <- find_table_file(table_code)

  if (is.null(f)) {
    next
  }

  df <- read_csv_flexible(f)

  imported_tables[[table_code]] <- df

  safe_save_rds(
    df,
    path("data", "interim", paste0("imported_", table_code, ".rds"))
  )

  safe_write_csv(
    tibble::tibble(
      table_code = table_code,
      source_file = f,
      n_rows = nrow(df),
      n_cols = ncol(df),
      columns = paste(names(df), collapse = "; ")
    ),
    path("outputs", "tables", paste0("import_summary_", table_code, ".csv"))
  )
}

# Import OA-to-region lookup.
lookup_hits <- all_source_files[
  stringr::str_detect(
    basename(all_source_files),
    stringr::regex(CONFIG$oa_lookup_pattern, ignore_case = TRUE)
  )
]

if (length(lookup_hits) == 0) {
  warning("No OA-region lookup found using pattern: ", CONFIG$oa_lookup_pattern)
} else {
  lookup_file <- lookup_hits[[1]]
  oa_lookup <- read_csv_flexible(lookup_file)

  safe_save_rds(
    oa_lookup,
    path("data", "interim", "imported_oa_region_lookup.rds")
  )

  log_message("Imported OA-region lookup: ", lookup_file)
}

log_message("Import step complete.")
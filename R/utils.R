# utils.R
# Shared helper functions for ethnic-density network analysis

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

log_message <- function(..., .sep = "") {
  msg <- paste0("[", timestamp(), "] ", paste(..., collapse = .sep))
  message(msg)
  invisible(msg)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

safe_write_csv <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_csv(x, path)
  log_message("Wrote: ", path)
  invisible(path)
}

safe_save_rds <- function(x, path) {
  ensure_dir(dirname(path))
  saveRDS(x, path)
  log_message("Wrote: ", path)
  invisible(path)
}

safe_read_csv <- function(path) {
  log_message("Reading: ", path)

  readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    guess_max = 100000
  ) |>
    janitor::clean_names()
}

read_csv_flexible <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("csv", "txt")) {
    return(safe_read_csv(path))
  }

  if (ext == "tsv") {
    log_message("Reading TSV: ", path)
    return(
      readr::read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        guess_max = 100000
      ) |>
        janitor::clean_names()
    )
  }

  stop("Unsupported file type: ", path)
}

normalise_chr <- function(x) {
  x |>
    as.character() |>
    stringr::str_squish() |>
    stringr::str_trim()
}

find_column <- function(df, candidates = character(), regex = NULL, required = TRUE) {
  nms <- names(df)

  direct <- intersect(candidates, nms)

  if (length(direct) > 0) {
    return(direct[[1]])
  }

  if (!is.null(regex)) {
    hit <- nms[stringr::str_detect(nms, regex)]
    if (length(hit) > 0) {
      return(hit[[1]])
    }
  }

  if (required) {
    stop(
      "Could not find required column. Candidates: ",
      paste(candidates, collapse = ", "),
      if (!is.null(regex)) paste0("; regex: ", regex) else ""
    )
  }

  NULL
}

detect_oa_column <- function(df, required = TRUE) {
  find_column(
    df,
    candidates = c(
      "oa21cd", "oa", "oa_code", "output_area_code",
      "geography_code", "geography", "date_code"
    ),
    regex = "^oa$|oa21|output_area|geography_code|oa_code",
    required = required
  )
}

standardise_oa_column <- function(df) {
  oa_col <- detect_oa_column(df, required = TRUE)

  df |>
    dplyr::rename(oa21cd = dplyr::all_of(oa_col)) |>
    dplyr::mutate(oa21cd = normalise_chr(.data$oa21cd))
}

assert_columns <- function(df, cols, df_name = deparse(substitute(df))) {
  missing <- setdiff(cols, names(df))

  if (length(missing) > 0) {
    stop(
      df_name,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
}

check_unique_key <- function(df, key, df_name = deparse(substitute(df))) {
  assert_columns(df, key, df_name)

  dupes <- df |>
    dplyr::count(dplyr::across(dplyr::all_of(key)), name = "n") |>
    dplyr::filter(.data$n > 1)

  if (nrow(dupes) > 0) {
    stop(
      df_name,
      " has duplicate key values for: ",
      paste(key, collapse = ", "),
      ". Number of duplicate key rows: ",
      nrow(dupes)
    )
  }

  invisible(TRUE)
}

left_join_checked <- function(x, y, by, x_name = "x", y_name = "y") {
  n_before <- nrow(x)

  y_key <- by
  check_unique_key(y, y_key, y_name)

  out <- dplyr::left_join(x, y, by = by)

  n_after <- nrow(out)

  if (n_after != n_before) {
    stop(
      "Join changed row count: ",
      x_name, " had ", n_before, " rows; output has ", n_after,
      " rows. Check join key uniqueness."
    )
  }

  unmatched <- out |>
    dplyr::filter(dplyr::if_any(
      dplyr::all_of(setdiff(names(y), by)),
      ~ !is.na(.x)
    )) |>
    nrow()

  log_message(
    "Join completed: ", x_name, " + ", y_name,
    "; rows before/after = ", n_before, "/", n_after
  )

  invisible(out)
}

safe_prop <- function(num, denom) {
  dplyr::if_else(
    is.na(denom) | denom <= 0,
    NA_real_,
    as.numeric(num) / as.numeric(denom)
  )
}

empirical_logit <- function(p, eps = 0.0001) {
  p2 <- pmin(pmax(p, eps), 1 - eps)
  log(p2 / (1 - p2))
}

range_check_prop <- function(df, vars) {
  vars <- intersect(vars, names(df))

  out <- purrr::map_dfr(vars, function(v) {
    x <- df[[v]]

    tibble::tibble(
      variable = v,
      n = length(x),
      n_missing = sum(is.na(x)),
      min = suppressWarnings(min(x, na.rm = TRUE)),
      max = suppressWarnings(max(x, na.rm = TRUE)),
      n_below_0 = sum(x < 0, na.rm = TRUE),
      n_above_1 = sum(x > 1, na.rm = TRUE)
    )
  })

  bad <- out |>
    dplyr::filter(.data$n_below_0 > 0 | .data$n_above_1 > 0)

  if (nrow(bad) > 0) {
    warning("Some proportion variables are outside [0, 1]. Check outputs/tables/proportion_range_checks.csv")
  }

  out
}

missingness_table <- function(df) {
  tibble::tibble(
    variable = names(df),
    n = nrow(df),
    n_missing = purrr::map_int(df, ~ sum(is.na(.x))),
    pct_missing = purrr::map_dbl(df, ~ mean(is.na(.x)) * 100),
    class = purrr::map_chr(df, ~ paste(class(.x), collapse = "/"))
  ) |>
    dplyr::arrange(dplyr::desc(.data$pct_missing), .data$variable)
}

numeric_summary_table <- function(df, vars = NULL) {
  if (is.null(vars)) {
    vars <- names(df)[purrr::map_lgl(df, is.numeric)]
  }

  vars <- intersect(vars, names(df))

  purrr::map_dfr(vars, function(v) {
    x <- df[[v]]

    tibble::tibble(
      variable = v,
      n = sum(!is.na(x)),
      n_missing = sum(is.na(x)),
      pct_missing = mean(is.na(x)) * 100,
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      p00 = stats::quantile(x, 0.00, na.rm = TRUE),
      p01 = stats::quantile(x, 0.01, na.rm = TRUE),
      p25 = stats::quantile(x, 0.25, na.rm = TRUE),
      p50 = stats::quantile(x, 0.50, na.rm = TRUE),
      p75 = stats::quantile(x, 0.75, na.rm = TRUE),
      p99 = stats::quantile(x, 0.99, na.rm = TRUE),
      p100 = stats::quantile(x, 1.00, na.rm = TRUE)
    )
  })
}

make_file_inventory <- function(paths) {
  tibble::tibble(path = paths) |>
    dplyr::mutate(
      file_name = basename(.data$path),
      extension = tools::file_ext(.data$path),
      size_mb = file.info(.data$path)$size / 1024^2,
      modified_time = file.info(.data$path)$mtime
    ) |>
    dplyr::arrange(.data$file_name)
}

save_ggplot <- function(plot, path, width = 8, height = 6, dpi = 300) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  log_message("Wrote: ", path)
  invisible(path)
}

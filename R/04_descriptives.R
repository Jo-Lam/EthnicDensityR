# 04_descriptives.R
# Descriptive summaries and exploratory plots

source(file.path("R", "00_setup.R"))

log_message("Starting descriptive analysis step.")

if (!file.exists(CONFIG$processed_dataset)) {
  stop("Processed dataset not found. Run R/03_construct_measures.R first.")
}

df <- readRDS(CONFIG$processed_dataset)

missingness <- missingness_table(df)

safe_write_csv(
  missingness,
  path("outputs", "tables", "descriptive_missingness.csv")
)

numeric_vars <- names(df)[purrr::map_lgl(df, is.numeric)]

numeric_summary <- numeric_summary_table(df, numeric_vars)

safe_write_csv(
  numeric_summary,
  path("outputs", "tables", "numeric_summary.csv")
)

ethnic_vars <- intersect(CONFIG$ethnic_density_vars, names(df))

if (length(ethnic_vars) > 0) {
  ethnic_summary <- numeric_summary_table(df, ethnic_vars)

  safe_write_csv(
    ethnic_summary,
    path("outputs", "tables", "ethnic_density_summary.csv")
  )

  for (v in ethnic_vars) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[v]])) +
      ggplot2::geom_histogram(bins = 50) +
      ggplot2::labs(
        title = paste("Distribution of", v),
        x = v,
        y = "Number of OAs"
      ) +
      ggplot2::theme_minimal()

    save_ggplot(
      p,
      path("outputs", "figures", paste0("hist_", v, ".png")),
      width = 7,
      height = 5
    )
  }
}

# Correlation matrix for numeric variables with manageable missingness.
corr_vars <- numeric_summary |>
  dplyr::filter(.data$pct_missing < 30) |>
  dplyr::pull(.data$variable)

corr_vars <- intersect(corr_vars, names(df))

if (length(corr_vars) >= 2) {
  corr_mat <- stats::cor(
    df[, corr_vars, drop = FALSE],
    use = "pairwise.complete.obs"
  )

  corr_long <- as.data.frame(as.table(corr_mat)) |>
    tibble::as_tibble() |>
    dplyr::rename(
      var1 = Var1,
      var2 = Var2,
      correlation = Freq
    )

  safe_write_csv(
    corr_long,
    path("outputs", "tables", "numeric_correlation_matrix_long.csv")
  )

  p_corr <- ggplot2::ggplot(
    corr_long,
    ggplot2::aes(x = .data$var1, y = .data$var2, fill = .data$correlation)
  ) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      title = "Pairwise correlations among numeric variables",
      x = NULL,
      y = NULL,
      fill = "r"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  save_ggplot(
    p_corr,
    path("outputs", "figures", "numeric_correlation_heatmap.png"),
    width = 10,
    height = 9
  )
}

# Regional summaries if region is present.
if ("region" %in% names(df)) {
  region_summary <- df |>
    dplyr::count(.data$region, name = "n_oa") |>
    dplyr::arrange(dplyr::desc(.data$n_oa))

  safe_write_csv(
    region_summary,
    path("outputs", "tables", "region_oa_counts.csv")
  )
}

log_message("Descriptive analysis step complete.")
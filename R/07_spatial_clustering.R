# 07_spatial_clustering.R
# Spatial clustering diagnostics for OA-level ethnic density and contextual variables

source(file.path("R", "00_setup.R"))

log_message("Starting spatial clustering step.")

spatial_packages <- c("jsonlite", "sf", "spdep")

missing_spatial_packages <- spatial_packages[
  !purrr::map_lgl(spatial_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_spatial_packages) > 0) {
  stop(
    "Install missing spatial packages before running spatial clustering:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_spatial_packages), collapse = ", "),
    "))"
  )
}

invisible(lapply(spatial_packages, library, character.only = TRUE))

if (!file.exists(CONFIG$processed_dataset)) {
  stop("Processed dataset not found. Run R/03_construct_measures.R first.")
}

df <- readRDS(CONFIG$processed_dataset) |>
  dplyr::filter(stringr::str_starts(.data$oa21cd, "E"))

if (nrow(df) == 0) {
  stop("No England OAs found in processed dataset.")
}

lisa_dir <- path("outputs", "figures", "lisa_maps")
ensure_dir(lisa_dir)

candidate_spatial_files <- list.files(
  path("."),
  pattern = "\\.(shp|gpkg|geojson|json|kml|gml|dbf|shx|prj|csv|rds)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

spatial_inventory <- make_file_inventory(candidate_spatial_files) |>
  dplyr::mutate(
    geometry_candidate = stringr::str_detect(
      stringr::str_to_lower(.data$extension),
      "shp|gpkg|geojson|kml|gml"
    ),
    coordinate_candidate = stringr::str_detect(
      stringr::str_to_lower(.data$file_name),
      "centroid|coordinate|coords|longitude|latitude|easting|northing|pwc"
    )
  )

safe_write_csv(
  spatial_inventory,
  path("outputs", "tables", "spatial_source_inventory.csv")
)

local_geometry <- spatial_inventory |>
  dplyr::filter(.data$geometry_candidate)

local_coordinate <- spatial_inventory |>
  dplyr::filter(.data$coordinate_candidate)

if (nrow(local_geometry) > 0) {
  geometry_assessment <- "Local polygon/geometry files were found. Queen contiguity would be preferred for OA polygons because it directly uses shared boundaries. This script currently uses centroid k-nearest neighbours unless a reviewed polygon workflow is added."
} else if (nrow(local_coordinate) > 0) {
  geometry_assessment <- "Local coordinate-like files were found, but no local polygon geometry was found. k-nearest neighbours are preferred for point coordinates."
} else {
  geometry_assessment <- "No local OA polygons, centroids, shapefiles, GeoPackages, GeoJSON, or coordinate columns were found. The script downloads official ONS OA 2021 England/Wales population-weighted centroids and uses k-nearest neighbours."
}

neighbour_definition <- "k-nearest neighbours, k = 8, based on ONS OA 2021 population-weighted centroids in British National Grid (EPSG:27700)."

safe_write_csv(
  tibble::tibble(
    geometry_available_locally = nrow(local_geometry) > 0,
    coordinate_candidate_available_locally = nrow(local_coordinate) > 0,
    recommended_neighbourhood = neighbour_definition,
    rationale = geometry_assessment
  ),
  path("outputs", "tables", "spatial_neighbourhood_decision.csv")
)

centroid_cache <- path("data", "interim", "oa21_ew_population_weighted_centroids_v4.rds")

download_ons_centroids <- function(cache_path) {
  base_url <- paste0(
    "https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services/",
    "OA_December_2021_EW_PWC_V4/FeatureServer/0/query"
  )

  page_size <- 2000
  offset <- 0
  chunks <- list()

  repeat {
    query <- list(
      where = "1=1",
      outFields = "OA21CD",
      returnGeometry = "true",
      outSR = "27700",
      resultRecordCount = page_size,
      resultOffset = offset,
      f = "json"
    )

    url <- paste0(
      base_url,
      "?",
      paste(
        paste0(
          utils::URLencode(names(query), reserved = TRUE),
          "=",
          utils::URLencode(as.character(query), reserved = TRUE)
        ),
        collapse = "&"
      )
    )

    log_message("Downloading centroid records offset: ", offset)
    response <- jsonlite::fromJSON(url, simplifyVector = FALSE)

    if (!is.null(response$error)) {
      stop("ONS centroid download failed: ", response$error$message)
    }

    features <- response$features

    if (length(features) == 0) {
      break
    }

    chunk <- purrr::map_dfr(features, function(feature) {
      tibble::tibble(
        oa21cd = feature$attributes$OA21CD,
        x = feature$geometry$x,
        y = feature$geometry$y
      )
    })

    chunks[[length(chunks) + 1]] <- chunk

    if (is.null(response$exceededTransferLimit) || !isTRUE(response$exceededTransferLimit)) {
      break
    }

    offset <- offset + page_size
  }

  centroids <- dplyr::bind_rows(chunks) |>
    dplyr::distinct(.data$oa21cd, .keep_all = TRUE)

  safe_save_rds(centroids, cache_path)

  centroids
}

centroids <- if (file.exists(centroid_cache)) {
  log_message("Using cached ONS OA centroids: ", centroid_cache)
  readRDS(centroid_cache)
} else {
  download_ons_centroids(centroid_cache)
}

eng_centroids <- centroids |>
  dplyr::filter(stringr::str_starts(.data$oa21cd, "E"))

analysis <- df |>
  dplyr::inner_join(eng_centroids, by = "oa21cd")

if (nrow(analysis) != nrow(df)) {
  stop(
    "Centroid join changed row count: processed England rows = ",
    nrow(df),
    "; joined rows = ",
    nrow(analysis),
    ". Check OA code coverage."
  )
}

coords <- as.matrix(analysis[, c("x", "y")])

k_neighbours <- 8
weights_cache <- path("data", "interim", "spatial_weights_knn8_england.rds")

if (file.exists(weights_cache)) {
  log_message("Using cached spatial weights: ", weights_cache)
  weights <- readRDS(weights_cache)
  nb <- weights$nb
  listw <- weights$listw
} else {
  knn <- spdep::knearneigh(coords, k = k_neighbours, longlat = FALSE)
  nb <- spdep::knn2nb(knn, sym = TRUE)
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  safe_save_rds(
    list(nb = nb, listw = listw, k = k_neighbours, oa21cd = analysis$oa21cd),
    weights_cache
  )
}

neighbour_summary <- tibble::tibble(
  n_areas = nrow(analysis),
  neighbour_definition = neighbour_definition,
  min_neighbours = min(spdep::card(nb)),
  median_neighbours = stats::median(spdep::card(nb)),
  max_neighbours = max(spdep::card(nb)),
  zero_neighbour_areas = sum(spdep::card(nb) == 0)
)

safe_write_csv(
  neighbour_summary,
  path("outputs", "tables", "spatial_neighbour_summary.csv")
)

ethnic_vars <- intersect(CONFIG$ethnic_density_vars, names(analysis))

variable_map <- readr::read_csv(CONFIG$variable_map_path, show_col_types = FALSE) |>
  janitor::clean_names() |>
  dplyr::mutate(analytic_var = janitor::make_clean_names(.data$analytic_var))

contextual_vars <- variable_map |>
  dplyr::pull(.data$analytic_var) |>
  unique() |>
  setdiff(ethnic_vars) |>
  intersect(names(analysis))

analysis_vars <- unique(c(ethnic_vars, contextual_vars))

if (length(analysis_vars) == 0) {
  stop("No analysis variables found for spatial clustering.")
}

zscore <- function(x) {
  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}

classify_lisa <- function(z, lag_z, p_value, alpha = 0.05) {
  dplyr::case_when(
    is.na(z) | is.na(lag_z) | is.na(p_value) ~ "Not evaluated",
    p_value > alpha ~ "Not significant",
    z > 0 & lag_z > 0 ~ "High-High",
    z < 0 & lag_z < 0 ~ "Low-Low",
    z > 0 & lag_z < 0 ~ "High-Low",
    z < 0 & lag_z > 0 ~ "Low-High",
    TRUE ~ "Not significant"
  )
}

plot_lisa_map <- function(lisa_df, variable) {
  cluster_levels <- c(
    "High-High",
    "Low-Low",
    "High-Low",
    "Low-High",
    "Not significant",
    "Not evaluated"
  )

  lisa_df <- lisa_df |>
    dplyr::mutate(cluster = factor(.data$cluster, levels = cluster_levels))

  p <- ggplot2::ggplot(
    lisa_df,
    ggplot2::aes(x = .data$x, y = .data$y, colour = .data$cluster)
  ) +
    ggplot2::geom_point(size = 0.08, alpha = 0.75) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_manual(
      values = c(
        "High-High" = "#b2182b",
        "Low-Low" = "#2166ac",
        "High-Low" = "#ef8a62",
        "Low-High" = "#67a9cf",
        "Not significant" = "grey82",
        "Not evaluated" = "grey50"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = paste("Local Moran cluster map:", variable),
      subtitle = "England OAs; k-nearest neighbours on population-weighted centroids",
      x = NULL,
      y = NULL,
      colour = "LISA cluster"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 11),
      plot.subtitle = ggplot2::element_text(size = 8)
    )

  save_ggplot(
    p,
    file.path(lisa_dir, paste0("lisa_", variable, ".png")),
    width = 8,
    height = 9,
    dpi = 300
  )
}

global_rows <- list()
lisa_summary_rows <- list()
rank_rows <- list()

for (variable in analysis_vars) {
  log_message("Calculating Moran statistics for: ", variable)

  values <- analysis[[variable]]
  valid <- !is.na(values)

  if (sum(valid) < 3 || stats::sd(values[valid], na.rm = TRUE) == 0) {
    warning("Skipping variable with too few valid observations or no variation: ", variable)
    next
  }

  moran <- spdep::moran.test(
    values,
    listw,
    zero.policy = TRUE,
    na.action = stats::na.exclude
  )

  global_rows[[length(global_rows) + 1]] <- tibble::tibble(
    variable = variable,
    n = sum(valid),
    moran_i = unname(moran$estimate[["Moran I statistic"]]),
    expected_i = unname(moran$estimate[["Expectation"]]),
    variance_i = unname(moran$estimate[["Variance"]]),
    statistic = unname(moran$statistic),
    p_value = moran$p.value,
    neighbour_definition = neighbour_definition
  )

  local <- spdep::localmoran(
    values,
    listw,
    zero.policy = TRUE,
    na.action = stats::na.exclude
  )

  z <- zscore(values)
  lag_z <- spdep::lag.listw(listw, z, zero.policy = TRUE, NAOK = TRUE)

  lisa_df <- analysis |>
    dplyr::select("oa21cd", "region", "x", "y") |>
    dplyr::mutate(
      variable = variable,
      value = values,
      z_value = z,
      lag_z = lag_z,
      local_i = local[, "Ii"],
      p_value = local[, "Pr(z != E(Ii))"],
      cluster = classify_lisa(.data$z_value, .data$lag_z, .data$p_value)
    )

  lisa_summary <- lisa_df |>
    dplyr::count(.data$variable, .data$cluster, name = "n_oa") |>
    dplyr::group_by(.data$variable) |>
    dplyr::mutate(pct_oa = .data$n_oa / sum(.data$n_oa) * 100) |>
    dplyr::ungroup()

  lisa_summary_rows[[length(lisa_summary_rows) + 1]] <- lisa_summary

  significant_pct <- lisa_df |>
    dplyr::filter(.data$cluster %in% c("High-High", "Low-Low", "High-Low", "Low-High")) |>
    dplyr::summarise(pct_significant = n() / nrow(lisa_df) * 100) |>
    dplyr::pull(.data$pct_significant)

  rank_rows[[length(rank_rows) + 1]] <- tibble::tibble(
    variable = variable,
    n = sum(valid),
    moran_i = unname(moran$estimate[["Moran I statistic"]]),
    pct_lisa_significant = significant_pct
  )

  safe_write_csv(
    lisa_df,
    path("outputs", "tables", paste0("lisa_", variable, ".csv"))
  )

  plot_lisa_map(lisa_df, variable)
}

moran_global <- dplyr::bind_rows(global_rows) |>
  dplyr::arrange(dplyr::desc(.data$moran_i))

lisa_summary <- dplyr::bind_rows(lisa_summary_rows) |>
  dplyr::arrange(.data$variable, .data$cluster)

spatial_clustering_rank <- dplyr::bind_rows(rank_rows) |>
  dplyr::arrange(dplyr::desc(.data$moran_i))

safe_write_csv(
  moran_global,
  path("outputs", "tables", "moran_global.csv")
)

safe_write_csv(
  lisa_summary,
  path("outputs", "tables", "lisa_summary.csv")
)

safe_write_csv(
  spatial_clustering_rank,
  path("outputs", "tables", "spatial_clustering_rank.csv")
)

top_moran_plot <- moran_global |>
  dplyr::slice_max(order_by = .data$moran_i, n = 20, with_ties = FALSE) |>
  dplyr::mutate(variable = stats::reorder(.data$variable, .data$moran_i))

p_moran <- ggplot2::ggplot(
  top_moran_plot,
  ggplot2::aes(x = .data$variable, y = .data$moran_i)
) +
  ggplot2::geom_col(fill = "#2b6cb0") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Variables ranked by global spatial clustering",
    subtitle = "Global Moran's I, England OAs, k-nearest-neighbour weights",
    x = NULL,
    y = "Global Moran's I"
  ) +
  ggplot2::theme_minimal()

save_ggplot(
  p_moran,
  path("outputs", "figures", "moran_global_rank_top20.png"),
  width = 9,
  height = 7
)

lisa_plot_data <- lisa_summary |>
  dplyr::filter(.data$cluster != "Not evaluated") |>
  dplyr::mutate(
    variable = factor(
      .data$variable,
      levels = rev(spatial_clustering_rank$variable)
    )
  )

p_lisa <- ggplot2::ggplot(
  lisa_plot_data,
  ggplot2::aes(x = .data$cluster, y = .data$variable, fill = .data$pct_oa)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(option = "magma", direction = -1) +
  ggplot2::labs(
    title = "Local spatial cluster composition by variable",
    subtitle = "Percentage of England OAs in each LISA class",
    x = NULL,
    y = NULL,
    fill = "% OAs"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

save_ggplot(
  p_lisa,
  path("outputs", "figures", "lisa_cluster_composition_heatmap.png"),
  width = 9,
  height = 10
)

ethnic_lisa_labels <- c(
  share_white = "White",
  share_black = "Black",
  share_asian = "Asian",
  share_bangladeshi = "Bangladeshi",
  share_indian = "Indian",
  share_pakistani = "Pakistani",
  share_mixed = "Mixed",
  share_other_ethnic = "Other ethnic"
)

ethnic_lisa_order <- names(ethnic_lisa_labels)

ethnic_lisa_files <- path(
  "outputs",
  "tables",
  paste0("lisa_", ethnic_lisa_order, ".csv")
)

if (all(file.exists(ethnic_lisa_files))) {
  ethnic_lisa_panel_data <- purrr::map_dfr(
    ethnic_lisa_files,
    ~ readr::read_csv(
      .x,
      col_select = c("oa21cd", "x", "y", "variable", "cluster"),
      show_col_types = FALSE
    )
  ) |>
    dplyr::mutate(
      variable = factor(.data$variable, levels = ethnic_lisa_order),
      cluster = factor(
        .data$cluster,
        levels = c(
          "High-High",
          "Low-Low",
          "High-Low",
          "Low-High",
          "Not significant",
          "Not evaluated"
        )
      )
    )

  p_ethnic_lisa_panel <- ggplot2::ggplot(
    ethnic_lisa_panel_data,
    ggplot2::aes(x = .data$x, y = .data$y, colour = .data$cluster)
  ) +
    ggplot2::geom_point(size = 0.015, alpha = 0.85) +
    ggplot2::coord_equal() +
    ggplot2::facet_wrap(
      ~ variable,
      labeller = ggplot2::as_labeller(ethnic_lisa_labels),
      ncol = 4
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "High-High" = "#b2182b",
        "Low-Low" = "#2166ac",
        "High-Low" = "#ef8a62",
        "Low-High" = "#67a9cf",
        "Not significant" = "grey86",
        "Not evaluated" = "grey55"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Local Moran's I cluster classes for ethnic-density measures",
      subtitle = "England OAs; k-nearest neighbours on population-weighted centroids",
      x = NULL,
      y = NULL,
      colour = "LISA class"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1.4, "lines"),
      legend.key.height = ggplot2::unit(1.4, "lines"),
      panel.spacing = ggplot2::unit(0.6, "lines"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 3.5, alpha = 1),
        nrow = 2,
        byrow = TRUE
      )
    )

  save_ggplot(
    p_ethnic_lisa_panel,
    path("outputs", "figures", "figure_lisa_all_ethnic_density_panel.png"),
    width = 13,
    height = 8,
    dpi = 300
  )

  save_ggplot(
    p_ethnic_lisa_panel,
    path("outputs", "figures", "figure_lisa_all_ethnic_density_panel.pdf"),
    width = 13,
    height = 8
  )
} else {
  warning(
    "Skipping all-ethnic-density LISA panel because one or more LISA CSV files are missing."
  )
}

top_global <- utils::head(spatial_clustering_rank, 10)

report_lines <- c(
  "# Spatial clustering interpretation",
  "",
  paste0("Analysis date: ", Sys.Date()),
  "",
  "## Spatial data and neighbourhood definition",
  "",
  geometry_assessment,
  "",
  paste0("The implemented neighbourhood definition is: ", neighbour_definition),
  "",
  "Queen or rook contiguity would require reviewed OA polygon geometry. Because no local polygon geometry was present in the repository, k-nearest neighbours on official ONS population-weighted centroids provides a reproducible local-neighbour definition and avoids isolated areas.",
  "",
  "## Interpretation",
  "",
  "Global Moran's I summarises whether similar values tend to occur near one another across England OAs. Local Moran's I identifies local clusters and spatial outliers. These results should be interpreted as evidence of spatial concentration and contextual co-location, not as causal effects or pathways.",
  "",
  "The variables with the strongest global spatial clustering were:",
  "",
  paste0(
    seq_len(nrow(top_global)),
    ". ",
    top_global$variable,
    ": Moran's I = ",
    round(top_global$moran_i, 3),
    "; significant LISA share = ",
    round(top_global$pct_lisa_significant, 1),
    "%."
  ),
  "",
  "High-High clusters indicate OAs with high values surrounded by nearby OAs with high values. Low-Low clusters indicate spatial concentration of low values. High-Low and Low-High areas are local spatial outliers. These patterns may reflect residential concentration, regional settlement histories, housing-market sorting, institutional geography, and other contextual processes, but the analysis does not identify causal mechanisms.",
  "",
  "## Outputs",
  "",
  "- outputs/tables/moran_global.csv",
  "- outputs/tables/lisa_summary.csv",
  "- outputs/tables/spatial_clustering_rank.csv",
  "- outputs/figures/lisa_maps/"
)

writeLines(
  report_lines,
  path("outputs", "tables", "spatial_clustering_interpretation.md")
)

log_message("Wrote: ", path("outputs", "tables", "spatial_clustering_interpretation.md"))
log_message("Spatial clustering step complete.")

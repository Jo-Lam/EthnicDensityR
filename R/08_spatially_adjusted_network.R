# 08_spatially_adjusted_network.R
# Compare original England MGM networks with spatially adjusted residual networks

source(file.path("R", "00_setup.R"))

log_message("Starting spatially adjusted network step.")

adjustment_packages <- c("mgm", "qgraph", "spdep")

missing_adjustment_packages <- adjustment_packages[
  !purrr::map_lgl(adjustment_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_adjustment_packages) > 0) {
  stop(
    "Install missing packages before running spatially adjusted networks:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_adjustment_packages), collapse = ", "),
    "))"
  )
}

invisible(lapply(adjustment_packages, library, character.only = TRUE))

if (!file.exists(CONFIG$processed_dataset)) {
  stop("Processed dataset not found. Run R/03_construct_measures.R first.")
}

centroid_cache <- path("data", "interim", "oa21_ew_population_weighted_centroids_v4.rds")

if (!file.exists(centroid_cache)) {
  stop(
    "Centroid cache not found: ",
    centroid_cache,
    ". Run R/07_spatial_clustering.R first."
  )
}

df <- readRDS(CONFIG$processed_dataset) |>
  dplyr::filter(stringr::str_starts(.data$oa21cd, "E"))

centroids <- readRDS(centroid_cache) |>
  dplyr::filter(stringr::str_starts(.data$oa21cd, "E"))

analysis <- df |>
  dplyr::inner_join(centroids, by = "oa21cd")

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

  if (!identical(weights$oa21cd, analysis$oa21cd)) {
    stop("Cached spatial weights OA order does not match analysis data.")
  }

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

if (!file.exists(CONFIG$variable_map_path)) {
  stop("Variable map not found: ", CONFIG$variable_map_path)
}

variable_map <- readr::read_csv(
  CONFIG$variable_map_path,
  show_col_types = FALSE
) |>
  janitor::clean_names() |>
  dplyr::mutate(
    analytic_var = janitor::make_clean_names(.data$analytic_var)
  )

all_mapped_vars <- intersect(variable_map$analytic_var, names(analysis))
ethnic_density_vars <- intersect(CONFIG$ethnic_density_vars, names(analysis))
contextual_vars <- setdiff(all_mapped_vars, ethnic_density_vars)
contextual_vars <- intersect(contextual_vars, names(analysis))

if (length(ethnic_density_vars) == 0 || length(contextual_vars) == 0) {
  stop("Could not identify ethnic density and contextual variables.")
}

base_vars <- unique(c(ethnic_density_vars, contextual_vars))

transform_for_residual <- function(x) {
  x_num <- as.numeric(x)
  looks_like_prop <- all(x_num >= 0 & x_num <= 1, na.rm = TRUE)

  if (looks_like_prop) {
    return(empirical_logit(x_num))
  }

  x_num
}

standardise <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / s
}

residualise_variable <- function(variable) {
  transformed <- transform_for_residual(analysis[[variable]])
  spatial_lag <- spdep::lag.listw(listw, transformed, zero.policy = TRUE, NAOK = TRUE)

  model_df <- tibble::tibble(
    y = transformed,
    spatial_lag = spatial_lag,
    region = factor(analysis$region)
  )

  fit <- stats::lm(y ~ spatial_lag + region, data = model_df, na.action = stats::na.exclude)
  standardise(stats::residuals(fit))
}

log_message("Residualising variables against region fixed effects and spatial lag.")

residual_df <- tibble::tibble(oa21cd = analysis$oa21cd)

for (variable in base_vars) {
  residual_df[[variable]] <- residualise_variable(variable)
}

safe_save_rds(
  residual_df,
  path("data", "processed", "oa_ethnic_density_network_spatial_residuals.rds")
)

safe_write_csv(
  residual_df,
  path("data", "processed", "oa_ethnic_density_network_spatial_residuals.csv")
)

safe_write_csv(
  missingness_table(residual_df),
  path("outputs", "tables", "spatial_residual_missingness.csv")
)

prepare_adjusted_mgm_data <- function(target_var) {
  vars <- unique(c(target_var, contextual_vars))
  vars <- intersect(vars, names(residual_df))

  model_df <- residual_df |>
    dplyr::select(dplyr::all_of(vars)) |>
    tidyr::drop_na()

  if (nrow(model_df) < CONFIG$model$min_complete_cases) {
    stop(
      "Too few complete cases for adjusted network model: ",
      nrow(model_df),
      ". Minimum required: ",
      CONFIG$model$min_complete_cases
    )
  }

  list(
    data = as.matrix(model_df),
    type = rep("g", ncol(model_df)),
    level = rep(1L, ncol(model_df)),
    var_names = names(model_df),
    n = nrow(model_df)
  )
}

extract_edge_table <- function(model, var_names, target_var) {
  wadj <- model$pairwise$wadj

  if (is.null(wadj)) {
    stop("Could not find weighted adjacency matrix in mgm model object.")
  }

  dimnames(wadj) <- list(var_names, var_names)

  as.data.frame(as.table(wadj)) |>
    tibble::as_tibble() |>
    dplyr::rename(
      node1 = Var1,
      node2 = Var2,
      weight = Freq
    ) |>
    dplyr::filter(
      as.character(.data$node1) < as.character(.data$node2),
      !is.na(.data$weight),
      .data$weight != 0
    ) |>
    dplyr::mutate(
      target = target_var,
      involves_target = .data$node1 == target_var | .data$node2 == target_var
    ) |>
    dplyr::arrange(dplyr::desc(abs(.data$weight)))
}

fit_adjusted_network <- function(target_var) {
  model_path <- path(
    "outputs",
    "models",
    paste0("mgm_spatial_adjusted_", target_var, ".rds")
  )

  edge_path <- path(
    "outputs",
    "tables",
    paste0("edges_spatial_adjusted_", target_var, ".csv")
  )

  if (file.exists(model_path) && file.exists(edge_path)) {
    log_message("Using existing spatially adjusted model output for target: ", target_var)
    saved <- readRDS(model_path)
    edge_table <- readr::read_csv(edge_path, show_col_types = FALSE)

    return(
      tibble::tibble(
        target = target_var,
        n_complete = saved$n_complete,
        n_nodes = length(saved$var_names),
        n_edges = nrow(edge_table),
        model_path = model_path
      )
    )
  }

  log_message("Fitting spatially adjusted residual network for target: ", target_var)

  prepared <- prepare_adjusted_mgm_data(target_var)

  model <- mgm::mgm(
    data = prepared$data,
    type = prepared$type,
    level = prepared$level,
    k = 2,
    lambdaSel = CONFIG$model$mgm_lambda_selection,
    lambdaGam = CONFIG$model$mgm_lambda_gamma,
    ruleReg = "OR",
    scale = TRUE,
    pbar = FALSE
  )

  safe_save_rds(
    list(
      target = target_var,
      model = model,
      var_names = prepared$var_names,
      type = prepared$type,
      level = prepared$level,
      n_complete = prepared$n,
      adjustment = "Residualised each variable against region fixed effects and its kNN spatial lag."
    ),
    model_path
  )

  edge_table <- extract_edge_table(
    model = model,
    var_names = prepared$var_names,
    target_var = target_var
  )

  safe_write_csv(edge_table, edge_path)

  png_path <- path(
    "outputs",
    "figures",
    paste0("network_spatial_adjusted_", target_var, ".png")
  )

  grDevices::png(png_path, width = 2400, height = 2000, res = 300)
  qgraph::qgraph(
    model$pairwise$wadj,
    labels = prepared$var_names,
    layout = "spring",
    title = paste("Spatially adjusted residual network:", target_var)
  )
  grDevices::dev.off()

  log_message("Wrote: ", png_path)

  tibble::tibble(
    target = target_var,
    n_complete = prepared$n,
    n_nodes = length(prepared$var_names),
    n_edges = nrow(edge_table),
    model_path = model_path
  )
}

adjusted_summary <- purrr::map_dfr(ethnic_density_vars, fit_adjusted_network)

safe_write_csv(
  adjusted_summary,
  path("outputs", "tables", "spatial_adjusted_network_model_summary.csv")
)

original_edges <- purrr::map_dfr(ethnic_density_vars, function(target_var) {
  edge_path <- path("outputs", "tables", paste0("edges_england_", target_var, ".csv"))

  if (!file.exists(edge_path)) {
    stop(
      "England-only original edge table not found: ",
      edge_path,
      ". Run R/05b_network_models_england.R first."
    )
  }

  readr::read_csv(edge_path, show_col_types = FALSE) |>
    dplyr::mutate(target = target_var)
}) |>
  dplyr::filter(.data$node1 != "region", .data$node2 != "region") |>
  dplyr::mutate(
    edge_id = purrr::map2_chr(
      as.character(.data$node1),
      as.character(.data$node2),
      ~ paste(sort(c(.x, .y)), collapse = " -- ")
    )
  )

adjusted_edges <- purrr::map_dfr(ethnic_density_vars, function(target_var) {
  edge_path <- path("outputs", "tables", paste0("edges_spatial_adjusted_", target_var, ".csv"))

  readr::read_csv(edge_path, show_col_types = FALSE) |>
    dplyr::mutate(target = target_var)
}) |>
  dplyr::mutate(
    edge_id = purrr::map2_chr(
      as.character(.data$node1),
      as.character(.data$node2),
      ~ paste(sort(c(.x, .y)), collapse = " -- ")
    )
  )

edge_persistence <- original_edges |>
  dplyr::select(
    "target",
    "edge_id",
    original_node1 = "node1",
    original_node2 = "node2",
    original_weight = "weight",
    original_involves_target = "involves_target"
  ) |>
  dplyr::full_join(
    adjusted_edges |>
      dplyr::select(
        "target",
        "edge_id",
        adjusted_node1 = "node1",
        adjusted_node2 = "node2",
        adjusted_weight = "weight",
        adjusted_involves_target = "involves_target"
      ),
    by = c("target", "edge_id")
  ) |>
  dplyr::mutate(
    original_present = !is.na(.data$original_weight),
    adjusted_present = !is.na(.data$adjusted_weight),
    persisted = .data$original_present & .data$adjusted_present,
    dropped_after_adjustment = .data$original_present & !.data$adjusted_present,
    appeared_after_adjustment = !.data$original_present & .data$adjusted_present,
    original_abs_weight = abs(.data$original_weight),
    adjusted_abs_weight = abs(.data$adjusted_weight),
    abs_weight_change = .data$adjusted_abs_weight - .data$original_abs_weight,
    involves_target = dplyr::coalesce(
      .data$original_involves_target,
      .data$adjusted_involves_target,
      FALSE
    )
  ) |>
  dplyr::arrange(.data$target, dplyr::desc(.data$persisted), dplyr::desc(.data$original_abs_weight))

safe_write_csv(
  edge_persistence,
  path("outputs", "tables", "spatial_adjusted_edge_persistence.csv")
)

target_edge_persistence <- edge_persistence |>
  dplyr::filter(.data$involves_target)

safe_write_csv(
  target_edge_persistence,
  path("outputs", "tables", "spatial_adjusted_target_edge_persistence.csv")
)

persistence_summary <- edge_persistence |>
  dplyr::group_by(.data$target) |>
  dplyr::summarise(
    original_edges = sum(.data$original_present),
    adjusted_edges = sum(.data$adjusted_present),
    persisted_edges = sum(.data$persisted),
    dropped_edges = sum(.data$dropped_after_adjustment),
    appeared_edges = sum(.data$appeared_after_adjustment),
    persistence_pct = .data$persisted_edges / .data$original_edges * 100,
    .groups = "drop"
  )

target_persistence_summary <- target_edge_persistence |>
  dplyr::group_by(.data$target) |>
  dplyr::summarise(
    original_target_edges = sum(.data$original_present),
    adjusted_target_edges = sum(.data$adjusted_present),
    persisted_target_edges = sum(.data$persisted),
    dropped_target_edges = sum(.data$dropped_after_adjustment),
    appeared_target_edges = sum(.data$appeared_after_adjustment),
    target_persistence_pct = .data$persisted_target_edges / .data$original_target_edges * 100,
    .groups = "drop"
  )

safe_write_csv(
  persistence_summary,
  path("outputs", "tables", "spatial_adjusted_edge_persistence_summary.csv")
)

safe_write_csv(
  target_persistence_summary,
  path("outputs", "tables", "spatial_adjusted_target_edge_persistence_summary.csv")
)

p_persist <- ggplot2::ggplot(
  target_persistence_summary,
  ggplot2::aes(x = .data$target, y = .data$target_persistence_pct)
) +
  ggplot2::geom_col(fill = "#2f855a") +
  ggplot2::coord_cartesian(ylim = c(0, 100)) +
  ggplot2::labs(
    title = "Persistence of target-node edges after spatial adjustment",
    subtitle = "Adjustment: region fixed effects plus kNN spatial lag residualisation",
    x = NULL,
    y = "% original target edges persisting"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

save_ggplot(
  p_persist,
  path("outputs", "figures", "spatial_adjusted_target_edge_persistence.png"),
  width = 9,
  height = 6
)

comparison_targets <- c("share_white", "share_black")

white_black_comparison <- target_edge_persistence |>
  dplyr::filter(.data$target %in% comparison_targets, .data$original_present) |>
  dplyr::mutate(
    neighbour = dplyr::case_when(
      !is.na(.data$original_node1) & .data$original_node1 == .data$target ~ as.character(.data$original_node2),
      !is.na(.data$original_node2) & .data$original_node2 == .data$target ~ as.character(.data$original_node1),
      !is.na(.data$adjusted_node1) & .data$adjusted_node1 == .data$target ~ as.character(.data$adjusted_node2),
      TRUE ~ as.character(.data$adjusted_node1)
    )
  ) |>
  dplyr::group_by(.data$target) |>
  dplyr::slice_max(order_by = .data$original_abs_weight, n = 8, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(
    "target",
    "neighbour",
    original = "original_abs_weight",
    spatially_adjusted = "adjusted_abs_weight"
  ) |>
  tidyr::pivot_longer(
    cols = c("original", "spatially_adjusted"),
    names_to = "model",
    values_to = "abs_weight"
  ) |>
  dplyr::mutate(
    model = dplyr::recode(
      .data$model,
      original = "England baseline",
      spatially_adjusted = "Spatially adjusted"
    ),
    abs_weight = dplyr::coalesce(.data$abs_weight, 0)
  )

p_white_black <- ggplot2::ggplot(
  white_black_comparison,
  ggplot2::aes(
    x = stats::reorder(.data$neighbour, .data$abs_weight),
    y = .data$abs_weight,
    fill = .data$model
  )
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ target, scales = "free_y") +
  ggplot2::labs(
    title = "White and Black ethnic-density target edges before and after spatial adjustment",
    subtitle = "England-only comparison; zero indicates that the target edge was not retained after adjustment.",
    x = NULL,
    y = "Absolute regularised edge weight",
    fill = NULL
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )

save_ggplot(
  p_white_black,
  path("outputs", "figures", "figure_white_black_target_edge_adjustment_comparison.png"),
  width = 10,
  height = 7
)

node_labels <- c(
  share_white = "White",
  share_asian = "Asian",
  share_indian = "Indian",
  share_pakistani = "Pakistani",
  share_bangladeshi = "Bangladeshi",
  share_black = "Black",
  share_mixed = "Mixed",
  share_other_ethnic = "Other ethnic",
  share_born_africa = "Africa-born",
  share_born_middle_east_asia = "Middle East/Asia-born",
  share_born_uk = "UK-born",
  share_children_0_19 = "Children aged 0-19",
  share_christian = "Christian",
  share_hindu = "Hindu",
  share_large_households_5_plus = "Large households",
  share_level4_plus = "Degree-level qualification",
  share_limited_english = "Limited English",
  share_multi_ethnic_households = "Multi-ethnic households",
  share_muslim = "Muslim",
  share_no_religion = "No religion",
  share_one_person_households = "One-person households",
  share_older_65_plus = "Aged 65+",
  share_other_religion = "Other religion",
  share_owned = "Owner occupation",
  share_young_adults_20_34 = "Young adults aged 20-34"
)

all_group_adjustment_comparison <- target_edge_persistence |>
  dplyr::filter(.data$original_present) |>
  dplyr::mutate(
    neighbour = dplyr::case_when(
      !is.na(.data$original_node1) & .data$original_node1 == .data$target ~ as.character(.data$original_node2),
      !is.na(.data$original_node2) & .data$original_node2 == .data$target ~ as.character(.data$original_node1),
      !is.na(.data$adjusted_node1) & .data$adjusted_node1 == .data$target ~ as.character(.data$adjusted_node2),
      TRUE ~ as.character(.data$adjusted_node1)
    ),
    target_label = dplyr::recode(.data$target, !!!node_labels, .default = .data$target),
    neighbour_label = dplyr::recode(.data$neighbour, !!!node_labels, .default = .data$neighbour)
  ) |>
  dplyr::group_by(.data$target) |>
  dplyr::slice_max(order_by = .data$original_abs_weight, n = 5, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(
    "target",
    "target_label",
    "neighbour",
    "neighbour_label",
    original = "original_abs_weight",
    spatially_adjusted = "adjusted_abs_weight"
  ) |>
  tidyr::pivot_longer(
    cols = c("original", "spatially_adjusted"),
    names_to = "model",
    values_to = "abs_weight"
  ) |>
  dplyr::mutate(
    model = dplyr::recode(
      .data$model,
      original = "Original England MGM",
      spatially_adjusted = "Spatially adjusted MGM"
    ),
    abs_weight = dplyr::coalesce(.data$abs_weight, 0),
    target_label = factor(
      .data$target_label,
      levels = c(
        "White",
        "Asian",
        "Indian",
        "Pakistani",
        "Bangladeshi",
        "Black",
        "Mixed",
        "Other ethnic"
      )
    )
  )

safe_write_csv(
  all_group_adjustment_comparison,
  path("outputs", "tables", "spatial_adjusted_all_group_top5_edge_comparison.csv")
)

p_all_groups_adjusted <- ggplot2::ggplot(
  all_group_adjustment_comparison,
  ggplot2::aes(
    x = stats::reorder(.data$neighbour_label, .data$abs_weight),
    y = .data$abs_weight,
    fill = .data$model
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.75),
    width = 0.7
  ) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ target_label, scales = "free_y", ncol = 4) +
  ggplot2::scale_fill_manual(
    values = c(
      "Original England MGM" = "#2b6cb0",
      "Spatially adjusted MGM" = "#dd6b20"
    )
  ) +
  ggplot2::labs(
    title = "Top target-neighbour edges before and after spatial adjustment",
    subtitle = "England-only comparison; zero indicates that the original target edge was not retained after adjustment.",
    x = NULL,
    y = "Absolute regularised edge weight",
    fill = NULL
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8)
  )

save_ggplot(
  p_all_groups_adjusted,
  path("outputs", "figures", "figure_all_groups_target_edge_adjustment_comparison.png"),
  width = 14,
  height = 9
)

save_ggplot(
  p_all_groups_adjusted,
  path("outputs", "figures", "figure_all_groups_target_edge_adjustment_comparison.pdf"),
  width = 14,
  height = 9
)

top_persistent <- target_edge_persistence |>
  dplyr::filter(.data$persisted) |>
  dplyr::group_by(.data$target) |>
  dplyr::slice_max(order_by = .data$adjusted_abs_weight, n = 5, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    neighbour = dplyr::if_else(
      .data$original_node1 == .data$target,
      .data$original_node2,
      .data$original_node1
    )
  ) |>
  dplyr::select("target", "neighbour", "original_weight", "adjusted_weight")

safe_write_csv(
  top_persistent,
  path("outputs", "tables", "spatial_adjusted_top_persistent_target_neighbours.csv")
)

report_lines <- c(
  "# Spatially adjusted network interpretation",
  "",
  paste0("Analysis date: ", Sys.Date()),
  "",
  "## Adjustment used",
  "",
  "Each bounded proportion was empirical-logit transformed. For each variable, a linear model removed English region fixed effects and the variable's own spatial lag based on k-nearest-neighbour centroid weights (k = 8). MGM networks were then estimated on standardised residuals and compared with the England-only baseline MGM.",
  "",
  "This is a spatially adjusted residual-network analysis. It asks whether conditional associations remain after removing broad regional structure and local spatial dependence in each variable. It does not identify causal pathways.",
  "",
  "## Interpretation guidance",
  "",
  "Edges that persist after adjustment are less likely to be solely artefacts of broad region or local spatial clustering. Edges that drop after adjustment may have been strongly spatially patterned or regionally embedded. Edges appearing only after adjustment should be treated cautiously, because residualisation changes the scale and covariance structure.",
  "",
  "Target-edge persistence summary:",
  "",
  paste0(
    "- ",
    target_persistence_summary$target,
    ": ",
    target_persistence_summary$persisted_target_edges,
    " of ",
    target_persistence_summary$original_target_edges,
    " original target edges persisted (",
    round(target_persistence_summary$target_persistence_pct, 1),
    "%)."
  ),
  "",
  "Use the persistence table to focus interpretation on robust ecological conditional associations, while keeping the original networks as the descriptive baseline. Do not interpret dropped edges as disproven relationships; they indicate sensitivity to spatial adjustment.",
  "",
  "## Outputs",
  "",
  "- outputs/tables/spatial_adjusted_edge_persistence.csv",
  "- outputs/tables/spatial_adjusted_target_edge_persistence.csv",
  "- outputs/tables/spatial_adjusted_edge_persistence_summary.csv",
  "- outputs/tables/spatial_adjusted_target_edge_persistence_summary.csv",
  "- outputs/tables/spatial_adjusted_top_persistent_target_neighbours.csv",
  "- outputs/figures/spatial_adjusted_target_edge_persistence.png",
  "- outputs/figures/figure_white_black_target_edge_adjustment_comparison.png"
)

writeLines(
  report_lines,
  path("outputs", "tables", "spatial_adjusted_network_interpretation.md")
)

log_message("Wrote: ", path("outputs", "tables", "spatial_adjusted_network_interpretation.md"))
log_message("Spatially adjusted network step complete.")

# 05_network_models.R
# Estimate UK-wide group-specific mixed graphical/network models

source(file.path("R", "00_setup.R"))

log_message("Starting UK-wide network modelling step.")

model_packages <- c("mgm", "qgraph", "igraph")

missing_model_packages <- model_packages[
  !purrr::map_lgl(model_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_model_packages) > 0) {
  stop(
    "Install missing model packages before running network models:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_model_packages), collapse = ", "),
    "))"
  )
}

if (!file.exists(CONFIG$processed_dataset)) {
  stop("Processed dataset not found. Run R/03_construct_measures.R first.")
}

df <- readRDS(CONFIG$processed_dataset)

log_message("Fitting UK-wide network models: ", nrow(df), " OAs.")

if (!file.exists(CONFIG$variable_map_path)) {
  stop("Variable map not found. Run R/03_construct_measures.R after creating variable_map.csv.")
}

variable_map <- readr::read_csv(
  CONFIG$variable_map_path,
  show_col_types = FALSE
) |>
  janitor::clean_names() |>
  dplyr::mutate(
    analytic_var = janitor::make_clean_names(.data$analytic_var)
  )

all_mapped_vars <- intersect(variable_map$analytic_var, names(df))

ethnic_density_vars <- intersect(CONFIG$ethnic_density_vars, names(df))

if (length(ethnic_density_vars) == 0) {
  stop(
    "No ethnic density variables found in processed data. Expected one or more of: ",
    paste(CONFIG$ethnic_density_vars, collapse = ", ")
  )
}

contextual_vars <- setdiff(all_mapped_vars, ethnic_density_vars)
contextual_vars <- intersect(contextual_vars, names(df))

# English region is not a harmonised UK-wide contextual node. Keep the primary
# UK network to harmonised continuous/contextual measures and use the separate
# England-only network for regional and spatial sensitivity analyses.
contextual_vars <- setdiff(contextual_vars, c("region", "rgn22cd"))
categorical_vars <- character()

prepare_mgm_data <- function(df, vars, categorical_vars = character()) {
  model_df <- df |>
    dplyr::select(dplyr::all_of(vars)) |>
    dplyr::mutate(
      dplyr::across(where(is.character), as.factor)
    )

  # Convert declared categorical variables to factor.
  for (v in intersect(categorical_vars, names(model_df))) {
    model_df[[v]] <- as.factor(model_df[[v]])
  }

  # Drop variables with no variation.
  n_unique <- purrr::map_int(model_df, ~ dplyr::n_distinct(.x, na.rm = TRUE))
  keep <- names(n_unique[n_unique > 1])
  model_df <- model_df[, keep, drop = FALSE]

  model_df <- model_df |>
    tidyr::drop_na()

  if (nrow(model_df) < CONFIG$model$min_complete_cases) {
    stop(
      "Too few complete cases for network model: ",
      nrow(model_df),
      ". Minimum required: ",
      CONFIG$model$min_complete_cases
    )
  }

  types <- character(ncol(model_df))
  levels <- integer(ncol(model_df))

  for (j in seq_along(model_df)) {
    v <- names(model_df)[[j]]
    x <- model_df[[j]]

    if (is.factor(x) || v %in% categorical_vars) {
      f <- as.factor(x)
      model_df[[j]] <- as.integer(f)
      types[[j]] <- "c"
      levels[[j]] <- nlevels(f)
    } else {
      x_num <- as.numeric(x)

      # Transform bounded proportions where appropriate.
      looks_like_prop <- all(x_num >= 0 & x_num <= 1, na.rm = TRUE)

      if (CONFIG$model$transform_proportions && looks_like_prop) {
        x_num <- empirical_logit(x_num)
      }

      model_df[[j]] <- x_num
      types[[j]] <- "g"
      levels[[j]] <- 1L
    }
  }

  list(
    data = as.matrix(model_df),
    type = types,
    level = levels,
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

  edges <- as.data.frame(as.table(wadj)) |>
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

  edges
}

fit_one_network <- function(target_var) {
  log_message("Fitting network for target: ", target_var)

  model_path <- path(
    "outputs",
    "models",
    paste0("mgm_uk_", target_var, ".rds")
  )

  model_path_main <- path(
    "outputs",
    "models",
    paste0("mgm_", target_var, ".rds")
  )

  edge_path <- path(
    "outputs",
    "tables",
    paste0("edges_uk_", target_var, ".csv")
  )

  edge_path_main <- path(
    "outputs",
    "tables",
    paste0("edges_", target_var, ".csv")
  )

  if (file.exists(model_path) && file.exists(edge_path)) {
    log_message("Using existing model output for target: ", target_var)
    saved <- readRDS(model_path)
    edge_table <- readr::read_csv(edge_path, show_col_types = FALSE)
    png_path <- path(
      "outputs",
      "figures",
      paste0("network_uk_", target_var, ".png")
    )
    png_path_main <- path(
      "outputs",
      "figures",
      paste0("network_", target_var, ".png")
    )

    safe_save_rds(saved, model_path_main)
    safe_write_csv(edge_table, edge_path_main)
    if (file.exists(png_path)) {
      file.copy(png_path, png_path_main, overwrite = TRUE)
      log_message("Wrote: ", png_path_main)
    }

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

  vars <- unique(c(target_var, contextual_vars, categorical_vars))
  vars <- intersect(vars, names(df))

  prepared <- prepare_mgm_data(
    df = df,
    vars = vars,
    categorical_vars = categorical_vars
  )

  if (!(target_var %in% prepared$var_names)) {
    stop("Target variable was dropped before modelling due to missingness or no variation: ", target_var)
  }

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
      n_complete = prepared$n
    ),
    model_path
  )

  safe_save_rds(
    list(
      target = target_var,
      model = model,
      var_names = prepared$var_names,
      type = prepared$type,
      level = prepared$level,
      n_complete = prepared$n
    ),
    model_path_main
  )

  edge_table <- extract_edge_table(
    model = model,
    var_names = prepared$var_names,
    target_var = target_var
  )

  safe_write_csv(
    edge_table,
    edge_path
  )

  safe_write_csv(
    edge_table,
    edge_path_main
  )

  png_path <- path(
    "outputs",
    "figures",
    paste0("network_uk_", target_var, ".png")
  )

  png_path_main <- path(
    "outputs",
    "figures",
    paste0("network_", target_var, ".png")
  )

  grDevices::png(png_path, width = 2400, height = 2000, res = 300)
  qgraph::qgraph(
    model$pairwise$wadj,
    labels = prepared$var_names,
    layout = "spring",
    title = paste("UK-wide conditional dependency network:", target_var)
  )
  grDevices::dev.off()

  log_message("Wrote: ", png_path)

  grDevices::png(png_path_main, width = 2400, height = 2000, res = 300)
  qgraph::qgraph(
    model$pairwise$wadj,
    labels = prepared$var_names,
    layout = "spring",
    title = paste("UK-wide conditional dependency network:", target_var)
  )
  grDevices::dev.off()

  log_message("Wrote: ", png_path_main)

  tibble::tibble(
    target = target_var,
    n_complete = prepared$n,
    n_nodes = length(prepared$var_names),
    n_edges = nrow(edge_table),
    model_path = model_path
  )
}

model_summary <- purrr::map_dfr(ethnic_density_vars, fit_one_network)

safe_write_csv(
  model_summary,
  path("outputs", "tables", "network_model_summary_uk.csv")
)

safe_write_csv(
  model_summary,
  path("outputs", "tables", "network_model_summary.csv")
)

log_message("UK-wide network modelling step complete.")

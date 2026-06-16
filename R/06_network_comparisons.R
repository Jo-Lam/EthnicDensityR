# 06_network_comparisons.R
# Compare group-specific network structures

source(file.path("R", "00_setup.R"))

log_message("Starting UK-wide network comparison step.")

edge_files <- list.files(
  DIRS$tables,
  pattern = "^edges_uk_.*\\.csv$",
  full.names = TRUE
)

if (length(edge_files) == 0) {
  stop("No UK-wide edge tables found. Run R/05_network_models.R first.")
}

edges <- purrr::map_dfr(edge_files, function(f) {
  readr::read_csv(f, show_col_types = FALSE) |>
    dplyr::mutate(source_file = basename(f))
})

safe_write_csv(
  edges,
  path("outputs", "tables", "all_network_edges.csv")
)

safe_write_csv(
  edges,
  path("outputs", "tables", "all_network_edges_uk.csv")
)

target_neighbours <- edges |>
  dplyr::filter(.data$involves_target) |>
  dplyr::mutate(
    neighbour = dplyr::if_else(
      .data$node1 == .data$target,
      as.character(.data$node2),
      as.character(.data$node1)
    )
  ) |>
  dplyr::select("target", "neighbour", "weight") |>
  dplyr::arrange(.data$target, dplyr::desc(abs(.data$weight)))

safe_write_csv(
  target_neighbours,
  path("outputs", "tables", "target_node_neighbours.csv")
)

safe_write_csv(
  target_neighbours,
  path("outputs", "tables", "target_node_neighbours_uk.csv")
)

target_neighbours_top5 <- target_neighbours |>
  dplyr::group_by(.data$target) |>
  dplyr::slice_max(
    order_by = abs(.data$weight),
    n = 5,
    with_ties = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(.data$target) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::select("target", "rank", "neighbour", "weight")

safe_write_csv(
  target_neighbours_top5,
  path("outputs", "tables", "table_main_top5_target_neighbours.csv")
)

safe_write_csv(
  target_neighbours_top5,
  path("outputs", "tables", "table_main_top5_target_neighbours_uk.csv")
)

p_top5 <- ggplot2::ggplot(
  target_neighbours_top5,
  ggplot2::aes(
    x = stats::reorder(.data$neighbour, abs(.data$weight)),
    y = abs(.data$weight)
  )
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ target, scales = "free_y") +
  ggplot2::labs(
    title = "Top five target-node neighbours in UK-wide networks",
    subtitle = "Bars show absolute regularised MGM edge weight; retained edges are selected by EBIC regularisation.",
    x = NULL,
    y = "Absolute edge weight"
  ) +
  ggplot2::theme_minimal()

save_ggplot(
  p_top5,
  path("outputs", "figures", "figure_main_top5_target_neighbours.png"),
  width = 12,
  height = 8
)

node_labels <- c(
  share_asian = "Asian density",
  share_bangladeshi = "Bangladeshi density",
  share_black = "Black density",
  share_indian = "Indian density",
  share_mixed = "Mixed density",
  share_other_ethnic = "Other ethnic density",
  share_pakistani = "Pakistani density",
  share_white = "White density",
  share_born_africa = "Africa-born",
  share_born_middle_east_asia = "Middle East/\nAsia-born",
  share_born_uk = "UK-born",
  share_christian = "Christian",
  share_level4_plus = "Degree-level\nqualification",
  share_hindu = "Hindu",
  share_multi_ethnic_households = "Multi-ethnic\nhouseholds",
  share_muslim = "Muslim",
  share_no_car = "No-car\nhouseholds",
  share_no_religion = "No religion",
  share_nssec_high = "High NS-SEC",
  share_older_65_plus = "Aged 65+",
  share_other_religion = "Other religion",
  share_professional_occupations = "Professional\noccupations"
)

target_order <- c(
  "share_white",
  "share_black",
  "share_asian",
  "share_bangladeshi"
)

target_panel_edges <- target_neighbours_top5 |>
  dplyr::filter(.data$target %in% target_order) |>
  dplyr::mutate(
    target = factor(.data$target, levels = target_order),
    target_label = dplyr::recode(
      as.character(.data$target),
      !!!node_labels
    ),
    neighbour_label = dplyr::recode(
      .data$neighbour,
      !!!node_labels,
      .default = .data$neighbour
    ),
    angle = dplyr::case_when(
      .data$rank == 1 ~ 90,
      .data$rank == 2 ~ 30,
      .data$rank == 3 ~ -30,
      .data$rank == 4 ~ -90,
      TRUE ~ 150
    ),
    x = cos(.data$angle * pi / 180),
    y = sin(.data$angle * pi / 180),
    edge_label = sprintf("%.2f", .data$weight)
  )

target_panel_nodes <- dplyr::bind_rows(
  target_panel_edges |>
    dplyr::distinct(.data$target, label = .data$target_label) |>
    dplyr::mutate(
      x = 0,
      y = 0,
      node_type = "Target"
    ),
  target_panel_edges |>
    dplyr::transmute(
      target = .data$target,
      label = .data$neighbour_label,
      x = .data$x,
      y = .data$y,
      node_type = "Neighbour"
    )
)

p_target_network_panel <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = target_panel_edges,
    ggplot2::aes(
      x = 0,
      y = 0,
      xend = .data$x,
      yend = .data$y,
      linewidth = abs(.data$weight)
    ),
    colour = "grey35",
    lineend = "round"
  ) +
  ggplot2::geom_label(
    data = target_panel_edges,
    ggplot2::aes(
      x = .data$x * 0.52,
      y = .data$y * 0.52,
      label = .data$edge_label
    ),
    size = 3,
    linewidth = 0.15,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = target_panel_nodes,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      fill = .data$node_type,
      size = .data$node_type
    ),
    shape = 21,
    colour = "grey20",
    stroke = 0.35
  ) +
  ggplot2::geom_label(
    data = target_panel_nodes,
    ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
    size = 3.2,
    lineheight = 0.9,
    linewidth = 0,
    fill = "white",
    label.padding = ggplot2::unit(0.12, "lines")
  ) +
  ggplot2::facet_wrap(
    ~ target,
    labeller = ggplot2::as_labeller(node_labels),
    ncol = 2
  ) +
  ggplot2::scale_fill_manual(
    values = c(Target = "#1f567d", Neighbour = "#f2c14e")
  ) +
  ggplot2::scale_size_manual(
    values = c(Target = 18, Neighbour = 12)
  ) +
  ggplot2::scale_linewidth(
    range = c(0.5, 3.5)
  ) +
  ggplot2::coord_equal(xlim = c(-1.55, 1.55), ylim = c(-1.35, 1.35)) +
  ggplot2::labs(
    title = "Top five retained target-neighbour edges in selected UK-wide MGM networks",
    subtitle = "Labels show regularised MGM edge weights; panels show the target node and its five strongest retained neighbours.",
    x = NULL,
    y = NULL,
    linewidth = "Edge weight",
    fill = NULL,
    size = NULL
  ) +
  ggplot2::guides(
    fill = "none",
    size = "none"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    legend.position = "bottom",
    panel.spacing = ggplot2::unit(1.2, "lines"),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA)
  )

save_ggplot(
  p_target_network_panel,
  path(
    "outputs",
    "figures",
    "figure_main_top5_network_panel_white_black_asian_bangladeshi.png"
  ),
  width = 10,
  height = 8
)

all_target_order <- c(
  "share_white",
  "share_black",
  "share_asian",
  "share_bangladeshi",
  "share_indian",
  "share_pakistani",
  "share_mixed",
  "share_other_ethnic"
)

all_target_panel_edges <- target_neighbours_top5 |>
  dplyr::filter(.data$target %in% all_target_order) |>
  dplyr::mutate(
    target = factor(.data$target, levels = all_target_order),
    target_label = dplyr::recode(
      as.character(.data$target),
      !!!node_labels
    ),
    neighbour_label = dplyr::recode(
      .data$neighbour,
      !!!node_labels,
      .default = .data$neighbour
    ),
    angle = dplyr::case_when(
      .data$rank == 1 ~ 90,
      .data$rank == 2 ~ 30,
      .data$rank == 3 ~ -30,
      .data$rank == 4 ~ -90,
      TRUE ~ 150
    ),
    x = cos(.data$angle * pi / 180),
    y = sin(.data$angle * pi / 180),
    edge_label = sprintf("%.2f", .data$weight)
  )

all_target_panel_nodes <- dplyr::bind_rows(
  all_target_panel_edges |>
    dplyr::distinct(.data$target, label = .data$target_label) |>
    dplyr::mutate(
      x = 0,
      y = 0,
      node_type = "Target"
    ),
  all_target_panel_edges |>
    dplyr::transmute(
      target = .data$target,
      label = .data$neighbour_label,
      x = .data$x,
      y = .data$y,
      node_type = "Neighbour"
    )
)

p_all_target_network_panel <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = all_target_panel_edges,
    ggplot2::aes(
      x = 0,
      y = 0,
      xend = .data$x,
      yend = .data$y,
      linewidth = abs(.data$weight)
    ),
    colour = "grey35",
    lineend = "round"
  ) +
  ggplot2::geom_label(
    data = all_target_panel_edges,
    ggplot2::aes(
      x = .data$x * 0.52,
      y = .data$y * 0.52,
      label = .data$edge_label
    ),
    size = 2.6,
    linewidth = 0.15,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = all_target_panel_nodes,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      fill = .data$node_type,
      size = .data$node_type
    ),
    shape = 21,
    colour = "grey20",
    stroke = 0.35
  ) +
  ggplot2::geom_label(
    data = all_target_panel_nodes,
    ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
    size = 2.7,
    lineheight = 0.88,
    linewidth = 0,
    fill = "white",
    label.padding = ggplot2::unit(0.10, "lines")
  ) +
  ggplot2::facet_wrap(
    ~ target,
    labeller = ggplot2::as_labeller(node_labels),
    ncol = 4
  ) +
  ggplot2::scale_fill_manual(
    values = c(Target = "#1f567d", Neighbour = "#f2c14e")
  ) +
  ggplot2::scale_size_manual(
    values = c(Target = 14, Neighbour = 10)
  ) +
  ggplot2::scale_linewidth(
    range = c(0.45, 3.2)
  ) +
  ggplot2::coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.25, 1.25)) +
  ggplot2::labs(
    title = "Top five retained target-neighbour edges in UK-wide MGM networks",
    subtitle = "Labels show regularised MGM edge weights; panels show each ethnic-density target and its five strongest retained neighbours.",
    x = NULL,
    y = NULL,
    linewidth = "Edge weight",
    fill = NULL,
    size = NULL
  ) +
  ggplot2::guides(
    fill = "none",
    size = "none"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    legend.position = "bottom",
    panel.spacing = ggplot2::unit(0.9, "lines"),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA)
  )

save_ggplot(
  p_all_target_network_panel,
  path(
    "outputs",
    "figures",
    "figure_main_top5_network_panel_all_groups.png"
  ),
  width = 14,
  height = 8.5
)

save_ggplot(
  p_all_target_network_panel,
  path(
    "outputs",
    "figures",
    "figure_main_top5_network_panel_all_groups.pdf"
  ),
  width = 14,
  height = 8.5
)

save_ggplot(
  p_target_network_panel,
  path(
    "outputs",
    "figures",
    "figure_main_top5_network_panel_white_black_asian_bangladeshi.pdf"
  ),
  width = 10,
  height = 8
)

sensitivity_thresholds <- c(0.05, 0.10)

threshold_sensitivity <- purrr::map_dfr(
  sensitivity_thresholds,
  function(threshold) {
    target_neighbours |>
      dplyr::filter(abs(.data$weight) >= threshold) |>
      dplyr::mutate(threshold = threshold)
  }
)

safe_write_csv(
  threshold_sensitivity,
  path("outputs", "tables", "supp_threshold_target_neighbours.csv")
)

safe_write_csv(
  threshold_sensitivity,
  path("outputs", "tables", "supp_threshold_target_neighbours_uk.csv")
)

threshold_counts <- threshold_sensitivity |>
  dplyr::count(.data$threshold, .data$target, name = "n_edges")

safe_write_csv(
  threshold_counts,
  path("outputs", "tables", "supp_threshold_target_neighbour_counts.csv")
)

safe_write_csv(
  threshold_counts,
  path("outputs", "tables", "supp_threshold_target_neighbour_counts_uk.csv")
)

p_threshold <- ggplot2::ggplot(
  threshold_counts,
  ggplot2::aes(
    x = .data$target,
    y = .data$n_edges,
    fill = factor(.data$threshold)
  )
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::labs(
    title = "Sensitivity: target-neighbour counts under edge-weight thresholds",
    x = NULL,
    y = "Number of target-node neighbours",
    fill = "|weight| threshold"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

save_ggplot(
  p_threshold,
  path("outputs", "figures", "supp_threshold_target_neighbour_counts.png"),
  width = 9,
  height = 6
)

# Edge-set similarity across ethnic-density networks.
edge_sets <- edges |>
  dplyr::mutate(
    edge_id = purrr::map2_chr(
      as.character(.data$node1),
      as.character(.data$node2),
      ~ paste(sort(c(.x, .y)), collapse = " -- ")
    )
  ) |>
  dplyr::group_by(.data$target) |>
  dplyr::summarise(
    edge_set = list(unique(.data$edge_id)),
    .groups = "drop"
  )

targets <- edge_sets$target

jaccard <- expand.grid(
  target1 = targets,
  target2 = targets,
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::rowwise() |>
  dplyr::mutate(
    similarity = {
      a <- edge_sets$edge_set[[match(.data$target1, edge_sets$target)]]
      b <- edge_sets$edge_set[[match(.data$target2, edge_sets$target)]]
      length(intersect(a, b)) / length(union(a, b))
    }
  ) |>
  dplyr::ungroup()

safe_write_csv(
  jaccard,
  path("outputs", "tables", "network_edge_jaccard_similarity.csv")
)

safe_write_csv(
  jaccard,
  path("outputs", "tables", "network_edge_jaccard_similarity_uk.csv")
)

p_jaccard <- ggplot2::ggplot(
  jaccard,
  ggplot2::aes(x = .data$target1, y = .data$target2, fill = .data$similarity)
) +
  ggplot2::geom_tile() +
  ggplot2::labs(
    title = "Network edge-set similarity across ethnic-density measures",
    x = NULL,
    y = NULL,
    fill = "Jaccard\nsimilarity"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

save_ggplot(
  p_jaccard,
  path("outputs", "figures", "network_edge_jaccard_similarity.png"),
  width = 8,
  height = 7
)

# Wide table of target-node neighbour weights.
target_neighbour_wide <- target_neighbours |>
  dplyr::mutate(weight = round(.data$weight, 4)) |>
  tidyr::pivot_wider(
    names_from = "target",
    values_from = "weight"
  )

safe_write_csv(
  target_neighbour_wide,
  path("outputs", "tables", "target_node_neighbour_weights_wide.csv")
)

safe_write_csv(
  target_neighbour_wide,
  path("outputs", "tables", "target_node_neighbour_weights_wide_uk.csv")
)

log_message("UK-wide network comparison step complete.")

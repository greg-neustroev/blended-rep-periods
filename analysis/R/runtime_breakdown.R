# runtime_breakdown.R -- Figure 3 (runtime_breakdown): wall-clock time (s) per pipeline stage vs.
# n_rep_periods, for all four realistic case studies (5-bus, GEP, P2X, 118-bus).
#
# One block per case study, stacked top to bottom, EACH case study named as a header above its own
# block (not a rotated side strip) -- block height is proportional to sqrt(that case study's own
# max time), so it is visually obvious 5-bus is far faster than GEP/P2X/118-bus rather than every
# block reading as equally "tall" under its own independently-scaled axis. Within a block, one
# panel per clustering method (k-means, Hierarchical, Chronological, Convex hull; horizontal strip
# labels on top), and within each panel, one pair of adjacent stacked columns per n_rep_periods
# value -- Dirac and Convex weights side by side -- each stacked by pipeline stage
# (read/preprocess, cluster, fit weights, formulate model, solve). A dashed horizontal line marks
# that case study's full-horizon (n_rep=1) reference solve time, so the reader can see directly
# whether the reduced model is actually cheaper.
#
# Source data: analysis/output/data/runtime_breakdown.csv, IN THIS REPO (produced by
# analysis/paper/runtime_breakdown.jl from the raw per-run outputs). This script is fully
# self-contained within the experiments repo -- it never reads from or writes into the sibling
# `clustering` repo; copying the finished figure over is a separate, manual step.
#
# Output: analysis/output/figures/runtime_breakdown.pdf (+ a PNG preview).

library(ggplot2)
library(dplyr)
library(patchwork)

# Bootstrap: walk up from the current working directory to find the repo root (marked by
# Project.toml), then source common.R from there. Deliberately NOT based on introspecting this
# script's own file path (`commandArgs()`'s `--file=`, or scanning `sys.frames()` for an `ofile`)
# -- that approach is runner-dependent and broke under at least one editor's "Run file" mechanism
# with "Could not determine this script's own file path". Only requires the R process's working
# directory to be at or below the repo root, true for every normal way of launching R here.
.repo_root_bootstrap <- local({
  dir <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "Project.toml"))) break
    parent <- dirname(dir)
    if (parent == dir) {
      stop("Could not find the repo root (looked for Project.toml walking up from ", getwd(), ")")
    }
    dir <- parent
  }
  dir
})
source(file.path(.repo_root_bootstrap, "analysis", "R", "common.R"))

df <- read.csv(data_path("runtime_breakdown.csv"))

case_titles <- c(`5bus` = "5-bus", gep = "GEP", p2x = "P2X", `118bus` = "118-bus")
case_order <- names(case_titles)

clustering_names <- c(k_means = "k-means", hierarchical = "Hierarchical",
                       chronological = "Chronological", convex_hull = "Convex hull")
clustering_order <- names(clustering_names)

weight_names <- c(dirac = "Dirac", convex = "Convex")
weight_order <- names(weight_names)

stage_names <- c(
  read_preprocess = "Read/preprocess",
  cluster = "Cluster",
  fit_weights = "Fit weights",
  formulate_model = "Formulate",
  solve = "Solve"
)
stage_order <- names(stage_names)

# Viridis (perceptually uniform, colorblind-safe, and stays distinguishable when printed in
# grayscale) for the 5 pipeline stages -- a different encoding axis from regret.R's
# clustering-method/weight-type palettes, so this figure never borrows a hue with a different
# meaning elsewhere in the paper.
stage_colors <- viridisLite::viridis(5, end = 0.92)

panels <- df %>%
  filter(clustering_type %in% clustering_order, weight_type %in% weight_order, stage %in% stage_order) %>%
  mutate(
    clustering_label = factor(clustering_names[clustering_type], levels = unname(clustering_names)),
    weight_label = factor(weight_names[weight_type], levels = unname(weight_names)),
    stage_label = factor(stage_names[stage], levels = rev(unname(stage_names))),
    n_rep_factor = factor(n_rep_periods, levels = c(10, 20, 40, 80))
  )
# x-axis: one pair of adjacent bars per n_rep value -- "D"/"C" (Dirac/Convex) tick labels, built
# directly as "n_rep letter" strings (NOT via interaction(), whose default arg order put weight
# first and silently produced all-NA levels against the n_rep-first level list below; NOT via
# ggh4x::guide_axis_nested() either -- deprecated as of ggh4x 0.3.0 and its replacement,
# legendry::guide_axis_nested(), isn't installed here -- so "# RP" is a plain shared axis title
# rather than a true nested per-pair label). A blank spacer level between consecutive n_rep groups
# (an empty factor level with no data, excluded from `breaks` so it draws no label of its own)
# opens a visible gap there, so the four pairs read as distinct groups rather than one
# undifferentiated run of 8 bars.
weight_letter <- c(Dirac = "D", Convex = "C")
panels$x_group <- paste(panels$n_rep_factor, weight_letter[as.character(panels$weight_label)])
n_rep_levels <- levels(panels$n_rep_factor)
x_levels <- unlist(lapply(seq_along(n_rep_levels), function(i) {
  pair <- paste(n_rep_levels[i], unname(weight_letter[levels(panels$weight_label)]))
  if (i < length(n_rep_levels)) c(pair, paste0(".spacer", i)) else pair
}))
panels$x_group <- factor(panels$x_group, levels = x_levels)

ref_lines <- df %>%
  filter(method_label == "full_reference") %>%
  transmute(case_study, ref_time_s = time_mean_s)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(8, "pt"),
    plot.margin = margin(2, 8, 2, 2),
    strip.text = element_text(face = "bold", size = 11, angle = 0),
    plot.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 11),
    axis.text.x = element_text(size = 8),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

y_labels <- function(x) format(round(x, 1), scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

panels <- panels %>% left_join(ref_lines, by = "case_study")

# One block per case study, case-study name as a plain header ABOVE its own 4-column strip (not a
# rotated side label); each block's own y-scale (via being its own ggplot, not a shared facet_grid
# axis) so the four case studies' very different time scales don't force one one shared axis. Only
# 3 y-breaks are shown -- 0 (with a solid black axis line, via base_theme's axis.line.x), the
# dashed full-horizon reference time, and the tallest ACTUAL METHOD bar (bar_max_time_s, NOT
# unioned with the reference -- the reference sits above the tallest bar in some blocks (5-bus,
# GEP/P2X/118-bus at low n_rep) and below it in others (GEP/P2X at n_rep=80), so both need their
# own break reported independently rather than collapsing to whichever happens to be larger) --
# rather than an arbitrary evenly-spaced tick sequence, since those are the only three values the
# figure is actually meant to be read off of (is a bar above or below the reference line, and how
# far below).
make_block <- function(case, show_x_title, ref_time_s, bar_max_time_s) {
  d <- panels %>% filter(case_study == case)
  y_breaks <- sort(unique(c(0, ref_time_s, bar_max_time_s)))
  ggplot(d, aes(x = x_group, y = time_mean_s, fill = stage_label)) +
    geom_col(width = 0.8) +
    geom_hline(aes(yintercept = ref_time_s), linetype = "dashed", linewidth = 0.6, color = "grey30",
               na.rm = TRUE) +
    facet_wrap(~clustering_label, nrow = 1) +
    scale_fill_manual(name = "Pipeline stage", values = stage_colors, breaks = unname(stage_names)) +
    scale_y_continuous(name = "Time [s]", breaks = y_breaks, labels = y_labels,
                       expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(name = if (show_x_title) "# RP, weight type (D/C)" else NULL,
                     breaks = x_levels[!grepl("^\\.spacer", x_levels)], drop = FALSE) +
    labs(title = case_titles[[case]]) +
    base_theme
}

# each block's own max STACKED bar height (sum of stage times per x_group x clustering cell, not
# the per-stage time_mean_s rows themselves) and reference time, needed both for the 3-value
# y-breaks (bar max and reference kept separate, see make_block) and for the proportional block
# heights below (which DO need the overall visual extent -- bars unioned with the reference line --
# since the panel must be tall enough to show the dashed line even where it sits above every bar).
stacked_totals <- panels %>%
  group_by(case_study, clustering_label, x_group) %>%
  summarise(total_s = sum(time_mean_s), .groups = "drop")
block_ref <- setNames(ref_lines$ref_time_s, ref_lines$case_study)[case_order]
block_bar_max <- sapply(case_order, function(cs) max(stacked_totals$total_s[stacked_totals$case_study == cs], na.rm = TRUE))
block_panel_max <- pmax(block_bar_max, block_ref, na.rm = TRUE)

blocks <- lapply(seq_along(case_order), function(i) {
  cs <- case_order[i]
  make_block(cs, i == length(case_order), block_ref[[cs]], block_bar_max[[cs]])
})

# Block heights proportional to sqrt(that case study's own max visual extent) -- a plain linear
# proportion would make 5-bus's block illegibly thin next to GEP's (max time spans ~50x across
# case studies); sqrt compresses that range while still making the "5-bus is far faster" and
# "GEP takes far longer" comparisons visually obvious at a glance, which a free, independently
# auto-scaled axis per block would otherwise hide (every block would just look equally "tall").
block_heights <- sqrt(block_panel_max)

final <- wrap_plots(blocks, ncol = 1, heights = block_heights) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(figure_path("runtime_breakdown.pdf"), final, width = 9, height = 9, device = cairo_pdf, limitsize = FALSE)
ggsave(figure_path("runtime_breakdown_preview.png"), final, width = 9, height = 9, dpi = 150, limitsize = FALSE)

# runtime_breakdown.R -- Figure 3 (runtime_breakdown): wall-clock time (s) per pipeline stage vs.
# n_rep_periods, for all four realistic case studies (5-bus, GEP, P2X, 118-bus).
#
# One grid: rows = case study, columns = clustering method (k-means, Hierarchical, Chronological,
# Convex hull). Within each panel, one pair of adjacent stacked columns per n_rep_periods value --
# Dirac and Convex weights side by side -- each stacked by pipeline stage (read/preprocess,
# cluster, fit weights, formulate model, solve). A dashed horizontal line marks that case study's
# full-horizon (n_rep=1) reference solve time, so the reader can see directly whether the reduced
# model is actually cheaper.
#
# Source data: analysis/output/data/runtime_breakdown.csv, IN THIS REPO (produced by
# analysis/paper/runtime_breakdown.jl from the raw per-run outputs). This script is fully
# self-contained within the experiments repo -- it never reads from or writes into the sibling
# `clustering` repo; copying the finished figure over is a separate, manual step.
#
# Output: analysis/output/figures/runtime_breakdown.pdf (+ a PNG preview).

library(ggplot2)
library(dplyr)

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
    case_label = factor(case_titles[case_study], levels = unname(case_titles)),
    clustering_label = factor(clustering_names[clustering_type], levels = unname(clustering_names)),
    weight_label = factor(weight_names[weight_type], levels = unname(weight_names)),
    stage_label = factor(stage_names[stage], levels = rev(unname(stage_names))),
    n_rep_factor = factor(n_rep_periods, levels = c(10, 20, 40, 80))
  )
# x-axis: one pair of adjacent bars (Dirac, Convex) per n_rep value -- built directly as
# "n_rep\nweight" strings (not via interaction(), whose default arg order put weight first and
# silently produced all-NA levels against the n_rep-first level list below) so each n_rep's pair
# stays adjacent and the four pairs read left to right in increasing n_rep order.
panels$x_group <- paste(panels$n_rep_factor, panels$weight_label, sep = "\n")
x_levels <- as.vector(t(outer(levels(panels$n_rep_factor), levels(panels$weight_label),
                               FUN = function(n, w) paste(n, w, sep = "\n"))))
panels$x_group <- factor(panels$x_group, levels = x_levels)

ref_lines <- df %>%
  filter(method_label == "full_reference") %>%
  transmute(case_label = factor(case_titles[case_study], levels = unname(case_titles)), ref_time_s = time_mean_s)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing = unit(8, "pt"),
    plot.margin = margin(2, 2, 2, 2),
    strip.text = element_text(face = "bold", size = 11),
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

y_labels <- function(x) format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

panels <- panels %>% left_join(ref_lines, by = "case_label")

p <- ggplot(panels, aes(x = x_group, y = time_mean_s, fill = stage_label)) +
  geom_col(width = 0.8) +
  geom_hline(aes(yintercept = ref_time_s), linetype = "dashed", linewidth = 0.6, color = "grey30",
             na.rm = TRUE) +
  facet_grid(case_label ~ clustering_label, scales = "free_y", switch = "y") +
  scale_fill_manual(name = "Pipeline stage", values = stage_colors, breaks = unname(stage_names)) +
  scale_y_continuous(name = NULL, labels = y_labels) +
  scale_x_discrete(name = "# RP, weight type") +
  base_theme

ggsave(figure_path("runtime_breakdown.pdf"), p, width = 11, height = 11, device = cairo_pdf, limitsize = FALSE)
ggsave(figure_path("runtime_breakdown_preview.png"), p, width = 11, height = 11, dpi = 150, limitsize = FALSE)

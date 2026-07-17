# runtime_breakdown.R -- Figure 3 (runtime_breakdown): wall-clock time (s) per pipeline stage vs.
# n_rep_periods, for all four realistic case studies (5-bus, GEP, P2X, 118-bus).
#
# Each case study is one 2-row block of panels, one panel per method (full clustering+weight names,
# not single-letter codes): PROPOSED (conical hull + convex weights) against a representative
# baseline spread that varies clustering method and weight type independently against it --
# k-means/Hierarchical/Chronological at Dirac (the conventional pairing), k-means/Hierarchical
# re-paired with PROPOSED's own (convex) weight type, and PROPOSED's own clustering method (conic
# hull) re-paired with Dirac -- see analysis/paper/runtime_breakdown.jl's RUNTIME_METHODS for the
# exact list. Within a panel, bars are stacked by pipeline stage (read/preprocess, cluster, fit
# weights, formulate model, solve) across n_rep_periods in {10,20,40,80}; a dashed horizontal line
# marks that case study's full-horizon (n_rep=1) reference solve time, so the reader can see
# directly whether the reduced model is actually cheaper.
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

# Full clustering+weight names (not single-letter codes), in the display order set by
# RUNTIME_METHODS in analysis/paper/runtime_breakdown.jl.
method_order <- c("k-means, Dirac", "Hierarchical, Dirac", "Chronological, Dirac", "Conic, Dirac",
                   "k-means, Convex", "Hierarchical, Convex", 
                   "Conic, Convex (proposed)")

stage_names <- c(
  read_preprocess = "Read/preprocess",
  cluster = "Cluster",
  fit_weights = "Fit weights",
  formulate_model = "Formulate",
  solve = "Solve"
)
stage_order <- names(stage_names)

# Okabe-Ito-derived 5-color subset (colorblind-safe), one per pipeline stage -- a different
# encoding axis from regret.R's clustering-method/weight-type palettes, so this figure never
# borrows a hue with a different meaning elsewhere in the paper.
stage_colors <- c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00")

panels <- df %>%
  filter(method_label %in% method_order, stage %in% stage_order) %>%
  mutate(
    case_label = factor(case_titles[case_study], levels = unname(case_titles)),
    method_label_disp = factor(method_label, levels = method_order),
    stage_label = factor(stage_names[stage], levels = rev(unname(stage_names))),
    n_rep_factor = factor(n_rep_periods, levels = c(10, 20, 40, 80))
  )

ref_lines <- df %>%
  filter(method_label == "full_reference") %>%
  transmute(case_label = factor(case_titles[case_study], levels = unname(case_titles)),
            ref_time_s = time_mean_s)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing = unit(10, "pt"),
    plot.margin = margin(2, 2, 2, 2),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

y_labels <- function(x) format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

# One 2-row block per case study (own y-axis scale, since total time spans two-plus orders of
# magnitude across case studies), one panel per method (full clustering+weight name, wrapped 4+3
# across two rows so the longer labels stay legible); bars stacked by pipeline stage. The n_rep=1
# full-horizon reference time is a dashed horizontal line spanning every method panel in that case
# study's block, so the "does reduction actually save time" comparison is direct.
make_row <- function(case) {
  d <- panels %>% filter(case_study == case)
  ref <- ref_lines %>% filter(case_label == case_titles[[case]])
  ref_y <- if (nrow(ref) > 0) ref$ref_time_s[1] else NA_real_

  ggplot(d, aes(x = n_rep_factor, y = time_mean_s, fill = stage_label)) +
    geom_col(width = 0.7) +
    { if (!is.na(ref_y)) geom_hline(yintercept = ref_y, linetype = "dashed", linewidth = 0.6, color = "grey30") } +
    facet_wrap(~method_label_disp, nrow = 2) +
    scale_fill_manual(name = "Pipeline stage", values = stage_colors, breaks = unname(stage_names)) +
    scale_y_continuous(name = paste0(case_titles[[case]], "\ntime [s]"), labels = y_labels) +
    scale_x_discrete(name = if (case == case_order[length(case_order)]) "# RP" else NULL) +
    base_theme
}

rows <- lapply(case_order, make_row)

# One shared legend (pipeline-stage fill), collected via patchwork's plot_layout(guides="collect")
# -- exactly the mechanism regret.R already relies on -- and de-duplicated to one copy at the
# bottom, rather than repeated once per case-study block.
final <- wrap_plots(rows, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(figure_path("runtime_breakdown.pdf"), final, width = 10, height = 18, device = cairo_pdf, limitsize = FALSE)
ggsave(figure_path("runtime_breakdown_preview.png"), final, width = 10, height = 18, dpi = 150, limitsize = FALSE)

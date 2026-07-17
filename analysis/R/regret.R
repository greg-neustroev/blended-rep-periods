# regret.R -- Figure 2 (central_regret_by_nrep): regret (%) vs. n_rep, economic normalization,
# rendered separately for each of the four realistic case studies (GEP, 5-bus, P2X, 118-bus).
#
# For each case study, two strips show every clustering x weight combination (28 cells), not just
# a single representative slice of the other axis -- the point is to let the reader see for
# themselves whether clustering method or weight type drives regret, and confirm it holds
# regardless of what the other axis is set to, rather than asserting it from one fixed pairing.
#   Strip A: facet by weight type (4 panels); each panel's lines are the 7 clustering methods.
#   Strip B: facet by clustering method (7 panels); each panel's lines are the 4 weight types.
# Stochastic clustering methods (k-means, k-medoids) get +-1 SD error bars over the 5 seeds, and
# propagate that same stochasticity to every weight type fit on top of their seed-dependent RPs;
# deterministic methods (regret_sd_pct == 0) show as a plain line+point, per R1g/R5.19/R5.20.
# Distinct color (colorblind-safe Okabe-Ito) + linetype + shape per series, thick lines, per R4.4.
#
# Source data: analysis/output/data/regret_summary.csv, IN THIS REPO (produced by
# analysis/paper/regret.jl from the raw per-run outputs; already the tidy, per-clustering x weight
# x n_rep x case_study x normalization summary with 5-seed mean/SD used throughout the paper's
# tables). This script is fully self-contained within the experiments repo -- it never reads from
# or writes into the sibling `clustering` repo; copying the finished figures over is a separate,
# manual step.
#
# Output: analysis/output/figures/central_regret_by_nrep_<case>.pdf, one per case study.

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggh4x)

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

df <- read.csv(data_path("regret_summary.csv"))
df <- df %>%
  filter(normalization == "economic") %>%
  mutate(
    n_rep_periods = as.numeric(n_rep_periods),
    is_deterministic = tolower(as.character(is_deterministic)) == "true",
    pareto_frontier = tolower(as.character(pareto_frontier)) == "true"
  )

clustering_names <- c(
  k_means = "k-Means",
  k_medoids = "k-Medoids",
  hierarchical = "Hierarchical",
  chronological = "Chronological",
  convex_hull = "Convex",
  convex_hull_with_null = "Convex with 0",
  conical_hull = "Conic"
)
clustering_levels <- names(clustering_names)

weight_names <- c(
  dirac = "Dirac",
  convex = "Convex",
  conical_bounded = "Sub-unit",
  conical = "Conic"
)
weight_levels <- names(weight_names)

case_titles <- c(
  `5bus` = "5-bus", gep = "GEP",
  p2x = "P2X", `118bus` = "118-bus"
)

# Okabe-Ito colorblind-safe palette (7 colors, one per clustering method -- all of Okabe-Ito's
# non-black hues are used here, so weight type below deliberately draws from a different family
# rather than reusing/stretching the same 7 hues or falling back to a narrow black-to-grey ramp
# (whose middle two shades, e.g. grey60/grey80, read as nearly identical at small sizes/print).
okabe_ito <- c("#E69F00", "#D55E00", "#817066", "#CC79A7", "#56B4E9", "#009E73", "#0072B2")
clustering_shapes <- c(16, 17, 15, 3, 7, 8, 18)
# 7 fully distinct linetypes (no repeats) -- a repeated linetype value breaks ggplot's
# legend merging across the color/linetype/shape scales and produces a duplicated legend.
clustering_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1213")

# Weight type uses a 4-color viridis subset: colorblind-safe like Okabe-Ito, but a visually
# distinct family (perceptually-uniform purple-to-yellow) so the two legends are never confused
# with each other at a glance, unlike reusing Okabe-Ito hues or a low-contrast grey ramp.
weight_colors <- viridisLite::viridis(4, end = 0.85)
weight_shapes <- c(16, 17, 15, 18)
weight_linetypes <- c("solid", "dashed", "dotted", "dotdash")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    # legend.position is NOT "none" here: individual panels keep their real (collectible) guide,
    # and plot_layout(guides = "collect") on the combined grid below deduplicates the 4 repeated
    # "Clustering" guides (one per strip_a instance) and 4 repeated "Weight" guides (one per
    # strip_b instance) down to one of each, shown once at the bottom of the whole figure --
    # no separate legend image needed.
    panel.grid.minor = element_blank(),
    panel.spacing = unit(20, "pt"),
    plot.margin = margin(2, 2, 2, 2),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 18),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

y_labels <- function(x) format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

REGRET_CAP <- 100  # fallback y-axis ceiling for the Pareto-front plot below (central_pareto_overview.pdf)

# Strip A: facet by weight type, lines/colors = clustering method (all 7).
# Regret is NOT clipped in the data. Each weight-type facet gets its own y-axis upper limit (via
# ggh4x::facetted_pos_scales, since plain coord_cartesian()/scale limits apply to the whole plot,
# not per facet).
# A rank-based cap (e.g. "1.5x the 3rd-worst method's value") was tried and rejected: it assumes
# only 2 of 7 methods are ever bad at N_RP=10, but for 5-bus/GEP specifically, 3-4 of the 4
# baselines (k-means, k-medoids, hierarchical, chronological) are ALL bad at N_RP=10 (not just
# k-means) -- so a rank-based cap ends up dominated by a baseline method anyway, still squashing
# the paper's own proposed hull methods together near zero. Instead, the cap is 1.5x the WORST
# (highest-regret) of just the three hull methods (convex, convex-with-null, conical) at N_RP=10,
# with a floor of 10 so a tame facet is never zoomed in past the 0-10% range: this guarantees all
# three hull methods are always comfortably visible with margin, directly matching what the figure
# needs to show (Section~\ref{sec:results}: the hull methods vs. the four baselines), regardless
# of how many baselines happen to be bad in a given facet -- baselines worse than this cap simply
# run off the top of the panel (line ends where it exits), rather than being flattened to the cap
# value as if that were their true regret.
HULL_CLUSTERING_TYPES <- c("convex_hull", "convex_hull_with_null", "conical_hull")
cap_from_hull_methods <- function(d) {
  at10 <- d$regret_mean_pct[d$n_rep_periods == 10 & d$clustering_type %in% HULL_CLUSTERING_TYPES]
  if (length(at10) == 0) return(max(d$regret_mean_pct))
  max(1.5 * max(at10), 10)
}

make_strip_a <- function(case, show_y_title = TRUE) {
  d <- df %>% filter(case_study == case)
  d$clustering_label <- factor(clustering_names[d$clustering_type], levels = clustering_names)
  d$weight_label <- factor(weight_names[d$weight_type], levels = weight_names)
  d$n_rep_factor <- factor(d$n_rep_periods, levels = c(10, 20, 40, 80))
  d$ymin <- pmax(d$regret_mean_pct - d$regret_sd_pct, 1e-3)
  d$ymax <- d$regret_mean_pct + d$regret_sd_pct

  # One scale_y_continuous per weight_label level, in factor-level order, each with its own
  # data-driven upper limit; facetted_pos_scales applies them positionally to facet_wrap's panels.
  y_scales <- lapply(levels(d$weight_label), function(w) {
    cap <- cap_from_hull_methods(d[d$weight_label == w, ])
    scale_y_continuous(labels = y_labels, limits = c(0, cap))
  })

  dodge <- position_dodge(width = 0.5)

  # geom_errorbar keeps every row (not just the stochastic ones) so position_dodge sees the same
  # 7 groups in every layer; deterministic rows are drawn fully transparent (zero-height anyway,
  # since their SD is 0) rather than dropped, which would otherwise change the dodge group count
  # for that layer alone and misalign its bars against the (7-group-dodged) lines/points.
  # x-axis title suppressed unconditionally: strip_a is never the bottom row of the combined grid
  # (strip_b always is), so its "# RP" title would just repeat.
  ggplot(d, aes(x = n_rep_factor, y = regret_mean_pct, color = clustering_label,
                linetype = clustering_label, shape = clustering_label, group = clustering_label)) +
    facet_wrap(~weight_label, ncol = 1, scales = "free_y") +
    facetted_pos_scales(y = y_scales) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax, alpha = is_deterministic),
      width = 0.3, linewidth = 0.8, linetype = "solid",
      position = dodge, show.legend = FALSE
    ) +
    geom_line(linewidth = 0.5, position = dodge) +
    geom_point(size = 3, position = dodge) +
    scale_alpha_manual(values = c(`TRUE` = 0, `FALSE` = 1), guide = "none") +
    scale_x_discrete(name = NULL) +
    ylab(if (show_y_title) "Regret [%]" else NULL) +
    scale_color_manual(name = "Clustering", values = okabe_ito, breaks = unname(clustering_names)) +
    scale_linetype_manual(name = "Clustering", values = clustering_linetypes, breaks = unname(clustering_names)) +
    scale_shape_manual(name = "Clustering", values = clustering_shapes, breaks = unname(clustering_names)) +
    labs(title = case_titles[[case]]) +
    base_theme
}

# Strip B: facet by clustering method, lines/colors = weight type (all 4).
# x-axis title kept: strip_b is always the bottom row of the combined grid, for every column.
make_strip_b <- function(case, show_y_title = TRUE) {
  d <- df %>% filter(case_study == case)
  d$clustering_label <- factor(clustering_names[d$clustering_type], levels = clustering_names)
  d$weight_label <- factor(weight_names[d$weight_type], levels = weight_names)
  d$n_rep_factor <- factor(d$n_rep_periods, levels = c(10, 20, 40, 80))
  d$ymin <- pmax(d$regret_mean_pct - d$regret_sd_pct, 1e-3)
  d$ymax <- d$regret_mean_pct + d$regret_sd_pct

  dodge <- position_dodge(width = 0.5)

  ggplot(d, aes(x = n_rep_factor, y = regret_mean_pct, color = weight_label,
                linetype = weight_label, shape = weight_label, group = weight_label)) +
    facet_wrap(~clustering_label, ncol = 1, scales = "free_y") +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax, alpha = is_deterministic),
      width = 0.25, linewidth = 0.8, linetype = "solid",
      position = dodge, show.legend = FALSE
    ) +
    #geom_line(linewidth = 1.05, position = dodge) +
    geom_point(size = 3, position = dodge) +
    scale_alpha_manual(values = c(`TRUE` = 0, `FALSE` = 1), guide = "none") +
    scale_x_discrete(name = "# RP") +
    scale_y_continuous(name = if (show_y_title) "Regret [%]" else NULL, labels = y_labels) +
    scale_color_manual(name = "Weight", values = weight_colors, breaks = unname(weight_names)) +
    scale_linetype_manual(name = "Weight", values = weight_linetypes, breaks = unname(weight_names)) +
    scale_shape_manual(name = "Weight", values = weight_shapes, breaks = unname(weight_names)) +
    base_theme
}

# One combined image: all 4 case studies stacked vertically (one row each), each row itself the
# strip_a/strip_b pair for that case study. No legends anywhere in this image -- see
# make_combined_legend() below for the single shared (one-line) legend, saved separately so it is
# not repeated 4x. This is tall by construction (four full case studies' worth of detail); it is
# meant to be placed as a dedicated full-page figure, split top/bottom via \includegraphics
# trim=/clip= in LaTeX rather than shrunk to fit one page, which would make it illegible.
# (The old separate central_regret_by_nrep_all.pdf detailed spread and the standalone
# central_regret_legend.pdf are no longer generated: the transposed grid below embeds its own
# legend directly via plot_layout(guides = "collect"), so a separate legend image is redundant.)

# ---------------------------------------------------------------------------------------------
# One-page grid, TRANSPOSED: 4 columns (one per case study, in the order GEP, P2X, 5-bus, 118-bus)
# x 11 rows (the first 4 are make_strip_a's weight-type facets, now stacked vertically: clustering
# varies, weight fixed per row; the last 7 are make_strip_b's clustering-type facets, also stacked
# vertically: weight varies, clustering fixed per row), with a small gap row between the two row
# groups. make_strip_a/make_strip_b themselves now facet with ncol=1 (stacked), not nrow=1/ncol=7
# (side by side), to match this transpose.
GRID_COL_ORDER <- c("5bus", "gep", "p2x", "118bus")

# Two block-heading banners and a divider make the row1-4/row5-11 switch in what varies explicit
# rather than left for the reader to infer from row names alone -- each is a blank (theme_void())
# ggplot spanning all 4 columns, placed via patchwork::area() so a single element can occupy a
# full row of the outer grid; this still keeps everything in ONE flat element list passed to ONE
# top-level wrap_plots(design=...) call (no nested plot_layout()), avoiding the same patchwork
# 1.3.2 nesting bug noted above.
make_block_banner <- function(text) {
  ggplot() +
    xlim(-1, 1) + ylim(-1, 1) +
    annotate("text", x = 0, y = 0, label = text, fontface = "bold", size = 4.8) +
    theme_void()
}
make_divider <- function() {
  ggplot() +
    xlim(0, 1) + ylim(0, 1) +
    annotate("segment", x = 0, xend = 1, y = 0.5, yend = 0.5, linewidth = 0.7, color = "grey40") +
    theme_void()
}

# show_y_title = TRUE only for the first column (GEP) -- the other three columns' y-axis titles
# would just repeat "Regret [%]" with no new information (per-panel scales="free_y" still gives
# each its own tick numbers, only the redundant title text is dropped).
grid_elements <- c(
  list(make_block_banner("Clustering varies, conditioned on weight type")),
  lapply(GRID_COL_ORDER, function(case) make_strip_a(case, show_y_title = (case == GRID_COL_ORDER[1]))),
  list(make_divider()),
  list(make_block_banner("Weight type varies, conditioned on clustering method")),
  lapply(GRID_COL_ORDER, function(case) make_strip_b(case, show_y_title = (case == GRID_COL_ORDER[1])))
)

grid_design <- c(
  area(t = 1, l = 1, b = 1, r = 4),
  area(t = 2, l = 1, b = 2, r = 1), area(t = 2, l = 2, b = 2, r = 2),
  area(t = 2, l = 3, b = 2, r = 3), area(t = 2, l = 4, b = 2, r = 4),
  area(t = 3, l = 1, b = 3, r = 4),
  area(t = 4, l = 1, b = 4, r = 4),
  area(t = 5, l = 1, b = 5, r = 1), area(t = 5, l = 2, b = 5, r = 2),
  area(t = 5, l = 3, b = 5, r = 3), area(t = 5, l = 4, b = 5, r = 4)
)

grid_combined <- wrap_plots(grid_elements, design = grid_design) +
  plot_layout(guides = "collect", heights = c(0.35, 4, 0.15, 0.35, 7)) &
  theme(
    legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(),
    legend.text = element_text(size = 13), legend.key.size = unit(0.6, "cm"),
    legend.spacing.x = unit(4, "pt")
  ) &
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  )

ggsave(figure_path("central_regret_overview.pdf"), grid_combined,
       width = 9, height = 21.5, device = cairo_pdf, limitsize = FALSE)
ggsave(figure_path("central_regret_overview_preview.png"), grid_combined,
       width = 9, height = 21.5, dpi = 150, limitsize = FALSE)

# ---------------------------------------------------------------------------------------------
# Pareto-front plot: regret vs. total time, one combined 4-panel figure (one panel per case study,
# 2x2), across all 112 clustering x weight x n_rep combinations per panel. Color = clustering
# method, shape = weight type (same encodings as above). The Pareto frontier is already flagged
# per-row in regret_summary.csv's `pareto_frontier` column (no other combination has both lower
# regret AND lower time); Pareto-optimal points are highlighted (colored/shaped, with error bars)
# but NOT connected by a line -- with only a handful of Pareto-optimal points per panel, a
# connecting line reads as a trend the data doesn't actually support (regret and time don't trade
# off smoothly along it), so the points are left to speak for themselves.
# Stats bars on both axes: vertical = +-1 SD regret over 5 seeds (regret_summary.csv); horizontal
# = +-1 SD wall-clock time over the same 5 seeds (computed here from regret_by_seed.csv, since
# time is not deterministic even for methods whose regret is -- floating point/system noise).
by_seed <- read.csv(data_path("regret_by_seed.csv")) %>%
  filter(normalization == "economic") %>%
  mutate(n_rep_periods = as.numeric(n_rep_periods))

time_stats <- by_seed %>%
  group_by(case_study, clustering_type, weight_type, n_rep_periods) %>%
  summarise(time_sd_s = sd(total_time_s), .groups = "drop")

# Non-frontier combinations are drawn as small, uncolored background context (no error bars --
# with ~100 combinations per case study, bars on every point made the plot unreadable); only the
# Pareto-optimal combinations are colored/shaped by method, carry error bars, and are connected by
# the frontier line, so the reader's eye goes to the trade-off that actually matters.
make_pareto_plot <- function(case) {
  d <- df %>%
    filter(case_study == case) %>%
    left_join(time_stats, by = c("case_study", "clustering_type", "weight_type", "n_rep_periods"))
  d$clustering_label <- factor(clustering_names[d$clustering_type], levels = clustering_names)
  d$weight_label <- factor(weight_names[d$weight_type], levels = weight_names)
  d$xmin <- pmax(d$time_mean_s - d$time_sd_s, 1e-3)
  d$xmax <- d$time_mean_s + d$time_sd_s
  d$ymin <- pmax(d$regret_mean_pct - d$regret_sd_pct, 0)
  d$ymax <- d$regret_mean_pct + d$regret_sd_pct

  pareto_pts <- d %>% filter(pareto_frontier) %>% arrange(time_mean_s)
  other_pts <- d %>% filter(!pareto_frontier)
  y_cap <- min(REGRET_CAP, max(d$regret_mean_pct))

  # Each panel gets its OWN legend, restricted to just the clustering/weight combinations that
  # actually appear on THIS case study's Pareto front (not the full 7/4-category list) -- a
  # shared/collected legend across all four panels would otherwise show categories that never
  # appear as Pareto-optimal here, and a global frontier-driven legend doesn't generalize across
  # panels with different Pareto-optimal methods anyway.
  clust_breaks <- clustering_names[names(clustering_names) %in% unique(as.character(pareto_pts$clustering_type))]
  weight_breaks <- weight_names[names(weight_names) %in% unique(as.character(pareto_pts$weight_type))]

  ggplot() +
    geom_point(data = other_pts,
               aes(x = time_mean_s, y = regret_mean_pct, color = clustering_label, shape = weight_label),
               size = 1, alpha = 0.35, stroke = 0.6) +
    geom_segment(data = pareto_pts,
                 aes(x = xmin, xend = xmax, y = regret_mean_pct, yend = regret_mean_pct),
                 color = "grey40", alpha = 0.6, linewidth = 0.4) +
    geom_errorbar(data = pareto_pts, aes(x = time_mean_s, ymin = ymin, ymax = ymax),
                  color = "grey40", width = 0, alpha = 0.6, linewidth = 0.4) +
    geom_point(data = pareto_pts,
               aes(x = time_mean_s, y = regret_mean_pct, color = clustering_label, shape = weight_label),
               size = 3.2, stroke = 1.2) +
    coord_cartesian(ylim = c(0, y_cap)) +
    scale_x_log10(name = "Total time [s]", labels = y_labels) +
    scale_y_continuous(name = "Regret [%]", labels = y_labels) +
    scale_color_manual(name = "Clustering (Pareto front)", values = okabe_ito,
                        breaks = unname(clust_breaks), limits = unname(clustering_names)) +
    scale_shape_manual(name = "Weight (Pareto front)", values = weight_shapes,
                        breaks = unname(weight_breaks), limits = unname(weight_names)) +
    labs(title = case_titles[[case]]) +
    guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1, stroke = 1.2)),
           shape = guide_legend(override.aes = list(size = 3.2, alpha = 1, stroke = 1.2))) +
    base_theme +
    theme(legend.position = "right", legend.box = "vertical")
}

# One combined 2x2 figure (case_titles' order: GEP, 5-bus, P2X, 118-bus); each panel keeps its OWN
# legend (see make_pareto_plot) rather than a shared/collected one, since the Pareto-optimal
# clustering/weight combinations differ across case studies.
pareto_panels <- lapply(names(case_titles), make_pareto_plot)
pareto_combined <- wrap_plots(pareto_panels, ncol = 2)

ggsave(figure_path("central_pareto_overview.pdf"), pareto_combined,
       width = 11, height = 9, device = cairo_pdf)
ggsave(figure_path("central_pareto_overview_preview.png"), pareto_combined,
       width = 11, height = 9, dpi = 150)

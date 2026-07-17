# Shared helpers for the paper's R plotting scripts (regret.R, runtime_breakdown.R). Fully
# self-contained within THIS repo: they read CSVs from, and write figures into,
# analysis/output/{data,figures}/ here -- never the sibling `clustering` repo directly. Copying the
# finished data/figures into `clustering` is a separate, deliberate step the user does by hand.
#
# Path resolution is anchored to the repo root, found by walking UP from the current working
# directory looking for Project.toml (this repo's own marker file) -- NOT by trying to introspect
# the running script's own file path. The latter (`commandArgs()`'s `--file=`, or scanning
# `sys.frames()` for an `ofile`) is exactly what broke under at least one runner that evaluates the
# script without a normal `Rscript`/`source()` call frame (some editor "Run file" integrations do
# this), throwing "Could not determine this script's own file path". Walking up from `getwd()`
# has no such runner-dependence: it only requires the R process's working directory to be
# somewhere at or below the repo root, true for every normal way of launching R here.

find_repo_root <- function(start = getwd(), marker = "Project.toml") {
  dir <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, marker))) return(dir)
    parent <- dirname(dir)
    if (parent == dir) {
      stop(
        "Could not find the repo root (looked for '", marker, "' walking up from ", start, "). ",
        "Run this script with the R process's working directory inside the experiments repo."
      )
    }
    dir <- parent
  }
}

.repo_root <- find_repo_root()
.output_dir <- file.path(.repo_root, "analysis", "output")
dir.create(file.path(.output_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(.output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

data_path <- function(...) file.path(.output_dir, "data", ...)
figure_path <- function(...) file.path(.output_dir, "figures", ...)

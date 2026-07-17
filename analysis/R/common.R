# Shared helpers for the paper's R plotting scripts (regret.R, runtime_breakdown.R). These scripts
# read CSVs from, and write figures into, the sibling `clustering` repo, so path resolution must be
# anchored to THIS file's own location, not the caller's working directory: a bare relative path
# like "../../clustering/data/x.csv" only resolves correctly when Rscript is launched from this
# exact directory, and silently breaks (wrong file, or "cannot open the connection") when run from
# analysis/ or analysis/R/ instead. `this_file()` recovers the running script's own path from
# Rscript's `--file=` argument (falling back to the `source()` call frame, for interactive/sourced
# use), so `data_path()`/`figure_path()` below work regardless of the invoker's cwd.

this_file <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1])))
  }
  # Rscript'd via source(): the enclosing frame knows its own file.
  frame <- sys.frames()[[1]]
  normalizePath(frame$ofile)
}

# analysis/R/ -> repo root (Tulipa/) is two levels up; clustering/ is a sibling of experiments/.
.repo_root <- normalizePath(file.path(dirname(this_file()), "..", ".."))
.clustering_root <- normalizePath(file.path(.repo_root, "..", "clustering"), mustWork = FALSE)

data_path <- function(...) file.path(.clustering_root, "data", ...)
figure_path <- function(...) file.path(.clustering_root, "figures", ...)

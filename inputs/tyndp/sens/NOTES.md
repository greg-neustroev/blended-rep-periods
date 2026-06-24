# PGD sensitivity sweep

Same data as `p2x` (every file is a symlink to `../p2x/`); a focused sweep that
varies the projected-gradient-descent tolerance `tol` (1e-1 … 1e-5) and the weight
initialization `init` (`auto` / `default` / `moore_penrose`) on hull/convex at 20
and 40 RPs. Purpose: justify the headline `tol = 1e-2` (fit time vs projection
error vs RP-selection stability vs regret) and show the initialization is not
load-bearing. P2X is the weight-fitting-dominated regime where `tol` has the most
effect; the recorded `kappa` per run confirms it is a high-condition-number case.

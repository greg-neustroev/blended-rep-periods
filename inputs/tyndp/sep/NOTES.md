# P2X Storage Expansion Planning (SEP) variant

The TYNDP P2X dispatch case with the short-term **batteries** (technology
`Battery`, `is_seasonal = false`) made investable, turning the dispatch problem
into a battery storage-expansion problem. Every file except `investments.csv` is
a symlink to `../p2x/`; only the investability differs.

- **Investable assets:** all 25 `*_Battery` assets.
- **Capital cost:** a stylised annualised value of `50.0` in the dataset's
  investment-cost units (where gas ≈ 23 and offshore wind ≈ 88), from scaling a
  real ≈ €90k/MW-yr battery CAPEX into that range. **PLACEHOLDER — review and
  re-tune after the first run** if investment is degenerate (all-zero or huge).
- Investing scales both battery **power** (charge/discharge limit) and **energy**
  (SoC cap) together at the battery's fixed duration.
- **Regret metric:** `investment_regret` (fix invested battery units, re-solve full).

# Replication Package for "Trade, Misallocation, and Capital Market Integration"

**Author:** Laszlo Tetenyi  
**Journal:** Journal of the European Economic Association  

---

## Data Availability Statement

The replication materials draw on both restricted administrative microdata and publicly available macroeconomic datasets.

**Restricted data.** The firm-level balance sheet data used in the empirical analysis originate from the tax return records of Hungary's National Tax and Customs Administration (NAV, formerly APEH), harmonized by the ELTE KRTK Databank (Hungarian Academy of Sciences, Centre for Economic and Regional Studies). Access is governed by data-sharing agreements with the data owner and can be obtained through the ELTE KRTK Databank upon application. Researchers wishing to replicate the firm-level empirical results should contact the ELTE KRTK Databank directly regarding access conditions. The restricted file is not included in this replication package. An exemption for restricted data was granted by the co-editor, and the dataset has been shared with the Data Editor temporarily for replication purposes only.

**Publicly available data.** All remaining data used in this paper are publicly accessible. See the Data Sources section below for individual citations and access links.

**Replication code.** The Julia code required to reproduce all model-based figures and tables is available in the online appendix at https://github.com/lt1245/trade_misallocation_capital_appendix. Stata do-files for the empirical analysis are included in this package; their execution requires access to the restricted ELTE KRTK firm-level dataset described above.

---

## Data Sources and Citations

| Dataset | Use in paper | Access |
|---|---|---|
| ELTE KRTK Databank — NAV/APEH firm-level tax records, Hungary | Firm-level empirical analysis (Sections 2–3) | Upon application to ELTE KRTK Databank |
| World Bank World Development Indicators (WDI) | Macro aggregates: GDP, trade, credit | https://databank.worldbank.org/source/world-development-indicators |
| World Bank World Integrated Trade Solution (WITS) | Import tariff data | https://wits.worldbank.org |
| Bank for International Settlements — Total Credit Statistics | Credit-to-GDP ratio | https://www.bis.org/statistics/totcredit.htm |
| OECD Trade in Value Added (TiVA) | Export value-added decomposition | https://www.oecd.org/en/data/datasets/trade-in-value-added.html |
| Penn World Table, version 7.1 (Heston, Summers, and Aten, 2012) | Calibration: consumption, investment shares | https://pwt.sas.upenn.edu/ |
| Penn World Table, version 9.0 (Feenstra, Inklaar, and Timmer, 2015) | Cross-country TFP panel (Figures 15–17) | https://www.rug.nl/ggdc/productivity/pwt/ |
| Chinn and Ito (2006) — Capital account openness index (KAOPEN) | Capital account liberalization dates | http://web.pdx.edu/~ito/Chinn-Ito_website.htm |
| World Inequality Database (WID) | Income share series | https://wid.world |
| Hungarian State Budget Execution Reports (*Magyar Állami Zárszámadás*) | Import tariff revenues | https://net.jogtar.hu |

**References:**  
Chinn, M. D., and H. Ito (2006). "What Matters for Financial Development? Capital Controls, Institutions, and Interactions." *Journal of Development Economics* 81(1): 163–192.  
Feenstra, R. C., R. Inklaar, and M. P. Timmer (2015). "The Next Generation of the Penn World Table." *American Economic Review* 105(10): 3150–3182.  
Heston, A., R. Summers, and B. Aten (2012). *Penn World Table Version 7.1.* Center for International Comparisons of Production, Income and Prices, University of Pennsylvania.

---

## Package Contents

```
replication_package/
├── Data/                          # Public datasets (see Data Sources above)
│   ├── complete_country_panel.dta      # Country-level panel, EU+EFTA, 1990–2008
│   ├── complete_country_panel_map.dta  # Same, merged with shapefile IDs for spmap
│   ├── Europe_coord.dta                # European shapefile polygon vertices
│   ├── Hungary_macrodata.xlsx          # Hungarian aggregate fiscal and macro data
│   └── importshare_WB.dta              # Import shares from World Bank WDI
├── Empirics/                      # Stata replication code
│   ├── empirics_replication.do         # Master do-file: all empirical figures and tables
│   ├── figures_15_18_table13.do        # Choropleth maps and cross-country regression
│   └── olleyp.ado                      # Custom Olley–Pakes productivity decomposition ado
├── Figures/                       # All figures in the paper (PDF, EPS, TikZ)
├── Model/                         # Julia model code
│   ├── code_new.jl                     # Core model: household and firm problems, equilibrium
│   ├── CompEcon.jl                     # Computational economics utilities (basis matrices etc.)
│   ├── main_script.jl                  # Solves the four main steady-state equilibria
│   ├── transition.jl                   # Solves transition path equilibria
│   ├── figures_tables.jl               # Generates all model-based figures and tables
│   ├── developed_steady_states.jl      # Steady states for developed-economy comparison
│   ├── fixed_cost_steady_states.jl     # Steady states for fixed-cost extension
│   ├── search_developed.jl             # Equilibrium search for developed-economy calibration
│   ├── Project.toml                    # Julia environment: direct package dependencies
│   └── Manifest.toml                   # Julia environment: pinned versions of all dependencies
├── Tables/                        # LaTeX table source files (generated by code)
├── Transition_files/              # Pre-computed transition price vectors (21 CSV files)
├── data_availability.pdf          # Data availability statement
├── LICENSE                        # CC BY 4.0 license
└── README.md                      # This file
```

> **Note:** The restricted firm-level dataset is not included in this package. It must be obtained from ELTE KRTK Databank and placed at the path specified in `Empirics/empirics_replication.do`.

---

## Software Requirements

### Julia (model)
- **Version:** [Julia 1.12, although the core code should work on any stable (1.0) Julia release]
- **Packages:**

| Package | Purpose |
|---|---|
| `DelimitedFiles` | Read/write CSV and delimited files |
| `Optim` | Numerical optimization |
| `DataFrames` | Tabular data manipulation |
| `CSV` | CSV I/O |
| `SharedArrays`, `Distributed` | Parallel computing |
| `LinearAlgebra`, `Statistics` | Standard numerical routines |
| `SparseArrays`, `StructArrays` | Sparse matrices and structured arrays |
| `BasisMatrices` | Chebyshev and spline basis functions |
| `Arpack` | Sparse eigenvalue decomposition |
| `LeastSquaresOptim` | Nonlinear least-squares optimization |
| `BenchmarkTools` | Timing utilities |
| `QuantEcon` | Grid generation and quadrature routines |
| `Plots`, `StatsPlots` | Figure generation |
| `StatsBase` | Summary statistics |
| `GeometryTypes` | Geometric primitives (used in plotting) |
| `GLM` | Linear regression (for table outputs) |
| `Printf` | Formatted string output |
| `XLSX` | Excel file I/O |

The Julia environment is fully specified by `Model/Project.toml` and `Model/Manifest.toml`. To restore the exact package versions, run `] activate Model/; instantiate` from the Julia REPL before executing any model code.

### Stata (empirics)
- **Version:** [Stata 19.5]
- **User-written packages** (installed automatically by `empirics_replication.do`):
  `matsort`, `gtools`, `egenmore`, `require`, `reghdfe`, `eventdd`
- **User-written packages** (installed automatically by `figures_15_18_table13.do`):
  `spmap`, `estout`
- **Custom ado:** `olleyp.ado` (included in `Empirics/`; copied to personal ado directory by the master do-file)

---

## Replication Instructions

### Step 1 — Empirical analysis (Stata)

1. Obtain access to "NAV-merleg" from ELTE KRTK Databank and note its path on your machine.
2. Open `Empirics/empirics_replication.do` and update the following paths at the top of the file to match your local environment:
   - Line 15: `cd` working directory
   - Line 18: `use` path to `firm_data_appended.dta`
   - Lines 20–21: `local figures` and `local tables` output directories
3. Run `Empirics/empirics_replication.do`. This produces all empirical figures (Figures 1–14) and tables in the `Figures/` and `Tables/` directories. It also creates a txt file `Tables/Hungary_firm_level_data.txt` that is used as an input for `Model/figures_tables.jl` to create tables containing mixed empirical and model variables.
4. Run `Empirics/figures_15_18_table13.do`. Update the data and output paths at the top of the file if needed. This produces Figures 15–18 and Table 13 using only public data from the `Data/` folder.

### Step 2 — Model computation (Julia)

All model outputs depend on pre-computed steady-state solutions and transition price vectors. The `Transition_files/` folder contains pre-computed transition paths so that `figures_tables.jl` can be run without solving for the full transition. Solving for the transition path requires at least 16GB RAM, preferably 32GB RAM. The code assumes at least 4 CPU cores (parallel execution required).

**To reproduce figures and tables only (fast path):**

```julia
include("Model/figures_tables.jl")
```

This reads the pre-computed solutions from `Transition_files/` and generates all model-based figures and tables in `Figures/` and `Tables/`. This is the only file that needs to be run to reproduce the results pertaining to the model, it calls all necessary data files and solutions.

**To load steady states for only the main text from scratch:**

```julia
include("Model/main_script.jl")
```

**To recompute transition paths from scratch (computationally intensive):**

```julia
# In transition.jl, set case_final and load_solution = 0
include("Model/transition.jl")
```

See comments at the top of `transition.jl` for the list of cases and their correspondence to the scenarios in the paper (some are not utilized in the paper).

### Step 3 — Compile the paper

The paper manuscript is not included in this replication package. After running Steps 1 and 2, all figures and tables referenced via `\input{Tables/...}` and `\includegraphics{Figures/...}` will be present in the `Tables/` and `Figures/` directories and ready for inclusion in the manuscript.

---

## Expected Running Times

Times were measured on an Intel Core Ultra 7 265U (12 cores, 32 GB RAM) running Windows 11 Enterprise 64-bit.

| Task | Approximate time (minutes) |
|---|---|
| `empirics_replication.do` (with restricted data) | 10 |
| `figures_15_18_table13.do` | 2 |
| `figures_tables.jl` (from pre-computed files) | 30 |
| `main_script.jl` (steady states only) | 5 |
| `transition.jl` (full transition, all cases) | 15 |

---

## License

The replication materials in this package — code, data, and generated tables — are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). See the `LICENSE` file for details.

The paper manuscript is not included in this package and remains subject to the copyright terms of the Journal of the European Economic Association.

---

## Notes

- The `Transition_files/` CSV files store equilibrium price vectors for the transition paths. Their filenames encode the scenario (e.g., `vec_price_clCM.csv` = closed capital markets).
- All LaTeX table files in `Tables/` are generated programmatically by the replication code; do not edit them by hand.

# Bayesian latent growth curve modeling of OGTT glucose and insulin: replication code

Code to reproduce the R-generated figures in

> Waldman, M. R., Hagen-Lillevik, S., Bosma, G., Schmiege, S. J., & Hernandez, T. L.
> *Bayesian Latent Growth Curve Modeling for Repeated Blood Biomarker Measurements.*
> Submitted to *Nurse Researcher*.

The manuscript fits a Bayesian latent growth curve model (B-LGCM) with a
Legendre-polynomial growth basis to oral glucose tolerance test (OGTT) glucose
and insulin trajectories from the CHOICE randomized trial (two visits, 31 and
36 weeks gestation), separates measurement error from true-score variation,
and compares measurement-error-corrected treatment effects on area under the
curve (AUC) with the conventional observed-data estimate.

## What is and is not here

Included:

- `R/` — every function the figures depend on (data preparation, Legendre
  basis, brms formula and prior construction, AUC treatment-effect extraction,
  trajectory summaries, plotting).
- `priors/` — the digitized literature OGTT curves (Carlson et al., 2025 for
  glucose; Mittendorfer et al., 2013 for insulin) from which the
  Legendre-coefficient priors are derived at fit time, with the audit scripts
  and READMEs documenting the derivation.
- `scripts/` — the numbered pipeline and a verification script.

Not included, by design:

- The participant-level OGTT data (`data/raw/`, see its README for the
  expected layout).
- The fitted models (`fits/`), because a `brmsfit` stores the model frame.
- Figure 3, the path diagram, which was drawn outside R.

Without the data the pipeline cannot execute end to end. The code is
published so that the model specification, priors, estimands and figure
construction are fully inspectable, and so that the analysis can be re-run by
anyone with authorized access to the data.

## Figures produced

| Figure | Script | Needs |
|---|---|---|
| Fig 1: true-score vs observed trajectory (conceptual) | `04_make_figures.R` | nothing |
| Fig 2: prior, likelihood, posterior for residual SD (conceptual) | `04_make_figures.R` | nothing |
| Fig 4: 36-week OGTT trajectories by arm with AUC treatment-effect insets | `04_make_figures.R` | primary fit |
| Fig 5: prior-sensitivity forest of the adjusted AUC treatment effect | `04_make_figures.R` | all five fits |
| Fig 6: priorsense prior power-scaling robustness | `05_make_figure6_priorsense.R` | primary fit |

## Running

From the repository root, with the workbook location and (optionally) an
existing fits directory supplied through environment variables:

```
# PowerShell
$env:BLGCM_DATA_XLSX = "C:\path\to\OGTTdata.xlsx"
$env:BLGCM_FITS_DIR  = "C:\path\to\fits"      # optional; default fits/
Rscript run_all.R
```

or step by step:

```
Rscript scripts/01_prepare_data.R            # -> derived/joined_01_long.rds
Rscript scripts/02_fit_models.R              # -> fits/*.rds  (hours; skipped if present)
Rscript scripts/03_derive_estimands.R        # -> derived/estimands.rds
Rscript scripts/04_make_figures.R            # -> figures/Fig1, Fig2, Fig4, Fig5
Rscript scripts/05_make_figure6_priorsense.R # -> figures/Fig6
Rscript scripts/verify_model_spec.R          # fits match the rebuilt specification
```

`verify_model_spec.R` rebuilds the formula, priors and data for each of the
five fits and checks that the generated Stan program and Stan data are
identical to those stored in the fit objects. This is how we confirmed that
the code in this repository is the code that produced the reported fits.

### Model fitting settings (`scripts/02_fit_models.R`)

Four chains, 1,000 warmup and 2,000 sampling iterations each, `adapt_delta =
0.99`, `max_treedepth = 25`, seed 42, cmdstanr backend with within-chain
threading. The sensitivity fits scale every prior SD by 0.1, 10 or 100 and,
for the independence model, drop the cross-response random-effect
correlations; their chains are initialized from the primary posterior.

## Software versions used for the reported results

| Software | Version |
|---|---|
| R | 4.5.1 |
| brms | 2.23.0 |
| cmdstanr | 0.9.0.9000 |
| CmdStan | 2.37.0 |
| rstan | 2.36.0.9000 |
| posterior | 1.7.0 |
| priorsense | 1.2.0 |
| bayestestR | 0.17.0 |
| mice | 3.18.0 |
| ggplot2 | 4.0.1 |
| patchwork | 1.3.2 |
| ggtext | 0.1.2 |
| orthopolynom | 1.0.6.1 |
| readxl | 1.4.5 |
| tidyverse | 2.0.0 |
| pbapply | 1.7.4 |

Package attach lists are at the top of each script.

## Layout

```
R/                      helper functions (sourced by R/_setup.R)
priors/                 literature OGTT reference curves + prior-derivation audit scripts
scripts/01..05_*.R      pipeline
scripts/verify_model_spec.R
run_all.R               runs scripts 01-05 in order
data/raw/               (empty) expected workbook layout in README
fits/                   (empty) fitted models land here
derived/                (empty) intermediate objects
figures/                (empty) figure outputs
```

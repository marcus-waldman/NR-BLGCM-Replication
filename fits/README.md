# Fitted models (not distributed)

`scripts/02_fit_models.R` writes five `brmsfit` objects here. They are excluded
by `.gitignore` because a `brmsfit` stores the model frame (`fit$data`), which
holds participant-level rows.

| Key | File | Specification |
|---|---|---|
| `primary` | `fit_joined_01_v2.rds` | Informative priors (prior SD scale 1) |
| `tenthx` | `fit_joined_tenthx_01.rds` | Prior SDs scaled 0.1x |
| `diffuse_10x` | `fit_diffuse_10x_01.rds` | Prior SDs scaled 10x |
| `diffuse_100x` | `fit_diffuse_100x_01.rds` | Prior SDs scaled 100x |
| `indep_10x` | `fit_indep_10x_01.rds` | Independence model (no cross-response random-effect correlations), prior SDs scaled 10x |

If you already have the fits elsewhere, set `BLGCM_FITS_DIR` to that directory
and step 2 will skip fitting. `scripts/verify_model_spec.R` checks that fits in
`BLGCM_FITS_DIR` match the specification rebuilt by this repository.

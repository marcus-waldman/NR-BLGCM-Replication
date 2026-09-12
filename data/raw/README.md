# Participant data (not distributed)

The OGTT data come from the CHOICE randomized controlled trial and contain
participant-level measurements. They are not part of this repository and are
excluded by `.gitignore`.

To run the pipeline, either place the workbook at `data/raw/OGTTdata.xlsx` or
point the environment variable `BLGCM_DATA_XLSX` at its location. The reader
(`R/read_ogtt_data.R`, `R/prepare_ogtt_long.R`) expects:

| Sheet | Visit | Suffix in the model |
|---|---|---|
| 3 | Baseline, 31 weeks gestation | `0` |
| 4 | Post-intervention, 36 weeks gestation | `1` |

Columns used on each sheet:

| Column | Meaning |
|---|---|
| `PTID` | Participant id |
| `Group` | Arm: `0` = LC/Conventional, `1` = CHOICE |
| `Completers` | Non-missing for study completers (not used: all rows on a sheet are read) |
| `Glucose-15`, `Fasting Glucose`, `Glucose30`, `Glucose60`, `Glucose90`, `Glucose120` | Plasma glucose (mg/dL) at -15, 0, 30, 60, 90, 120 min |
| `Insulin -15`, `Fasting Insulin`, `Insulin30`, `Insulin60`, `Insulin90`, `Insulin120` | Plasma insulin at the same times |

The analytic sample is every participant with at least one 36-week value
(N = 46); `scripts/01_prepare_data.R` stops if it finds a different count.

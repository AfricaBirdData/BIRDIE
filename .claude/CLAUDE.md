# BIRDIE Pipeline – Session Context

## Project
- **Repo:** `2026/BIRDIE` on Desktop — this is the active working branch (branch of BIRDIE-main on GitHub)
- **R package:** BIRDIE is an R package. Source lives in `2026/BIRDIE/R/`. Pipeline scripts live in `analysis/scripts/`.
- **Do NOT edit** `BIRDIE-main` or `BIRDIE_OLD` — all edits go to `2026/BIRDIE/`

## Active Pipeline
- **Script being tested:** `analysis/scripts/pipeline_script_dst.R`
- **Module:** `dst` (distribution indicators, occupancy modelling)
- **Config:** `year = 2024`, `dur = 3`, `package = "spOccupancy"`, `server = FALSE`
- **Output dir:** `analysis/output/` (relative to project root)

## Bug Fixed This Session
- **Error:** `cannot open file 'analysis/output/4/occu_det_dat_sa_22_24.csv': No such file or directory`
- **Cause:** `ppl_run_pipe_dst1.R` tried to write to `analysis/output/<sp_code>/` but that subdirectory didn't exist
- **Fix applied:** Added `dir.create(file.path(config$out_dir, sp_code), recursive = TRUE, showWarnings = FALSE)` in `R/ppl_run_pipe_dst1.R` before the first write to the species subdirectory
- **Fix already present** in `2026/BIRDIE/R/ppl_run_pipe_dst1.R:59` (comment: `# create directory?`)

## Package Load Issue
- `library(BIRDIE)` threw: `lazy-load database ... is corrupt`
- **Cause:** Previously installed version of BIRDIE in the R library is corrupt
- **Fix:** Run in R console:
  ```r
  remove.packages("BIRDIE")
  install.packages("~/Desktop/2026/BIRDIE", repos = NULL, type = "source")
  ```
- **Alternative (no install needed):** `devtools::load_all("~/Desktop/2026/BIRDIE")`
- User note: previously ran pipeline without reinstalling — likely used RStudio Build tab or had a working install already

## Key File Locations
| File | Path |
|------|------|
| Pipeline script | `analysis/scripts/pipeline_script_dst.R` |
| Main pipeline function | `R/ppl_run_pipe_dst1.R` |
| Site/visit data prep | `R/ppl_create_site_visit.R` |
| Config function | `R/configPipeline.R` |
| Utility functions | `R/utils.R` |
| Output directory | `analysis/output/` |
| Log files | `analysis/output/reports/pipe_log_dst_*.csv` |

## Pipeline Flow (dst module)
1. `pipeline_script_dst.R` loops over all species in `config$species`
2. Calls `ppl_run_pipe_dst1()` for each species × year
3. Inside `ppl_run_pipe_dst1`:
   - **data:** calls `ppl_create_site_visit()` → downloads ABAP/GEE data → writes site/visit/det CSVs
   - **fit:** calls `ppl_fit_occu_model()` → saves `.rds`
   - **diagnose:** calls `ppl_diagnose_occu()` → saves diagnostics CSV + PPC rds
   - **summary:** calls `ppl_summarise_occu()` or `createPredFromAbap()` → saves predictions CSV + plot
4. Species-specific files go to `analysis/output/<sp_code>/`
5. Shared files (GEE site/visit data) go to `analysis/output/`

## Notes
- `force_site_visit = TRUE` is set — always re-prepares site/visit data
- `force_abap_dwld = FALSE` — uses cached ABAP download from `tempdir()` if available
- GEE downloads only triggered if files missing or `force_gee_dwld = TRUE`

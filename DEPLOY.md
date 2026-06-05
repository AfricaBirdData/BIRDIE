# Deploying / Updating BIRDIE on the Server

**Two rules that cause almost every server problem:**

1. The pipeline runs the *installed* BIRDIE package (`library(BIRDIE)` + the `furrr`
   workers load `packages = "BIRDIE"`), **not** the source files in the git checkout.
   `git pull` updates the source only — **you must reinstall after every pull.**
2. **Never use `sudo`.** BIRDIE and all its dependencies live in the `birdie` user
   library (`/home/birdie/R/x86_64-pc-linux-gnu-library/4.3`). `sudo` runs as root,
   whose `.libPaths()` does **not** include that library, so under `sudo` BIRDIE and
   its deps are invisible (`there is no package called 'BIRDIE'`, or
   `dependencies ... are not available`). Run everything as the `birdie` user.

```
git pull  →  reinstall (no sudo)  →  run (no sudo)
```

## Steps

```bash
# 1. Get the latest code (server tracks origin/main)
cd /drv_birdie/Working/git/BIRDIE
git pull origin main

# 2. Reinstall the package — NO sudo, skip dependency reinstall (deps already exist)
Rscript -e 'install.packages("/drv_birdie/Working/git/BIRDIE", repos = NULL, type = "source", dependencies = FALSE)'

# 3. Verify the installed package is the new one (check a known fix)
Rscript -e 'library(BIRDIE); cat("loaded from:", find.package("BIRDIE"), "\n"); print(args(BIRDIE:::fitSpOccu))'

# 4. Run the pipeline — NO sudo
Rscript /drv_birdie/Working/git/BIRDIE/analysis/scripts/pipeline_script_dst_parall.R
```

## Output directory permissions

The pipeline writes to `config$out_dir` = **`/drv_birdie/birdie_ftp`** (logs,
per-species results, etc.). If a *previous* `sudo` run created files there as root,
a no-sudo run fails with:

```
cannot open file '/drv_birdie/birdie_ftp/reports/pipe_log_...csv': Permission denied
```

Fix ownership once (do **not** switch back to running R with sudo):

```bash
sudo chown -R birdie:birdie /drv_birdie/birdie_ftp
```

(`sudo` is fine here — it's a one-off `chown`, not running R.)

## Common pitfalls

- **Forgot to reinstall after `git pull`** → pipeline still runs the old code and a
  fixed bug reappears. Always reinstall.
- **`there is no package called 'BIRDIE'` / `dependencies ... are not available`** →
  you used `sudo`. Drop it. (A *failed* `sudo R CMD INSTALL` also deletes the
  existing BIRDIE — another reason to avoid it.)
- **`Permission denied` writing outputs** → `chown` the output dir as above, then
  rerun without sudo.
- **A dependency genuinely needs (re)installing** → install it as the `birdie` user,
  no sudo: `Rscript -e 'install.packages("<pkg>")'`, then reinstall BIRDIE.

## Local development (laptop)

No reinstall needed each time — use `devtools::load_all("~/Desktop/2026/BIRDIE")`
to run against the live source, or reinstall with
`install.packages("~/Desktop/2026/BIRDIE", repos = NULL, type = "source")`.

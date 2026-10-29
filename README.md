# NHATS digital disconnection — analysis code

Version: v1.0.0 (29 October 2026)

Study: Complete Digital Disconnection and Mortality or Residential-Care Transition in Older Adults With Cardiovascular Risk: An Adjusted Association Study in NHATS

Study authors: Peng Liu, Ruilong Gao, Yaliu Yang, Zhengqin Zhai, and Yifeng Zhou.
Submission: 64acdf5b-c7ad-4902-94ec-b9fa55000aed
Repository: https://github.com/aimar1215/nhats-digital-disconnection

This repository contains analysis code, documentation, and aggregate results. NHATS participant-level data are not redistributed. No DOI has been assigned.

## Inputs
No participant-level data are included. Obtain NHATS public-use files under their applicable terms from https://www.nhats.org/data-access/nhats . Keep the following filenames together in a separate raw-data directory:
  NHATS_R12_Int_Inc_Imp_File_V2.dta
  NHATS_R13_Int_Inc_Imp_File.dta
  NHATS_Round_10B_SP_File.dta
  NHATS_Round_11_Int_Inc_Imp_File.dta
  NHATS_Round_11_SP_File_V2.dta
  NHATS_Round_12_SP_File.dta
  NHATS_Round_13_SP_File.dta
  NHATS_Round_14_SP_File.dta
  NHATS_Round_14_Tracker_File.dta

## Run

PowerShell example; replace paths:

```powershell
  $env:NHATS_RAW_DIR = 'D:/NHATS_raw'
  $env:RSCRIPT = 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe'
  $env:MI_WORKERS = '4'
  python run_analysis.py
```

Python uses only its standard library. R packages: haven, dplyr, tidyr, survey, survival, mice. This run used Windows and Times New Roman for plots. See documentation/R_environment_verified.txt for versions.
Start from a new copy of this package without generated data/baseline/tables/logs. Default MI_WORKERS=4; supported 1–8, with 1 running sequentially. Imputations use fixed seeds 20261028 + imputation number. Parallel chunks change execution order only. Cached imputation objects may be reused only for the same inputs and code. Do not mix chunk files from different worker counts or runs. For interruption recovery with unchanged inputs/settings, --start-at SCRIPT_NAME resumes at a named stage; prerequisites must already exist. A full fresh run is the preferred verification route.

## Outputs And Checking
Generated tables/ holds aggregate results, figures/ holds plots. expected_aggregate_results/ contains the verified reference CSVs. documentation/reproduction_comparison.csv records comparisons to the prior revision. Numeric tolerance was 1e-9. All 33 tables matched numerically; two Figure*_data tables differed only in label capitalization and dash character. Baseline derived objects also matched the submitted analytic objects at all.equal tolerance 1e-12 including attributes.
An end-to-end raw-data run completed in a fresh output directory with new R processes on the existing local R 4.5.1 installation. It was not a clean-machine or cross-platform test. All 20 MI replicates completed. Initial replicates 1–2 were reused from this same fresh run after switching to parallel execution. Warnings concerning startup locale and packages built under R 4.5.3 did not stop the run; scripts set a Windows UTF-8 locale.

## Scope And Interpretation
Baseline reconstruction preserves submitted definitions. Revision analyses address reviewer requests: survey variance, repeated entries, M4, first-entry and interval sensitivities, MI, censoring-aware risks, IPCW, and component outcomes. They do not establish causal effects or validate a prediction model. MI uses official income draws plus chained equations for residual item nonresponse; structural categories remain. Standardized risks condition on target covariates; IPCW inference conditions on estimated weights. Fine-Gray uses discrete-round expansion and a survey Taylor exposure-contrast variance; the aliased nuisance indicator is documented in the supplement. The conceptual DAG is illustrative, not a verified sufficient adjustment set.

## Data access and reuse
Obtain the raw data directly from NHATS after registration and agreement to its [Conditions of Use](https://www.nhats.org/conditions-of-use). The official data provider's terms apply to the data independently of this repository. Generated participant-level datasets and imputation objects must remain outside the public repository.

Aggregate estimates in `MI_individual_estimates.csv` are estimates for individual IMPUTATIONS, not individual participants. The public package excludes the manuscript, reviewer correspondence, internal audit reports, and participant-level data.

No open-source license is granted by this release. Copyright and reuse rights remain with the respective rights holders; public visibility is not a grant of a separate reuse license.

# WA coexistence model workflow

This folder contains the R scripts used to analyse the WA *Trachymene* coexistence experiment. The main analyses use `fullWAdata.csv` as the central analysis dataset. The workflow estimates fecundity responses, compares competition-model structures, calculates low-density growth rates, and summarises predicted coexistence outcomes for `TRCY` (*T. cyanopetala*) and `TROR` (*T. ornata*) across sites, cover treatments, and facilitator densities.

## Main data file

- `fullWAdata.csv` — primary dataset for the coexistence, fecundity, GLMM, and plotting analyses.

Key variables used across scripts include:

- `focsp` — focal species, `TRCY` or `TROR`
- `Site` — site along the wet-to-dry gradient: Ben, GH, Nam, PJ, CS
- `Cov` — cover treatment, shade or sun
- `BLK` — block identity
- `FOCAL_GERM` — focal germination status
- `fitness` — post-germination fecundity/fitness response
- `TRCY`, `TROR` — neighbour densities used to calculate intra- and interspecific competition
- `NO_GRASS` — facilitator density, used as `fac_sc` in the coexistence models
- `dens_trt` — density/facilitation treatment label

Most model scripts filter to germinated focal individuals and analyse post-germination fecundity. Intra- and interspecific densities are assigned according to the focal species: for `TRCY`, conspecific density is `TRCY` and heterospecific density is `TROR`; for `TROR`, conspecific density is `TROR` and heterospecific density is `TRCY`.

## Germination and seed bank defaults

The primary coexistence analysis is a post-germination fitness analysis. By default, germination is **not** included in the outcome calculations:

```r
use_germination <- FALSE
```

With this default, germination is set to `g = 1` for both species, so the coexistence outcomes are driven by the competition experiment in `fullWAdata.csv`. Germination priors are loaded and validated in the model-selection scripts for bookkeeping, but they are not included in the Stage 0 or Stage 1 model likelihoods.

A germination-weighted sensitivity analysis can be run by setting:

```r
use_germination <- TRUE
```

When germination is included, the scripts draw site- and species-specific germination values from Beta priors. Direct germination evidence is used for PJ and Ben, while weakly informed priors are used for GH, Nam, and CS. The germination-prior helper script uses `collated_data_18_02_2020.csv` only to build these optional priors; it is not the main analysis dataset.

No seed bank is included in either the primary or germination-sensitivity version. Seed survival is fixed at:

```r
s = 0
```

## General approach

The analysis uses a staged workflow:

1. Build optional germination priors.
2. Prepare the `fullWAdata.csv` competition dataset.
3. Compare Beverton-Holt and Ricker growth models across alternative error structures.
4. Select the best-supported alpha/competition structure using LOO cross-validation.
5. Test whether facilitators modify competition coefficients as a sensitivity analysis.
6. Use posterior draws from the selected model to calculate lambda, alpha values, niche differences, fitness differences, low-density growth rates, and coexistence outcomes.
7. Fit complementary Bayesian GLMMs to visualise fecundity responses across site, cover, competition, and facilitator density.

Fitted Bayesian models and LOO objects are saved to disk so they can be reloaded without refitting.

## Script overview

### `00_germination_priors_CoexModels.R`

Builds optional germination priors for `TRCY` and `TROR`. Direct Beta posteriors are estimated for PJ and Ben, and weakly informed Beta priors are generated for GH, Nam, and CS. Germination is calculated as:

```r
field_number_germinants / field_number_seeds
```

Nonviable and lost seeds are retained in the denominator, zero-germination rows are dropped, and Hao Ran rows are excluded. Outputs are saved in `00_germination_priors/`.

### `00_model_selection_error_selection_CoexModels.R`

Prepares the main competition dataset from `fullWAdata.csv`, filters to germinated focal individuals, creates intra- and interspecific density variables, and compares Beverton-Holt versus Ricker models under four observation/error models: Gaussian raw fitness, Poisson count fitness, negative binomial count fitness, and Gaussian log-transformed fitness. LOO comparison outputs and the selected family/error combination are saved in `00_family_error_selection/`.

### `01_alpha_structure_selection_CoexModels.R`

Uses the selected family/error structure from Stage 0 and compares alternative alpha structures. The main comparison asks whether competition coefficients are species-only, site- and cover-varying, or cover-varying. Germination priors are loaded and carried forward, but they are not inserted into the model likelihood. Outputs are saved in `01_alpha_structure_selection/`.

### `01_alpha_facilitation_sensitivity_CoexModels.R`

Runs a targeted sensitivity analysis asking whether facilitator density modifies competition coefficients. It compares the selected baseline model, where facilitators affect lambda only, against models where facilitators also modify interspecific competition or both intra- and interspecific competition. Outputs are saved in `01b_alpha_facilitation_sensitivity/`.

### `02_plot_outcomes__CoexModels.R`

Calculates and plots coexistence outcomes from the selected Beverton-Holt model. It extracts posterior draws for lambda and alpha, optionally applies germination draws, calculates niche and fitness differences, classifies coexistence outcomes, and produces posterior support and coexistence-plane figures. By default, `use_germination <- FALSE` and `s = 0`. Outputs are saved in `02_outcomes/`.

### `02_BH_Bayesian_LDGR_Plot.R`

Generates additional Beverton-Holt posterior summaries and low-density growth rate plots. It uses the selected model to calculate lambda, alpha values, niche and fitness differences, and invasion growth rates under none versus high facilitator conditions. By default, germination is excluded and there is no seed bank.

### `03_invasion_growth_by_site_cover_fac_CoexModels.R`

Calculates invasion growth rates by site, cover, and facilitator condition. The “without facilitator” condition is evaluated at `fac_sc = 0`; the “with facilitator” condition is evaluated using the mean positive facilitator density at or above the 50th percentile within each site-by-cover combination. Outputs are saved in `03_invasion_growth/`.

### `04_GLMM_Analysis_Bayesian.R`

Fits complementary Bayesian GLMMs for fecundity using `fullWAdata.csv`. The script checks germination adequacy, builds the fecundity dataset, compares site and cover effects, tests linear and quadratic facilitator effects, runs CS inclusion/exclusion sensitivity analyses, computes moment-matched LOO comparisons, and saves predicted GLMM summaries. Outputs are saved in `env_fac_model_comparison_v1/`.

### `04_GLMM_Plots_Bayesian.R`

Creates figures from the selected GLMM, including predicted fecundity under lambda and high-competition conditions, coefficient plots, and response-surface figures showing effects of facilitator density and competitor density across sites and cover treatments.

## Recommended run order

For the full workflow, run:

```r
source("00_germination_priors_CoexModels.R")
source("00_model_selection_error_selection_CoexModels.R")
source("01_alpha_structure_selection_CoexModels.R")
source("01_alpha_facilitation_sensitivity_CoexModels.R")
source("02_plot_outcomes__CoexModels.R")
source("02_BH_Bayesian_LDGR_Plot.R")
source("03_invasion_growth_by_site_cover_fac_CoexModels.R")
source("04_GLMM_Analysis_Bayesian.R")
source("04_GLMM_Plots_Bayesian.R")
```

The GLMM scripts are complementary to the Beverton-Holt coexistence workflow and can be run separately after `fullWAdata.csv` is available.

## Main output folders

- `00_germination_priors/` — optional germination priors and summaries
- `00_family_error_selection/` — Stage 0 model fits, LOO objects, and family/error selection outputs
- `01_alpha_structure_selection/` — Stage 1 alpha-structure fits and selected model information
- `01b_alpha_facilitation_sensitivity/` — sensitivity models where facilitators modify alpha values
- `02_outcomes/` — niche/fitness summaries, posterior outcome classifications, and figures
- `03_invasion_growth/` — site-by-cover invasion growth summaries and figures
- `env_fac_model_comparison_v1/` — Bayesian GLMM fits, summaries, model comparisons, and figures

## Notes

Site order follows the wet-to-dry gradient:

```r
Ben, GH, Nam, PJ, CS
```

The core interpretation is based on posterior uncertainty rather than point estimates alone. Coexistence outcomes, niche differences, fitness differences, and invasion growth rates are all calculated from posterior draws of the selected model.

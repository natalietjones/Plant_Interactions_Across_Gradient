**WA annual plants and positive interactions across aridity gradient**

```
**The workflow has three main parts:
1) Background figures and raw-data summaries
2) Frequentist GLMM analysis using `glmmTMB`
3) Bayesian coexistence analysis**
---
The focal species are:
`TRCY` = Trachymene cyanopetala
`TROR` = Trachymene ornata
The main field dataset is:
```text
data/field/fullWAdata.csv
```
Sites are usually ordered from wetter to drier:
```text
Ben, GH, Nam, PJ, CS
```
---
Folder structure
```text
project/
├── data/
│   └── field/
├── raw_data_summaries_map/
│   └── scripts/
├── glmmTMB/
│   ├── scripts/
│   ├── glmm_outputs/
│   └── figures_final_glmm/
├── coexistence_models/
│   ├── scripts/
│   ├── 00_germination_priors/
│   ├── 00_family_error_selection/
│   ├── 01_alpha_structure_selection/
│   ├── 01b_alpha_facilitation_sensitivity/
│   └── 02_outcomes/
├── figures/
└── README.md
```
---
1. Background figures and raw-data summaries
Scripts in `raw_data_summaries_map/scripts/` generate maps, environmental summaries, germination plots, raw fecundity plots, and facilitator-density summaries.
`SM_SiteMap.R`
Creates the WA site map using 30-year April-November rainfall.

Main outputs:

```text
figures/Mean_Germination_Cover_TRT.png
figures/MeanSeedProduction_TRT.png
figures/NO_GRASS_frequency_by_site_cover.png
```

---
2. Frequentist GLMM analysis
Scripts in `glmmTMB/scripts/` fit and plot frequentist GLMMs for fecundity.
Models are run separately for `TRCY` and `TROR`. Fecundity is treated as a count response. Predictor variables include site, cover, conspecific density, heterospecific density, and facilitator density.
Facilitator density is represented by `NO_GRASS`.
`glmmTMB_analysis.R`
Fits the frequentist GLMMs and performs model selection.
The script:
reads `data/field/fullWAdata.csv`
filters to germinated focal individuals
defines fecundity from `fitness`
calculates conspecific, heterospecific, and facilitator densities
compares raw and `log1p` density transformations
compares Poisson, negative binomial, zero-inflated Poisson, and zero-inflated negative binomial models
compares fixed-effect structures using AIC
tests targeted facilitation interactions
runs DHARMa diagnostics
saves final model objects for plotting
The random-effect structure uses site × cover × block:
```r
(1 | SiteCovBlock)
```
This keeps block numbers separate across sites and cover treatments.
Main outputs:
```text
glmmTMB/glmm_outputs/
glmmTMB/glmm_outputs/model_objects/
```

3. Bayesian coexistence analysis
---
   Scripts in `coexistence_models/scripts/` estimate competition, facilitation, invasion growth, and coexistence outcomes for `TRCY` and `TROR`.
`00_germination_priors_CoexModels.R`
Creates site- and species-specific germination priors.
Main input:
```text
data/field/collated_data_18_02_2020.csv
```
Germination fraction is calculated as:
```text
field_number_germinants / field_number_seeds
```
Direct germination evidence is used for PJ and BEN. Weak priors are created for GH, NAM, and CS.
Main output:
```text
coexistence_models/00_germination_priors/selection/germination_priors.rds
```
The main coexistence workflow can be run with or without these germination priors.
`02_model_selection_error_selection_CoexModels.R`
Runs the main model-selection step for the coexistence analysis.
The script compares Beverton-Holt and Ricker population-growth models crossed with four error models:
`gaussian_raw`
`poisson_count`
`negbin_count`
`gaussian_logfit`
Conspecific and heterospecific densities are kept as raw counts so that competition coefficients are interpretable per competitor individual.
Facilitator density is also kept raw, so:
```text
NO_GRASS = 0
```
means no facilitator.
Main outputs:
```text
coexistence_models/00_family_error_selection/
```
`01b_alpha_facilitation_sensitivity.R`
Tests whether facilitator density affects competition coefficients, rather than only low-density growth.
The script compares:
`m1` = facilitator affects lambda only
`m4` = facilitator affects lambda + interspecific competition
`m5` = facilitator affects lambda + intra- and interspecific competition
Main outputs:
```text
coexistence_models/01b_alpha_facilitation_sensitivity/
```
Equivalent Ricker files are produced if the Ricker model is selected.

### `03_plot_outcomes_BH_Ricker.R`

Calculates coexistence outcomes and generates final coexistence figures for Beverton-Holt and Ricker model families.

The script:

- loads selected Bayesian model fits
- calculates lambda and alpha posterior summaries
- calculates invasion growth
- estimates niche and fitness differences
- classifies coexistence outcomes
- compares BH and Ricker outcome probabilities
- saves final tables and figures

Default setting:

```r
use_germination <- FALSE
```
This means the main outcomes do not apply germination correction unless this option is changed.

Main outputs:

```text
coexistence_models/02_outcomes/tables/
coexistence_models/02_outcomes/draws/
coexistence_models/02_outcomes/figures/
coexistence_models/02_outcomes/summaries/
```
## Recommended run order

```r
# Bayesian coexistence analysis
source("coexistence_models/scripts/00_germination_priors_CoexModels.R")
source("coexistence_models/scripts/02_model_selection_error_selection_CoexModels.R")
source("coexistence_models/scripts/01b_alpha_facilitation_sensitivity.R")
source("coexistence_models/scripts/03_plot_outcomes_BH_Ricker.R")
```

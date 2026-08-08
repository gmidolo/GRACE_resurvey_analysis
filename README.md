# Data and R code for *Warming and protected areas shape long-term plant biodiversity change in Central European grasslands*

------------------------------------------------------------------------

This repository includes the **data** and **R code** used in our study.
It allows for the reproduction of the main analyses and supplementary
materials.

We analyzed the GRACE resurvey dataset ([Haring et al., 2026](#citation))
consisting of 411 survey-resurvey vegetation plot pairs of grassland
vegetation in Austria, Czechia, and Slovakia covering 95 years of
timespan. We test how climate warming, landscape change, and
protected-area status shape long-term change in plant diversity.


## Table of Contents

-   [Authors and Contact](#authors-and-contact)
-   [Acknowledgements](#acknowledgements)
-   [Data Availability](#data-availability)
-   [Data (`data` folder)](#data-data-folder)
-   [R Code (`src` folder)](#r-code-src-folder)
    -   [1. Data preparation](#1-data-preparation)
    -   [2. Analysis](#2-analysis)
    -   [Helper scripts](#helper-scripts)
-   [Metadata:
    `data/model_data/model_input_data.csv`](#metadata-datamodel_datamodel_input_datacsv)
-   [License](#license)
-   [Citation](#citation)

## Authors and Contact

### Contact

**Gabriele Midolo** Department of Spatial Sciences Faculty of
Environmental Sciences Czech University of Life Sciences Prague, Kamýcká
129, 165 00 Praha - Suchdol, Czech Republic ORCID:
<a href="https://orcid.org/0000-0003-1316-2546" target="_blank">0000-0003-1316-2546</a>
Email: `midolo [at] fzp.czu [dot] cz`

### Authors

Gabriele Midolo<sup>1</sup>
<a href="https://orcid.org/0000-0003-1316-2546" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Verena Haring<sup>2</sup>
<a href="https://orcid.org/0009-0007-2564-9315" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Adam T. Clark<sup>2</sup>
<a href="https://orcid.org/0000-0002-8843-3278" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Emma Shih Mendez<sup>3</sup>
<a href="https://orcid.org/0009-0000-6031-5416" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Hana Skokanová<sup>4</sup>
<a href="https://orcid.org/0000-0001-9677-2264" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Martina Sychrová<sup>4</sup>
<a href="https://orcid.org/0000-0002-9503-2544" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Franz Essl<sup>3</sup>
<a href="https://orcid.org/0000-0001-8253-2112" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Milan Chytrý<sup>5</sup>
<a href="https://orcid.org/0000-0002-8122-3075" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Stefan Dullinger<sup>6</sup>
<a href="https://orcid.org/0000-0003-3919-0887" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Jozef Šibík<sup>7</sup>
<a href="https://orcid.org/0000-0002-5949-862X" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Dariia Borovyk<sup>8,9</sup>
<a href="https://orcid.org/0000-0001-7140-7201" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Klára Friesová<sup>5</sup>
<a href="https://orcid.org/0000-0002-1644-2140" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Georg Hörmann<sup>2</sup>
<a href="https://orcid.org/0009-0000-7451-2074" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Eliška Jagošová<sup>1</sup>, Veronika Kalusová<sup>5</sup>
<a href="https://orcid.org/0000-0002-4270-321X" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Matouš Marek<sup>1</sup>
<a href="https://orcid.org/0009-0006-4911-5361" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Julian Pleyer<sup>6</sup>, Zdenka Preislerová<sup>5</sup>
<a href="https://orcid.org/0000-0003-1288-7609" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Helena Schwaiger<sup>10</sup>, Vojtěch Sobotka<sup>5</sup>
<a href="https://orcid.org/0000-0002-0668-3258" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Marie Vymazalová<sup>4</sup>
<a href="https://orcid.org/0000-0002-8842-9156" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Friederike J. R. Wölke<sup>1,11,12</sup>
<a href="https://orcid.org/0000-0001-9034-4883" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Petr Keil<sup>1</sup>
<a href="https://orcid.org/0000-0003-3017-1858" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>

### Affiliations

1.  Department of Spatial Sciences, Faculty of Environmental Sciences,
    Czech University of Life Sciences Prague, Prague - Suchdol, Czech
    Republic
2.  Department of Biology, University of Graz, Graz, Austria
3.  Division of BioInvasions, Global Change & Macroecology, Department
    of Botany and Biodiversity Research, University of Vienna, Vienna,
    Austria
4.  Landscape Research Institute, Průhonice, Czech Republic
5.  Department of Botany and Zoology, Faculty of Science, Masaryk
    University, Brno, Czech Republic
6.  Division of Biodiversity Dynamics and Conservation, Department of
    Botany and Biodiversity Research, University of Vienna, Vienna,
    Austria
7.  Plant Science and Biodiversity Center, Slovak Academy of Sciences,
    Bratislava, Slovak Republic
8.  Theoretical Ecology, Institute of Biology, Free University of
    Berlin, Berlin, Germany
9.  Department of Geobotany and Ecology, M.G. Kholodny Institute of
    Botany, National Academy of Sciences of Ukraine, Kyiv, Ukraine
10. Ingenieurbüro für Biologie Schwaiger, Linz, Austria
11. Department of Ecology, Environment and Plant Sciences, Stockholm
    University, Stockholm, Sweden
12. Bolin Centre for Climate Research, Stockholm University, Stockholm, Sweden

## Acknowledgements

The GRACE project is supported by the Czech Science Foundation grant
24-14299L and the Austrian Science Foundation (FWF) grant I 6578 (DOI:
[10.55776/I6578](https://doi.org/10.55776/I6578)).

## Data Availability

The current repository can be used to reproduce the analyses presented
in this paper. The resurvey dataset of the GRACE project (Haring et al., 2026) 
is openly available under a CC BY license at the following GitHub repository: 
[adamtclark/GRACE_resurveydata](https://github.com/adamtclark/GRACE_resurveydata).

------------------------------------------------------------------------

## Data (`data` folder)

The `data` folder contains the raw, intermediate, and processed data
used in the analyses. Some large raster/netCDF files listed below (e.g.
CHELSA climate rasters, per-plot landscape buffer rasters) are
intermediate products generated by the scripts in
[`src/1.data.preparation`](src/1.data.preparation) and may not be
distributed in full in the public repository due to file size - yet,
they can be regenerated from source using the corresponding scripts.

-   **`data/raw_veg_data/`** - Contains the shared vegetation-plot data
    of the GRACE resurvey dataset.
    -   `EVA&GRACE_header_20260202.csv`: Plot-level metadata
        (coordinates, sampling year, plot size, surveyor, dataset,
        country) for historic and resurvey plots. Refer to Haring et al.
        (2026) for additional details.
    -   `EVA&GRACE_species_20260202.csv`: Species occurrence/cover
        records in the vegetation layers per plot and survey (historical
        and resurvey). Refer to Haring et al. (2026) for additional
        details.
    -   `surveyor_name_IDs.csv`: Lookup table linking surveyor names to
        surveyor IDs.
    -   `Indicator.values-tables-2022-11-07-Zenodo.v2.xlsx`: Tichý et
        al. 2023 (DOI:
        [10.1111/jvs.13168](https://doi.org/10.1111/jvs.13168)).
        Elleberg-like indicator values for light [L], temperature [T],
        moisture [M], reaction [R], nitrogen [N]) used to compute
        community-mean indicator values (CM<sub>EIV</sub>).
-   **`data/landscape/`** - Landscape-structure source and intermediate
    data.
    -   `buffer_union_shape.rds`, `buffer_total.union_shape.rds`:
        Merged buffer polygons per plot/land-cover class and per plot,
        respectively.
    -   `buffer_rasters/`: 10 m resolution rasterized land-cover
        buffers, one multi-layer GeoTIFF per plot (layers correspond to
        19th century vs current land-cover periods).
-   **`data/climate/`** - CHELSA-derived climate data.
    -   `chelsa_tas_gs_annual_1980_2018.tif`: Stacked annual
        growing-season temperature raster (all years, 1979–2021).
    -   `chelsa_tas_trend_1980_2018.tif`: Per-cell linear warming trend
        (°C/year) from 1979 to 2021.
-   **`data/protected_areas/`** - Protected-area (WDPA) data.
    -   `protected_areas.rds`: Cropped protected-area polygons for the
        study region.
    -   `pa.file_GRACE_edit.csv`: Manually edited table used to define
        the "strict" protected-area subset (exclusion of weaker
        protection designations).
-   **`data/predictors/`** - Derived, plot-level predictor tables
    (compiled from the above sources; used for final data preparation).
    -   `bioclim_20260202.csv.gz`: Elevation, growing-season warming
        trend, and CHELSA bioclimatic variables at plot resurvey coordinates.
    -   `class.level_landscapemetrics_20260202.csv.gz`: Class-level
        landscape metrics (core area, mesh size, aggregation index,
        patch number, edge density) per land-cover class and historic
        period, computed from the buffer rasters.
    -   `pa_20260202.csv.gz`: Plot-level protected-area status (including all
        protection designations and strict subset used in the main analysis).
-   **`data/model_data/model_input_data.csv`** - Core dataset used for the main 
    analyses and modeling (411 plot pairs × 39 variables); it combines species
    richness/compositional change metrics with environmental predictors. This is the main output of
    [`src/1.data.preparation/1.F_finalize.data.for.modeling.R`](src/1.data.preparation/1.F_finalize.data.for.modeling.R).
    See the [metadata](#metadata-datamodel_datamodel_input_datacsv) for more detail.
-   **`data/output_brms_models/`** - Contains fitted model objects
    (from `brms` package; `.rds` format) for all main and sensitivity models. 
    In addition, it contains the following:
    -   `main_results_table.csv`: Posterior summaries (median, 95% HDI,
        probability of direction) of fixed effects for the four main
        models.
    -   `table_s1_morans_I.csv`: Moran's I test statistics for residual
        spatial autocorrelation of the main models.
-   **`data/output_pseudoturnover_sim/`** - Output of the
    observation-error (pseudoturnover) sensitivity simulation for the
    Jaccard dissimilarity metric.
    -   `jaccard_error_inflation.csv`: Mean Jaccard inflation,
        regression intercept/slope, and RMSE across a grid of
        false-negative/misidentification error rates.
    -   `jaccard_coef_sensitivity.csv`: Regression coefficients of the
        Jaccard model refitted under each simulated error scenario,
        compared against the baseline (observed data) model.

`fig/` contains the main-text and supplementary figures (`.pdf` or `.png`) generated by the scripts in `src/2.analysis`.

------------------------------------------------------------------------

## R Code (`src` folder)

The `src` folder contains the R scripts organized by their purpose, following the numbering `<stage>.<step>_<description>.R`.

### 1. Data preparation

Scripts in [`src/1.data.preparation`](src/1.data.preparation) build the plot-level predictor dataset from initial spatial and partially pre-processed vegetation data.

-   [`1.A_preprocess.landscape.macrostructure.R`](src/1.data.preparation/1.A_preprocess.landscape.macrostructure.R):
    Processes the 1.5 km land-cover-change buffer polygons around each
    site, unions polygons by plot/class/period, reassigns buffers to
    the nearest GRACE/EVA plot, and rasterizes each plot's buffer land
    cover (10 m resolution) for use in landscape-metric calculations.
-   [`1.B_calc.landscape.metrics.R`](src/1.data.preparation/1.B_calc.landscape.metrics.R):
    Calculates class-level landscape metrics from the
    rasterized buffers for each historic land-cover period and exports
    the resulting predictor table (in `data/predictors`).
-   [`1.C_get.clim.change.R`](src/1.data.preparation/1.C_get.clim.change.R):
    Downloads CHELSA v2.1 monthly temperature data for the study region,
    computes annual growing-season (Apr–Sep) means (1979–2021), derives
    a per-cell linear warming-trend raster, and produces supplementary
    visualizations of the spatial/temporal warming pattern.
-   [`1.D_get.bioclim.variables.R`](src/1.data.preparation/1.D_get.bioclim.variables.R):
    Extracts elevation, the CHELSA growing-season warming trend, and
    CHELSA bioclimatic variables at each plot's coordinates, and exports
    the combined bioclimatic predictor table (in `data/predictors`).
-   [`1.E_protected.areas.R`](src/1.data.preparation/1.E_protected.areas.R):
    Intersects plot locations with WDPA protected-area polygons to
    derive binary protection status under two definitions (all
    designations vs. a stricter subset), store results as a predictor table (in `data/predictors`).
-   [`1.F_finalize.data.for.modeling.R`](src/1.data.preparation/1.F_finalize.data.for.modeling.R):
    Assembles the final analysis-ready dataset - computes species
    richness change (logRR<sub>S</sub>), Jaccard dissimilarity, community-mean
    indicator value change (ΔCM<sub>N,M</sub>), and geographic
    relocation distance between historic and resurvey plots. The output is the dataset used in the main analyses
    (`data/model_data/model_input_data.csv`).

### 2. Analysis

Scripts in [`src/2.analysis`](src/2.analysis) fit and evaluate the
statistical models and produce the main text and supplementary
figures/tables.

-   [`2.A_model_fitting.R`](src/2.analysis/2.A_model_fitting.R): Fits
    the Bayesian generalized linear (mixed) models (`brms`) for the four
    main response variables - species richness change
    (logRR<sub>S</sub>), compositional dissimilarity (Jaccard), and community-mean nutrient and moisture
    indicator change (ΔCM<sub>N</sub>, ΔCM<sub>M</sub>) -
    together with their timespan-interaction versions and various supplementary analyses.
-   [`2.B_model_diagnostics_and_tables.R`](src/2.analysis/2.B_model_diagnostics_and_tables.R):
    Runs model adequacy checks (posterior predictive checks), tests
    multicollinearity among predictors (VIF), tests residual spatial
    autocorrelation (Moran's I), and exports the main results table
    (posterior medians, 95% HDI, probability of direction).
-   [`2.C_model_visualization.R`](src/2.analysis/2.C_model_visualization.R):
    Produces the main-text and supplementary figures.
-   [`2.D_modeling_homogenization_temporal.R`](src/2.analysis/2.D_modeling_homogenization_temporal.R):
    Analyzes temporal patterns of biotic homogenization/differentiation
    using Whittaker decomposition as proposed by [Blowes et al. 2024](https://doi.org/10.1126/sciadv.adj9395),
    applied across historic sampling-period time slices and country subgroups.
-   [`2.E_Jaccard_pseudoturnover_simulation.R`](src/2.analysis/2.E_Jaccard_pseudoturnover_simulation.R):
    Quantifies how much of the observed Jaccard dissimilarity could be
    attributable to observation error (false negatives and species
    misidentification) using simulation.
-   [`2.F_Ndeposition_EMEP_trends.R`](src/2.analysis/2.F_Ndeposition_EMEP_trends.R):
    Supplementary-only analysis comparing atmospheric nitrogen
    deposition trends (1990-2024) between Austria and
    Czechia-Slovakia using [EMEP model data](https://www.emep.int/mscw/mscw_moddata.html)).

### Helper scripts

-   [`src/utils.R`](src/utils.R): Contains utility functions / helpers

------------------------------------------------------------------------

## Metadata: `data/model_data/model_input_data.csv`

Final analysis-ready dataset: 411 rows (one per historic–resurvey plot
pair), 39 columns.

| Column Name                    | Description                                                                                                                                                                  |
|:-----------------|:-----------------------------------------------------|
| `plot_id`                      | Unique plot identifier (European Vegetation Archive database; data request no. 147)                                                                                          |
| `rsrv_S`                       | Vascular plant species richness at resurvey                                                                                                                                  |
| `historic_S`                   | Vascular plant species richness at historic survey                                                                                                                           |
| `logRR_S`                      | Log<sub>10</sub> response ratio of species richness: log10(rsrv_S / historic_S)                                                                                              |
| `jaccard`                      | Jaccard dissimilarity between historic and resurvey plant community composition on binary (presence/absence) data                                                            |
| `diff.CM_L`                    | Change in community-mean indicator value for Light (resurvey − historic)                                                                                                     |
| `diff.CM_T`                    | Change in community-mean indicator value for Temperature (resurvey − historic)                                                                                               |
| `diff.CM_M`                    | Change in community-mean indicator value for Moisture (resurvey − historic)                                                                                                  |
| `diff.CM_R`                    | Change in community-mean indicator value for Reaction/soil pH (resurvey − historic)                                                                                          |
| `diff.CM_N`                    | Change in community-mean indicator value for Nutrients/Nitrogen (resurvey − historic)                                                                                        |
| `country`                      | Country grouping: `AT` (Austria) or `CZSK` (Czech Rep & Slovakia)                                                                                                            |
| `rsrv_lon`, `rsrv_lat`         | Resurvey plot coordinates (WGS84)                                                                                                                                            |
| `historic_lon`, `historic_lat` | Historic plot coordinates (WGS84)                                                                                                                                            |
| `dist_meters`                  | Geographic distance (m) between historic and resurvey plot locations                                                                                                         |
| `rsrv_sampling_year`           | Year of resurvey sampling                                                                                                                                                    |
| `historic_sampling_year`       | Year of historic sampling                                                                                                                                                    |
| `diff_plot_size`               | Difference in plot size (m²): rsrv_plot_size − historic_plot_size                                                                                                            |
| `timespan`                     | Years elapsed between historic and resurvey sampling                                                                                                                         |
| `rsrv_plot_size`               | Plot size (m²) at resurvey                                                                                                                                                   |
| `historic_plot_size`           | Plot size (m²) at historic survey                                                                                                                                            |
| `rsrv_surveyor_name`           | Name of surveyor(s) at resurvey                                                                                                                                              |
| `historic_surveyor_name`       | Name of surveyor(s) at historic survey                                                                                                                                       |
| `EVA_dataset`                  | Source EVA dataset name                                                                                                                                                      |
| `EVA_country`                  | Country as recorded in the EVA database                                                                                                                                      |
| `elevation`                    | Elevation above sea level (m), extracted from the Copernicus GLO-90 DEM at resurvey coordinates                                                                              |
| `gs_temp_change_decade`        | Growing-season temperature warming trend (°C/decade), 1979–2021, extracted from CHELSA v2.1 at resurvey coordinates                                                           |
| `grassland_ai_1840`            | Aggregation index (%) of the grassland land-cover class within the 1.5 km buffer, historic period \~1840                                                                     |
| `grassland_ai_1870`            | Aggregation index (%) of the grassland land-cover class within the 1.5 km buffer, historic period \~1870                                                                     |
| `grassland_ai_2023`            | Aggregation index (%) of the grassland land-cover class within the 1.5 km buffer, current period (2023)                                                                      |
| `grassland_ca_1840`            | Total class (core) area (ha) of grassland within the 1.5 km buffer, \~1840                                                                                                   |
| `grassland_ca_1870`            | Total class (core) area (ha) of grassland within the 1.5 km buffer, \~1870                                                                                                   |
| `grassland_ca_2023`            | Total class (core) area (ha) of grassland within the 1.5 km buffer, 2023                                                                                                     |
| `grassland_np_1840`            | Number of grassland patches within the 1.5 km buffer, \~1840                                                                                                                 |
| `grassland_np_1870`            | Number of grassland patches within the 1.5 km buffer, \~1870                                                                                                                 |
| `grassland_np_2023`            | Number of grassland patches within the 1.5 km buffer, 2023                                                                                                                   |
| `pa_site`                      | Binary protected-area indicator (0/1): plot falls within a WDPA protected area, any designation                                                                              |
| `pa_site_restrict`             | Binary protected-area indicator (0/1): plot falls within a WDPA protected area under a stricter subset of designations (used as the main protection covariate in the models) |

N.B. In the modeling scripts
([`2.A_model_fitting.R`](src/2.analysis/2.A_model_fitting.R) onward),
`delta_CA` and `delta_NP` (change in grassland core area and patch
number, 2023 − 1870) are derived directly from `grassland_*`
columns above in the script.

## License

**Vegetation data** are available under the terms of the Creative Commons Attribution 4.0 International (CC BY) - see [Data
Availability](#data-availability).

**Code** in this repository is available under the terms of the GNU
General Public License v3.0 (GPL-3.0)
(<https://www.gnu.org/licenses/gpl-3.0.html>).

## Citation

#### Paper citation

Midolo, G., Haring, V., Clark, A.T., Shih Mendez, E., Skokanová, H.,
Sychrová, M., Essl, F., Chytrý, M., Dullinger, S., Šibík, J., Borovyk,
D., Friesová, K., Hörmann, G., Jagošová, E., Kalusová, V., Marek, M.,
Pleyer, J., Preislerová, Z., Schwaiger, H., Sobotka, V., Vymazalová, M.,
Wölke, F.J.R., Keil, P. Warming and protected areas shape long-term
plant biodiversity change in Central European grasslands. *In preparation/Under review*.

#### Data citation

Haring, V., Midolo, G., Keil, P., Chytrý, M., Essl, F., Mendez, E. S., ... & Clark, A. T. (2026).
GRACE Resurvey of Central European Grassland Vegetation to study the Impacts of Country-Level Landscape Dynamics on Biodiversity.
*EcoEvoRxiv*. <https://doi.org/10.32942/X2164H>
# Bioaccumulation in Aquatic Organisms: From One-Compartment Toxicokinetic Models to Real-World Complexity

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Code: GPL v3](https://img.shields.io/badge/code-GPL%20v3-blue.svg)](LICENSE.md)
[![Content: CC BY 4.0](https://img.shields.io/badge/content-CC%20BY%204.0-lightgrey.svg)](LICENSE.md)

R scripts, curated datasets and teaching material for a half-day, hands-on workshop on **toxicokinetic (TK) modelling** of organic contaminants in aquatic invertebrates.

**University of Waterloo, 12 November 2026** · Johannes Raths, Eawag

Toxicokinetic processes (uptake, elimination and biotransformation) set the internal concentration of an exposed organism and so link exposure to effects. Through short lectures and guided exercises on real experimental data, the workshop starts from the classic one-compartment model and extends it step by step towards real-world complexity.

## What you will learn

- Fit TK models and derive parameters such as kinetic rate constants and bioconcentration factors (BCFs)
- Model biotransformation with parent–metabolite models
- Explore multi-compartment kinetics, such as elimination resistance through receptor binding, and lipid partitioning
- Add dietary uptake and your own processes to a TK model
- Account for the effect of temperature on toxicokinetic processes

No prior experience with bioaccumulation or programming is needed.

## Exercises

| # | Script | Model | Dataset | Run time* |
|---|--------|-------|---------|-----------|
| 1 | `TK_1_one_compartment.R` | Standard one-compartment model (ku, ke) | 12 contaminants in *Gammarus pulex* and *Hyalella azteca* at 6–21 °C (Raths et al. 2023, *Glob. Change Biol.*) | ~6 s |
| 2 | `TK_2_biotransformation.R` | Parent + biotransformation product (km, kem) | Diflufenican and its metabolite AE B107137 in *Asellus aquaticus* (Raths et al., submitted) | ~30 s |
| 3a | `TK_3a_receptor_binding.R` | Two compartments: structure + saturable receptor binding (kon, koff, Bmax) | Thiacloprid in *G. pulex* (Raths et al. 2023, *Environ. Sci. Technol.*) | ~12 min (simplex only: ~2 min) |
| 3b | `TK_3b_lipid_partitioning.R` | Two compartments: structure + lipid | Dimethylnaphthalene in a marine copepod (Jager et al. 2017, *Environ. Sci. Technol.*) | ~30 s |
| 4 | `TK_4_dietary_uptake.R` | One compartment with uptake from water and food | Fluopyram in caged vs. uncaged *G. pulex* (Raths et al. 2023, *Environ. Toxicol. Chem.*) | ~21 s |

\*Default settings, on a standard laptop.

Each script has its matching data file (`*_data.xlsx`) in the same folder. `TK_Template.xlsx` shows the data format if you want to fit **your own data**.

`Scripts/Solutions/` contains reference results with the default settings (plots, logs and Excel output), so you can check your own runs against them.

## Getting started

### 1. Install R and RStudio

- R (version 4.3 or newer): <https://cran.r-project.org>
- RStudio Desktop: <https://posit.co/download/rstudio-desktop/>

### 2. Download the material

- **Easiest:** on the [Zenodo record](https://doi.org/10.5281/zenodo.XXXXXXX), download the ZIP file and unzip it.
- **Or** on GitHub: click **Code → Download ZIP**, or clone the repository:
  ```bash
  git clone https://github.com/HollenderGroup/Teaching_TK_Workshop_Waterloo_12.11.2026.git
  ```

### 3. Install the packages

Open RStudio and run `Scripts/install_packages.R` once (open it, then click **Source**). The scripts use base R graphics and only need four small packages:

| Package | Used for |
|---------|----------|
| `readxl` | reading the Excel data files |
| `writexl` | writing the Excel result files |
| `deSolve` | numerical ODE solver (receptor binding, and your own model extensions) |
| `rstudioapi` | the file browser for selecting the data file |

### 4. Run an exercise

1. Open a script, e.g. `Scripts/TK_1_one_compartment.R`, in RStudio.
2. Change settings in **SECTION 1** only (compound, scenarios, fitting options). Nothing below needs editing.
3. Click **Source**. A file browser opens. Select the matching `*_data.xlsx` file.
4. The results are written to a time-stamped folder **next to the data file**:
   - `..._medium_tissue_concentrations.png`: data and model fit
   - `..._parameter_space.png`: parameter-space explorer and profile likelihoods
   - `..._log.txt`: everything printed in the console
   - `<RUN_NAME>_<time stamp>_TK_results.xlsx`: all parameters with 95% CIs, goodness-of-fit metrics, and the data and model curves

How the models are fitted, how to read the parameter-space plot, and how to add your own uptake pathway are explained in plain language in [`docs/Model_mathematics.md`](docs/Model_mathematics.md).

## Workshop schedule (12 November 2026)

| Time | Session |
|------|---------|
| 13:00–13:15 | Welcome: why bioaccumulation matters |
| 13:15–13:30 | Bioaccumulation as a kinetic problem (OECD 305/321, kinetic vs. steady-state BCF) |
| 13:30–14:15 | **Exercise 1:** The standard one-compartment TK model, plus discussion |
| 14:15–15:00 | **Exercise 2:** Biotransformation, plus discussion |
| 15:00–15:15 | Coffee break |
| 15:15–16:00 | **Exercise 3:** Multi-compartment kinetics (receptor binding, lipid partitioning), plus discussion |
| 16:00–16:30 | **Outlook:** Dietary uptake, temperature effects and your own ideas and data |
| 16:30–17:00 | Final discussion: choosing the right model, linking TK to toxicity, feedback |

## Repository structure

```
├── README.md
├── TK UW - All Slides.pdf        # lecture slides (PDF)
├── TK UW - All Slides.pptx       # lecture slides (editable)
├── TK UW - Workshop Schedule.docx
├── TK UW - Model mathematics for dummies.docx
├── Scripts/
│   ├── install_packages.R        # run once to install the R packages
│   ├── TK_1_one_compartment.R    # exercise scripts …
│   ├── TK_1_one_compartment_data.xlsx   # … each with its dataset
│   ├── …
│   ├── TK_Template.xlsx          # data template for your own data
│   └── Solutions/                # reference results
├── docs/
│   └── Model_mathematics.md      # how the scripts work, in plain language
├── CITATION.cff                  # citation metadata (GitHub "Cite this repository")
├── .zenodo.json                  # metadata for the Zenodo archive
├── LICENSE.md                    # GPL v3 (code)
└── BYOM_LICENSE.txt              # original BYOM license (MIT)
```

## Background and credits

The scripts are a simplified R re-implementation of the BYOM model **TK-metabolism** (<https://github.com/HollenderGroup/Bioaccumulation_TK-metabolism>), built on the BYOM v6.2 framework by Tjalling Jager (<https://www.debtox.info/byom.html>). They include the parameter-space explorer from openGUTS (<https://www.openguts.info>) as implemented in BYOM, and the simplex optimisation. Fitting and confidence intervals follow BYOM. The goodness-of-fit criteria are TK analogues of the EFSA (2018) validation criteria.

### Data sources

Please cite the original studies when you use the datasets:

- Raths, J., Švara, V., Lauper, B., Fu, Q., & Hollender, J. (2023). Speed it up: How temperature drives toxicokinetics of organic contaminants in freshwater amphipods. *Global Change Biology*, 29(5), 1390–1406. <https://doi.org/10.1111/gcb.16542>
- Raths, J., Stäuble, C., Hermann, M., & Hollender, J. From laboratory to mesocosm: Temperature and species dependent toxicokinetics of the fluorinated herbicide diflufenican. *ES&T Toxicology*, submitted.
- Raths, J., Schinz, L., Mangold-Döring, A., & Hollender, J. (2023). Elimination resistance: Characterizing multi-compartment toxicokinetics of the neonicotinoid thiacloprid in the amphipod *Gammarus pulex* using bioconcentration and receptor-binding assays. *Environmental Science & Technology*. <https://doi.org/10.1021/acs.est.3c01891>
- Jager, T., Øverjordet, I. B., Nepstad, R., & Hansen, B. H. (2017). Dynamic links between lipid storage, toxicokinetics and mortality in a marine copepod exposed to dimethylnaphthalene. *Environmental Science & Technology*, 51(13), 7707–7713. <https://doi.org/10.1021/acs.est.7b02212>
- Raths, J., Schnurr, J., Bundschuh, M., Pinto, F. E., Janfelt, C., & Hollender, J. (2023). Importance of dietary uptake for in situ bioaccumulation of systemic fungicides using *Gammarus pulex* as a model organism. *Environmental Toxicology and Chemistry*, 42(9), 1993–2006. <https://doi.org/10.1002/etc.5615>

## How to cite this material

> Raths, J. (2026). *Bioaccumulation in Aquatic Organisms: From One-Compartment Toxicokinetic Models to Real-World Complexity. Workshop material, University of Waterloo* (Version 1.0.0) [Lesson]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

GitHub shows the citation under **"Cite this repository"** (from `CITATION.cff`).

## License

- **Code** (R scripts): [GNU General Public License v3.0 or later](LICENSE.md). The scripts re-implement the parameter-space explorer of openGUTS/BYOM, which is GPL v3, so they carry the same license. Parts derived from the general BYOM framework are also covered by the BYOM MIT license, reproduced in [`BYOM_LICENSE.txt`](BYOM_LICENSE.txt) (Copyright (c) 2012–2020 Tjalling Jager).
- **Teaching content** (slides, schedule, documentation): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- **Data:** compiled from the publications listed above. Please cite them.

## Contact

Questions, errors or ideas: please [open an issue](https://github.com/HollenderGroup/Teaching_TK_Workshop_Waterloo_12.11.2026/issues) or contact Johannes Raths (Eawag, Department of Environmental Chemistry).

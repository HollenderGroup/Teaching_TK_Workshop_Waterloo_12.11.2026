###############################################################################
###############################################################################
##
##   TK WORKSHOP - EXERCISE 2
##   ONE-COMPARTMENT TK MODEL WITH BIOTRANSFORMATION (TK-metabolism)
##
##   medium  --ku-->  parent  --km-->  biotransformation product (BTP)
##                      | ke                 | kem
##                      v                    v
##                    medium               medium
##
##   TK Workshop "Bioaccumulation in Aquatic Organisms", University of
##   Waterloo, 12 November 2026 (Johannes Raths, Eawag).
##
##   R re-implementation of the BYOM model "TK-metabolism":
##       https://github.com/HollenderGroup/Bioaccumulation_TK-metabolism
##   (one-compartment TK model with time-varying exposure and phase-one
##    biotransformation; BYOM v6.2 framework by Tjalling Jager; model by
##    Johannes Raths, 2026; for the Matlab version see Wehrli et al. 2026,
##    "Pesticides Never Walk Alone: Mixture Toxicity and Toxicokinetics of
##    Five Formulated Pesticide Products in the Soil Arthropod Folsomia candida").
##   It includes the parameter-space explorer (openGUTS/BYOM "parspace"
##   algorithm) and the simplex optimisation.
##
##   Set up for: diflufenican (DFF) and its metabolite AE B107137 in
##               Asellus aquaticus, experiment AE1
##   Data file : TK_2_biotransformation_data.xlsx
##   Data from : Raths, J.; Stauble, C.; Hermann, M.; Hollender, J. From
##               laboratory to mesocosm: Temperature and species dependent
##               toxicokinetics of the fluorinated herbicide Diflufenican.
##               ES&T Toxicology, submitted October 2026.
##
###############################################################################
###############################################################################
##
##  THE MODEL (identical to BYOM derivatives.m)
##  -------------------------------------------
##      dCp/dt  = ku * Cw_p(t)  -  ke * Cp  -  km * Cp
##      dCm/dt  = km * Cp       -  kem * Cm   [ + kum * Cw_m(t) , optional ]
##
##   Cp     internal concentration of the parent        [umol/kg ww]
##   Cm     internal concentration of the BTP           [umol/kg ww]
##   Cw_p   parent concentration in the medium          [umol/L]
##   Cw_m   BTP concentration in the medium (only used if
##          BTP_UPTAKE_FROM_MEDIUM = TRUE; not part of the BYOM model)
##   ku     uptake rate constant of the parent          [L/kg ww/d]
##   ke     elimination rate constant of the parent     [1/d]
##   km     biotransformation rate constant             [1/d]
##   kem    elimination rate constant of the BTP        [1/d]
##
##   The parent is now lost in two ways, so at steady state
##      BCF = Cp / Cw = ku / (ke + km)                  [L/kg ww]
##
##   WHY MOL? One molecule of parent gives one molecule of BTP, but the BTP
##   is lighter (283.20 instead of 394.29 g/mol). The mass balance between
##   parent and BTP is therefore only correct in molar units: the script
##   converts all concentrations from ug to umol before fitting
##   (USE_MOLAR_UNITS, SECTION 1.4). The plots show the units of the data file.
##
##   Without biotransformation (USE_BIOTRANSFORMATION = FALSE): km = 0, the
##   BTP is not modelled and ke describes the total elimination of the parent
##   (= the one-compartment model of Exercise 1).
##
##   The medium concentration Cw(t) is interpolated linearly between the
##   measured time points, and kept constant after the last measurement.
##   Repeated time points in the medium data (e.g. 1.15 and 0 at
##   t = 2.00125 d) are read as a STEP change in file order (end of
##   exposure, transfer to clean medium).
##
##  HOW TO USE THIS SCRIPT
##  ----------------------
##  1. All data are in ONE file (long format, one row per sample), with
##     experiment code, matrix, species, temperature/treatment, time,
##     replicate, compound, parent and BTP concentration
##     (TK_2_biotransformation_data.xlsx).
##  2. In SECTION 1.3 choose the data. A scenario is one combination of
##     "Experiment code" and "Temperature [C] / Treatment".
##     SEPARATE_MODELS: every selected scenario is calibrated as its OWN model.
##     POOLED_MODELS: scenarios calibrated TOGETHER with one parameter set.
##     The script prints a list of all available scenarios at the start.
##  3. Run the whole script (in RStudio: "Source").
##     With DATA_FILE <- "" a file browser opens to select the data file.
##     Progress bars show how far the fitting is.
##  4. Results are written to OUTPUT_DIR (next to the data file):
##      - ONE folder PER MODEL (time-stamped) with the plots and the log:
##          ..._medium_tissue_concentrations.png   data and model fit
##          ..._parameter_space.png                parameter-space explorer
##          ..._compartment_concentrations.png    each model compartment in its
##                                                own panel (only models with
##                                                more than one curve)
##          ..._log.txt                            everything printed in the console
##      - ONE Excel file for the WHOLE RUN, <RUN_NAME>_<time stamp>_TK_results.xlsx:
##          'TK_parameters'    all fitted and derived parameters, with 95% CIs
##          'quality_of_fit'   MLL, AIC, R2, PPC, NRMSE, PE per model
##          'data_and_curves'  observations and model curves, for your own plots
##
##  You do NOT need to change anything below SECTION 1.
##
##  HOW THE PARAMETERS ARE FITTED (identical to BYOM transfer.m)
##  ------------------------------------------------------------
##   For each fitted state (parent, BTP): data and model are transformed
##   (default square root), residuals are summed over all replicates, time
##   points and pooled scenarios, and the residual variance is profiled out:
##       -logL_state = (n/2) * log(SSQ) + n/2
##   The total minus-log-likelihood (MLL) is the sum over the states.
##   The 95% confidence intervals come from the profile likelihoods (see the
##   parameter-space plot; how to read it is explained in SECTION 8.5).
##
##  GOODNESS-OF-FIT METRICS (per state; SECTION 8.2)
##  ------------------------------------------------
##   R2      1 - SSQ/SST on the transformed scale (as BYOM)
##   AIC     2*MLL + 2*(number of fitted parameters); to compare models
##           fitted to the SAME data (lower = better)
##   PPC     % of observed tissue concentrations within the 95% prediction
##           interval (acceptable >= 50%)
##   NRMSE   root-mean-square error / mean of observations x 100%
##           (acceptable <= 50%)
##   PE      prediction error at the end of exposure:
##           (predicted - observed) / observed x 100%  (acceptable if <= 20%)
##   (toxicokinetic analogues of the three validation criteria of EFSA, 2018)
##
##  DIFFERENCES TO THE BYOM (Matlab) VERSION
##  ----------------------------------------
##   * The ODEs are solved exactly (analytical solution, SECTION 5) instead
##     of with ode113/ode15s; same result, ~100x faster, no solver noise.
##   * Fixed a BYOM 6.2 bug in the extra sampling after profiling (the jump
##     size was min - min = 0, so "extra" sets were copies of existing ones).
##   * Data come from one long-format file. Exposure profiles are built from
##     the medium rows (mean of the replicates per time point). Text entries
##     such as "<LOQ" are set to 0.
##   * Initial internal concentrations are taken from the tissue data at
##     t = 0 (BYOM: set by hand in X0mat, usually 0).
##   * Optional: BTP uptake from the medium (kum), not in BYOM.
##   * Added: BCF with CI, PPC, NRMSE, PE.
##
###############################################################################


###############################################################################
###############################################################################
##                                                                           ##
##   SECTION 1 :  USER SETTINGS   (this is the only part you need to edit)   ##
##                                                                           ##
###############################################################################
###############################################################################

## ---------------------------------------------------------------------------
## 1.1  Files and folders
## ---------------------------------------------------------------------------
# DATA_FILE: ""  = open a file browser to select the file (recommended)
#            or a path, e.g. "C:/Users/me/DFF/raths_dff_SI_B01.xlsx"

DATA_FILE     <- ""         # "" = choose with file browser; .xlsx / .xls / .csv
OUTPUT_DIR    <- "TK_2_biotransformation_output"   # results folder. A relative name is created
                                    # NEXT TO THE DATA FILE; or give a full path.
EXCEL_SHEET   <- 1          # Excel only: sheet number (1, 2, ...) or name;
# "all" = combine all sheets with the columns below
CSV_SEPARATOR <- "auto"     # CSV only: "auto", "," or ";"
CSV_DECIMAL   <- "."        # CSV only: "." or ","

RUN_NAME      <- "DFF"              # prefix for all output files


## ---------------------------------------------------------------------------
## 1.2  Column names in the data file
## ---------------------------------------------------------------------------
COL_EXPERIMENT <- "Experiment code"               # e.g. TK, AE1, I
COL_MATRIX     <- "Matrix"                        # e.g. Medium, Tissue
COL_SPECIES    <- "Species"                       # e.g. A. aquaticus
COL_TREATMENT  <- "Temperature [C] / Treatment"   # e.g. 16, 20.4, 50, Caged
COL_TIME       <- "Time [d]"                      # time in DAYS
COL_REPLICATE  <- "Replicate"
COL_COMPOUND   <- "Compound"                      # name of the compound (COMPOUND, SECTION 1.4)
COL_PARENT     <- "Parent concentration"          # concentration of the parent
COL_BTP        <- "BTP concentration"             # BTP (metabolite) concentration
COL_UNIT       <- "Unit"                          # optional (labels only)


## ---------------------------------------------------------------------------
## 1.3  Data selection: MODEL RUNS (which scenarios calibrate which model)
## ---------------------------------------------------------------------------
# EXPOSURE_MATRIX : value in the Matrix column for the exposure medium (ug/L)
# ORGANISM_MATRIX : value in the Matrix column for the internal (tissue)
#                   concentrations that are fitted (ug/kg ww)
# Other matrices (e.g. Macrophyte, Leaf) are ignored.

EXPOSURE_MATRIX <- "Medium"
ORGANISM_MATRIX <- "Tissue"

# A SCENARIO is one combination of experiment code and temperature/treatment,
# written as "CODE | TREATMENT". If several species share a code and
# treatment, the species is added: "CODE | TREATMENT | SPECIES".
# The list of all scenarios is printed at the start of the run.
#
# Scenarios in TK_2_biotransformation_data.xlsx:
#   "AE0 | 20.1" (pretest, 3 d uptake) and "AE1 | 20.4" (2 d uptake, 2 d elimination)
#
# Which models are calibrated? All models below are fitted one after the
# other in ONE run of the script (each with its own output folder).
#
# SEPARATE_MODELS: every scenario matched here becomes its OWN model.
#   - "AE1"                 : all scenarios of experiment AE1
#   - "AE1 | 20.4"          : one scenario
#   - "all"                 : every scenario in the file, each separate
#   - "choose"              : pick the scenarios in a selection window
#   - c()                   : no separate models
# Scenarios without parent in the medium (e.g. controls) are skipped automatically.
SEPARATE_MODELS <- c("AE1")

# POOLED_MODELS: each entry is ONE model; all its scenarios are calibrated
# together with one parameter set. The entry name is used for the files.
#   An experiment code alone pools all its scenarios.
#   POOLED_MODELS <- "choose" : pick the scenarios of ONE pooled model in a
#                               selection window;  list() = none
POOLED_MODELS <- list(
  # "AE0_AE1_pooled" = c("AE0 | 20.1", "AE1 | 20.4")
)

# Values below the limit of quantification: every text entry in a
# concentration column (e.g. "<LOQ", "<LOD", "n.d.") is set to
# BELOW_LOQ_VALUE. Entries listed in MISSING_CODES and empty cells are
# treated as missing (not used).
BELOW_LOQ_VALUE <- 0
MISSING_CODES   <- c("NA", "n.a.", "na", "-", "missing")


## ---------------------------------------------------------------------------
## 1.4  Compounds and molar masses
## ---------------------------------------------------------------------------
# COMPOUND: the compound to model = its name in the Compound column of the
# data file. Its BTP is in the BTP concentration column; BTP_NAME is its name
# (used for the molar mass and in plots and tables).
#
# WHY THE MOLAR MASS?
# The data file contains the concentrations as measured, in ug/L and ug/kg.
# The model, however, counts MOLECULES: right after reading the data, every
# concentration is divided by the molar mass (ug / (g/mol) = umol). The model
# is fitted in umol/L and umol/kg, and the results are multiplied by the molar
# mass again, so all plots are in ug/L and ug/kg.
# Molar units matter whenever molecules are converted or counted:
#   - biotransformation: 1 molecule of parent gives 1 molecule of BTP, but
#     the BTP has a different mass (Exercise 2);
#   - receptor binding: a fixed NUMBER of binding sites (Exercise 3a).
# In a one-compartment model the units cancel: ku, ke and the BCF are the
# same in ug or umol. Try it: USE_MOLAR_UNITS <- FALSE.
# Please verify the molar masses!

COMPOUND <- "Diflufenican"
BTP_NAME <- "AE B107137"

MOLAR_MASSES <- c(                  # g/mol
  "Diflufenican" =  394.29,   # C19H11F5N2O2, parent
  "AE B107137"   =  283.20    # C13H8F3NO3, BTP
)

USE_MOLAR_UNITS   <- TRUE           # TRUE = fit in umol/L and umol/kg
DATA_UNIT_MEDIUM  <- "ug/L"         # units in the data file (axis labels of the plots)
DATA_UNIT_TISSUE  <- "ug/kg ww"


## ---------------------------------------------------------------------------
## 1.5  Model switches
## ---------------------------------------------------------------------------
# USE_BIOTRANSFORMATION = TRUE  : full model: parent -> BTP (km, kem fitted;
#                                 parent AND BTP tissue data are fitted)
# USE_BIOTRANSFORMATION = FALSE : one-compartment TK of the parent only.
#                                 km = kem = 0, BTP data are ignored, ke is
#                                 the total elimination (incl. biotransf.)
USE_BIOTRANSFORMATION <- TRUE      # c(TRUE, FALSE) = run EVERY model in both
                                    # versions (file names get _bio / _nobio)

# Extra option, NOT in the BYOM model: direct uptake of the BTP from the
# medium (+ kum * Cw_m). Requires BTP data in the medium.
BTP_UPTAKE_FROM_MEDIUM <- FALSE

# Internal (tissue) concentrations at t = 0 = start values of the ODEs:
#   "data"   : mean of the tissue data at t = 0 of each scenario (DEFAULT;
#              0 if there are no tissue data at t = 0)
#   "manual" : use X0_PARENT and X0_BTP below
X0_SOURCE <- "data"
X0_PARENT <- 0
X0_BTP    <- 0


## ---------------------------------------------------------------------------
## 1.6  MODEL PARAMETERS: start values, fit flags and search ranges
## ---------------------------------------------------------------------------
#   name  : parameter name (do not change)
#   start : start value. Simplex starts here. The parameter-space explorer
#           does not need it (it searches the whole min-max range), but it is
#           used for parameters that are NOT fitted (fit = 0).
#   fit   : 1 = estimate this parameter, 0 = keep it fixed at given value
#   min   : lower bound of the search range
#   max   : upper bound of the search range
#   scale :  "lin" = linear scale; "log" = search on log10 scale (recommended when the range spans
#           more than ~1 order of magnitude)
#           (BYOM: 0 = log, 1 = lin)
#
#   Units (time in days; medium in umol/L):
#   ku   uptake rate constant                [L / kg ww / d]
#   ke   elimination rate constant (parent)  [1/d]
#   km   biotransformation rate constant     [1/d]  (parent -> BTP)
#   kem  elimination rate constant of BTP    [1/d]
#   kum  uptake rate constant of BTP from the medium [L / kg ww / d]
#        (only used if BTP_UPTAKE_FROM_MEDIUM = TRUE)
#
#   The parameter-space explorer searches the full min-max range, so keep
#   the ranges wide enough but realistic. If a CI runs into a bound, widen it.

PARAMETERS <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
  name    start     fit    min       max      scale
  ku      100       1      1         5000    log
  ke      1         1      0.01      100     log
  
  km      0.5       1      0.01      100     log
  kem     1         1      0.01      100     log
  
  kum     1         0      0.001     1000    log
")


## ---------------------------------------------------------------------------
## 1.7  Fitting settings
## ---------------------------------------------------------------------------
FIT_PARAMETERS <- TRUE      # FALSE = no fitting, just plot the model with the
                            #         start values (useful to check start values)

FIT_METHOD <- "parspace"    # "parspace" : parameter-space explorer (DEFAULT, as in
                            #              the BYOM scripts; slower, robust, gives
                            #              95% CIs on parameters and model curves)
                            # "simplex"  : Nelder-Mead simplex from the start values
                            #              (fast; best fit + profile-likelihood CIs).
                            #              CAUTION: can get stuck in a local optimum
                            #              (depends on start values); use for quick
                            #              checks, report results from "parspace".

DATA_TRANSFORM <- 0.5       # transformation of data AND model before residuals:
                            #   0.5 = square root (BYOM scripts default)
                            #   0   = log-transformation (zeros are set to 1e-10!)
                            #   1   = no transformation

# Parameter-space explorer (FIT_METHOD = "parspace")
PS_ROUGHNESS <- 1           # 0 = settings as in openGUTS (slow, most thorough)
                            # 1 = "rough" settings (BYOM scripts default)
                            # 2 = very rough (best fit ok, CIs not reliable)
PS_PROFILES  <- TRUE        # TRUE = profile likelihoods + extra sampling (BYOM
                            #        default, needed for reliable 95% CIs)

# Simplex (FIT_METHOD = "simplex")
SIMPLEX_RUNS       <- 2     # number of simplex runs, each restarting from the
                            # previous best (BYOM opt_optim.simno)
SIMPLEX_PROFILE_CI <- TRUE  # TRUE = profile likelihoods for 95% CIs after the
                            #        simplex (takes extra time)

RANDOM_SEED <- 42           # makes the parameter-space explorer reproducible

SOLVER <- "analytical"      # "analytical" : exact solution (fast, DEFAULT)
                            # "ode"        : numerical ODE solver (package deSolve;
                            #                ~100x slower; use to check or when you
                            #                change the model equations)


## ---------------------------------------------------------------------------
## 1.8  Goodness-of-fit metrics
## ---------------------------------------------------------------------------
PPC_MIN   <- 50     # acceptable if >= PPC_MIN % of observations are in the 95% PI
NRMSE_MAX <- 50     # acceptable if NRMSE <= NRMSE_MAX %
PE_MAX    <- 20     # acceptable if |PE| <= PE_MAX % (prediction error at the end of exposure)
# End of exposure for the PE: NULL = automatic (see header), or set it per
# scenario in days, e.g. END_OF_EXPOSURE <- c("AE0 | 20.1" = 3, "AE1 | 20.4" = 2.00125)
END_OF_EXPOSURE <- NULL


## ---------------------------------------------------------------------------
## 1.9  Output settings
## ---------------------------------------------------------------------------
SHOW_PROGRESS   <- TRUE     # show progress bars in the console
SHOW_PLOTS      <- TRUE     # also draw the plots in the R/RStudio plot window
SAVE_PLOTS      <- TRUE     # save plots to OUTPUT_DIR
N_TIME_POINTS_CURVE <- 500  # resolution of the model curves
CI_MAX_SETS     <- Inf      # max. parameter sets used for the 95% CI bands
                            # (Inf = all, as BYOM; e.g. 1000 for speed)
EXCEL_DIGITS    <- 2        # significant digits SHOWN in the Excel output (the cells
                            # keep the full precision; MLL and AIC: 1 decimal)

###############################################################################
###############################################################################
##                                                                           ##
##             END OF USER SETTINGS - no changes needed below                ##
##                                                                           ##
###############################################################################
###############################################################################




###############################################################################
## SECTION 2 : Packages and general helpers
###############################################################################

T_START <- Sys.time()
RUN_STAMP <- format(T_START, "%Y-%m-%d_%H%M%S")   # time stamp of this run (output folders)

# Only base R is required. Optional packages:
#   readxl     : to read Excel files            install.packages("readxl")
#   rstudioapi : nicer file browser in RStudio  install.packages("rstudioapi")
#   openxlsx   : Excel output with readable column widths (recommended);
#                without it, writexl writes plain sheets
#   deSolve    : only for SOLVER = "ode" and for the solver self-check
if (SOLVER == "ode" && !requireNamespace("deSolve", quietly = TRUE)) {
  stop("SOLVER = \"ode\" needs the package 'deSolve': install.packages(\"deSolve\")")
}

## --- Output files: <OUTPUT_DIR>/<experiment>_<time stamp>/<RUN_NAME>_<model>_<suffix>
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  if (!requireNamespace("writexl", quietly = TRUE))
    stop("Package 'openxlsx' is needed to write the Excel output files: install.packages(\"openxlsx\")")
  cat("NOTE: install.packages(\"openxlsx\") for Excel output with readable column widths.\n")
}
TRT_TAG   <- ""                      # set for every model in SECTION 9
MODEL_DIR <- ""                      # set for every model in SECTION 9
out_file <- function(suffix) {
  file.path(MODEL_DIR, paste0(RUN_NAME, if (nzchar(TRT_TAG)) paste0("_", TRT_TAG), "_", suffix))
}
# Excel output. The numbers are stored with FULL precision; only what Excel SHOWS
# is rounded (cell number format): EXCEL_DIGITS significant digits (SECTION 1.9),
# MLL and AIC with 1 decimal, times as they are. Click on a cell to see the exact
# value. With the package openxlsx every column also gets a width that fits its
# (shown) content, and the bold header row stays visible when scrolling.
excel_shown <- function(x, cl) {     # number format and shown text of each value
  fmt <- rep("General", length(x)); txt <- ifelse(is.na(x), "", as.character(x))
  if (cl %in% c("MLL", "AIC")) {
    fmt[] <- "0.0"; txt <- ifelse(is.na(x), "", formatC(x, format = "f", digits = 1))
  } else if (tolower(cl) != "time") {
    ok  <- is.finite(x) & x != 0
    e   <- floor(log10(abs(x[ok])))                  # order of magnitude
    d   <- pmax(0, EXCEL_DIGITS - 1 - e)             # decimals for EXCEL_DIGITS significant digits
    sci <- e < -3 | e >= 6                           # very small/large: scientific notation
    fmt[ok] <- ifelse(sci, paste0("0.", strrep("0", EXCEL_DIGITS - 1), "E+00"),
                      ifelse(d == 0, "0", paste0("0.", strrep("0", d))))
    txt[ok] <- ifelse(sci, sprintf("%.*e", as.integer(EXCEL_DIGITS - 1), x[ok]),
                      sprintf("%.*f", as.integer(d), x[ok]))
  }
  list(fmt = fmt, txt = txt)
}
write_excel <- function(sheets, path) {           # sheets: named list of data frames
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  if (!requireNamespace("openxlsx", quietly = TRUE)) {   # plain sheets, full precision
    writexl::write_xlsx(sheets, path)
    return(invisible(path))
  }
  wb  <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DDEBF7", border = "bottom",
                               wrapText = TRUE, valign = "center")
  styles <- list()                                       # one style per number format
  for (nm in names(sheets)) {
    df <- sheets[[nm]]
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, df, headerStyle = hdr)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
    w <- integer(ncol(df))
    for (j in seq_along(df)) {
      x <- df[[j]]
      if (is.numeric(x)) {
        s <- excel_shown(x, names(df)[j])
        for (f in setdiff(unique(s$fmt), "General")) {
          if (is.null(styles[[f]])) styles[[f]] <- openxlsx::createStyle(numFmt = f)
          openxlsx::addStyle(wb, nm, styles[[f]], rows = which(s$fmt == f) + 1, cols = j, stack = TRUE)
        }
        txt <- s$txt
      } else {
        txt <- as.character(x); txt[is.na(txt)] <- ""
      }
      w[j] <- max(nchar(names(df)[j]), nchar(txt))       # longest entry (header or shown value)
    }
    openxlsx::setColWidths(wb, nm, cols = seq_along(df), widths = pmin(60, pmax(8, w + 2)))
    long <- which(w > 60)                                # longer texts wrap within the column
    if (length(long) && nrow(df))
      openxlsx::addStyle(wb, nm, openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                         rows = seq_len(nrow(df)) + 1, cols = long, gridExpand = TRUE, stack = TRUE)
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

## --- Plot colours (CLR_..., NOT COL_...: those are the column names in 1.2) ---------------------------------------------------------------
CLR_PARENT   <- "#3B78B0"    # first curve: parent compound / whole organism (blue)
CLR_SECOND   <- "#6E9A1F"    # second curve (green)
CLR_THIRD    <- "#D98E2B"    # third curve (orange)
CLR_PS_INNER <- "green"      # parameter space: sets within the 95% CI of single parameters
CLR_PS_JOINT <- "cyan"       # parameter space: sets within the joint 95% region

## --- Log: everything printed with say() is also written to the log file ----
.LOG <- new.env()
.LOG$lines <- character(0)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  .LOG$lines <- c(.LOG$lines, txt)
  invisible(txt)
}
rule <- function(ch = "=") say(strrep(ch, 78))

## --- Progress bar helpers ---------------------------------------------------
pb_new <- function(n, title) {
  if (!SHOW_PROGRESS) return(NULL)
  cat("  ", title, "\n", sep = "")
  pb <- utils::txtProgressBar(min = 0, max = max(1, n), style = 3, width = 50)
  attr(pb, "nmax") <- max(1, n)
  pb
}
pb_set <- function(pb, i) {
  if (!is.null(pb)) utils::setTxtProgressBar(pb, min(i, attr(pb, "nmax")))
}
pb_close <- function(pb) {
  if (!is.null(pb)) { utils::setTxtProgressBar(pb, attr(pb, "nmax")); close(pb) }
}

## --- small utilities ----------------------------------------------------------
fmt <- function(x, d = 4) {               # nice number formatting
  ifelse(is.na(x), "NA", formatC(signif(x, d), format = "g", digits = d))
}
elapsed <- function(t0) {
  s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (s < 60) sprintf("%.1f s", s) else if (s < 3600) sprintf("%.1f min", s / 60)
  else sprintf("%.2f h", s / 3600)
}
norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))    # for name matching
sort_coll <- function(m) m[order(m[, ncol(m)]), , drop = FALSE]  # sort on MLL
rows_ab   <- function(m, a, b) {           # safe m[a:b, ] (empty if a > b)
  a <- max(1, a); b <- min(nrow(m), b)
  if (a > b) m[0, , drop = FALSE] else m[a:b, , drop = FALSE]
}


###############################################################################
## SECTION 3 : Checks of the user settings
###############################################################################

# Molar mass(es) of the compound(s) chosen in SECTION 1.4
if (!COMPOUND %in% names(MOLAR_MASSES))
  stop("No molar mass for COMPOUND = '", COMPOUND, "' in MOLAR_MASSES (SECTION 1.4).")
PARENT_NAME       <- COMPOUND
PARENT_MOLAR_MASS <- MOLAR_MASSES[[COMPOUND]]
if (!BTP_NAME %in% names(MOLAR_MASSES))
  stop("No molar mass for BTP_NAME = '", BTP_NAME, "' in MOLAR_MASSES (SECTION 1.4).")
BTP_MOLAR_MASS    <- MOLAR_MASSES[[BTP_NAME]]

check_settings <- function() {
  stopifnot(FIT_METHOD %in% c("parspace", "simplex"),
            SOLVER %in% c("analytical", "ode"),
            X0_SOURCE %in% c("data", "manual"),
            DATA_TRANSFORM >= 0, DATA_TRANSFORM <= 1,
            PS_ROUGHNESS %in% c(0, 1, 2))
  if (USE_MOLAR_UNITS && (!(PARENT_MOLAR_MASS > 0) || !(BTP_MOLAR_MASS > 0)))
    stop("Molar masses must be positive numbers (SECTION 1.4).")

  P <- PARAMETERS
  P$name  <- trimws(P$name)
  P$scale <- tolower(trimws(P$scale))
  need <- c("ku", "ke", "km", "kem", "kum")
  if (!setequal(P$name, need))
    stop("PARAMETERS must contain exactly these names: ", paste(need, collapse = ", "))
  P <- P[match(need, P$name), ]
  rownames(P) <- P$name
  if (any(!P$scale %in% c("log", "lin"))) stop("scale must be \"log\" or \"lin\".")
  if (any(!P$fit %in% c(0, 1)))           stop("fit must be 0 or 1.")

  # --- model switches overrule the parameter table ---
  if (!USE_BIOTRANSFORMATION) {
    for (nm in c("km", "kem", "kum")) { P[nm, "fit"] <- 0; P[nm, "start"] <- 0 }
  }
  if (!BTP_UPTAKE_FROM_MEDIUM || !USE_BIOTRANSFORMATION) { P["kum", "fit"] <- 0; P["kum", "start"] <- 0 }

  # --- plausibility of ranges ---
  for (nm in P$name[P$fit == 1]) {
    r <- P[nm, ]
    if (!(r$max > r$min))
      stop("Parameter ", nm, ": max must be larger than min.")
    if (r$start < r$min || r$start > r$max)
      stop("Parameter ", nm, ": start value is outside the min-max range.")
    if (r$scale == "log" && (r$min <= 0 || r$start <= 0))
      stop("Parameter ", nm, ": on log scale, min and start must be > 0.")
    if (r$min < 0) stop("Parameter ", nm, ": rate constants cannot be negative.")
  }
  if (any(P$start < 0)) stop("Parameter start values cannot be negative.")
  if (P["ke", "start"] <= 0 || (P["ke", "fit"] == 1 && P["ke", "min"] <= 0))
    stop("ke must be > 0 (also its lower bound).")
  if (sum(P$fit) > 10) stop("The parameter-space explorer can fit 10 parameters or less.")
  P
}
BIO_VERSIONS <- unique(as.logical(USE_BIOTRANSFORMATION))   # model versions to run
if (!length(BIO_VERSIONS) || any(is.na(BIO_VERSIONS)))
  stop("USE_BIOTRANSFORMATION must be TRUE, FALSE or c(TRUE, FALSE).")
BTP_UPTAKE_SETTING <- BTP_UPTAKE_FROM_MEDIUM
USE_BIOTRANSFORMATION <- BIO_VERSIONS[1]
PAR <- check_settings()               # checked for every version in setup_version()



###############################################################################
## SECTION 4 : Data import and preparation
###############################################################################

## 4.1 Select the data file (file browser if DATA_FILE is "") -------------------
choose_data_file <- function() {
  if (!interactive())
    stop("DATA_FILE is empty, but R is not running interactively.\n",
         "  Please give the path of the data file in DATA_FILE (SECTION 1.1).")
  # The dialog starts in the folder of THIS script. (Base file.choose() cannot
  # be given a start folder, so the script folder is found via RStudio or
  # source(), and passed to a dialog that accepts one.)
  start_dir <- tryCatch({
    p <- NULL
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
      p <- rstudioapi::getSourceEditorContext()$path          # script open in RStudio
    if (is.null(p) || !nzchar(p)) for (fr in rev(sys.frames()))
      if (!is.null(fr$ofile)) { p <- fr$ofile; break }        # source("...R")
    if (!is.null(p) && nzchar(p)) dirname(normalizePath(p, mustWork = FALSE)) else getwd()
  }, error = function(e) getwd())
  f <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    f <- tryCatch(rstudioapi::selectFile(caption = "Select the data file (Excel or CSV)",
                                         path = start_dir,
                                         filter = "Data files (*.xlsx *.xls *.csv)"),
                  error = function(e) NULL)
  } else if (.Platform$OS.type == "windows") {
    cat("Please select the data file in the dialog window",
        "(it may open behind other windows) ...\n")
    f <- utils::choose.files(default = file.path(start_dir, "*.*"), multi = FALSE)
    if (!length(f)) f <- ""
  } else {
    cat("Please select the data file in the dialog window",
        "(it may open behind other windows) ...\n")
    f <- file.choose()
  }
  if (is.null(f) || !nzchar(f)) stop("No data file selected.")
  f
}
if (!nzchar(DATA_FILE)) DATA_FILE <- choose_data_file()
if (!file.exists(DATA_FILE))
  stop("Data file not found: ", normalizePath(DATA_FILE, mustWork = FALSE),
       "\n  Working directory is: ", getwd())
DATA_FILE <- normalizePath(DATA_FILE)
if (grepl("\\.xlsx?$", tolower(DATA_FILE)) && !requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is needed to read Excel files.\n",
       "  Install it with: install.packages(\"readxl\")\n",
       "  or save your data as .csv.")
}
# Output folder: a relative OUTPUT_DIR is placed next to the data file
if (!grepl("^(/|~|[A-Za-z]:|\\\\\\\\)", OUTPUT_DIR)) OUTPUT_DIR <- file.path(dirname(DATA_FILE), OUTPUT_DIR)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

## 4.2 Read the file (all Excel sheets, or one sheet, or a CSV) --------------------
drop_empty_rows <- function(df) {
  if (!nrow(df)) return(df)
  filled <- as.matrix(as.data.frame(lapply(df, function(x) !is.na(x) & trimws(as.character(x)) != "")))
  df[rowSums(filled) > 0, , drop = FALSE]
}
read_tables <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("xlsx", "xls")) {
    sheets <- if (identical(EXCEL_SHEET, "all")) readxl::excel_sheets(file) else EXCEL_SHEET
    # If the table does not start in row 1 (e.g. a title cell above it, so the
    # columns are read as '...2', '...3', ...), the row that contains the
    # column names of SECTION 1.2 is searched in the first 25 rows.
    read_one <- function(s) {
      df <- as.data.frame(readxl::read_excel(file, sheet = s, guess_max = 1e6))
      if (any(norm_name(names(df)) == norm_name(COL_EXPERIMENT))) return(df)
      peek <- suppressMessages(as.data.frame(
        readxl::read_excel(file, sheet = s, col_names = FALSE, n_max = 25)))
      hdr <- which(apply(peek, 1, function(r) any(norm_name(as.character(r)) == norm_name(COL_EXPERIMENT))))
      if (!length(hdr)) return(df)
      say("  Sheet '", s, "': column names found in row ", hdr[1], "; the rows above are skipped.")
      as.data.frame(readxl::read_excel(file, sheet = s, skip = hdr[1] - 1, guess_max = 1e6))
    }
    tabs <- lapply(sheets, read_one)
    names(tabs) <- paste("sheet", sheets)
  } else if (ext %in% c("csv", "txt")) {
    sep <- CSV_SEPARATOR
    if (sep == "auto") {
      first <- readLines(file, n = 1, warn = FALSE)
      sep <- if (nchar(gsub("[^;]", "", first)) > nchar(gsub("[^,]", "", first))) ";" else ","
    }
    rd <- function(enc) read.csv(file, sep = sep, dec = CSV_DECIMAL, check.names = FALSE,
                                 stringsAsFactors = FALSE, strip.white = TRUE,
                                 fileEncoding = enc)
    tabs <- list(csv = tryCatch(rd("UTF-8-BOM"), warning = function(w) rd("latin1")))
  } else {
    stop("Unknown file type '.", ext, "'. Use .xlsx, .xls or .csv")
  }
  lapply(tabs, drop_empty_rows)
}

## 4.3 Find a column by (tolerant) name -------------------------------------
find_col <- function(df, wanted, required = TRUE) {
  nms <- names(df)
  i <- which(nms == wanted)                                    # 1) exact
  if (!length(i)) i <- which(norm_name(nms) == norm_name(wanted))   # 2) normalised
  if (!length(i)) {                                            # 3) part before "["
    stem <- norm_name(sub("\\[.*$", "", wanted))
    i <- which(norm_name(sub("\\[.*$", "", nms)) == stem)
  }
  if (length(i) > 1) i <- i[1]
  if (!length(i)) {
    if (!required) return(NA_character_)
    stop("Column '", wanted, "' not found in the data file.\n  Available columns: ",
         paste0("'", nms, "'", collapse = ", "), "\n  Adjust SECTION 1.2.")
  }
  nms[i]
}

## 4.4 Text to numbers: "0,42" -> 0.42; "<LOQ" etc. -> BELOW_LOQ_VALUE ------------
to_num <- function(x, what) {
  if (is.numeric(x)) return(as.numeric(x))
  x    <- trimws(as.character(x))
  y    <- suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
  txt  <- !is.na(x) & x != "" & is.na(y)                     # text entries
  miss <- txt & tolower(x) %in% tolower(MISSING_CODES)       # -> missing
  loq  <- txt & !miss                                        # -> below LOQ
  y[loq] <- BELOW_LOQ_VALUE
  if (any(loq))
    say("  NOTE: ", sum(loq), " entries in '", what, "' set to ", BELOW_LOQ_VALUE,
        " (below LOQ): ", paste(head(unique(x[loq]), 6), collapse = ", "))
  if (any(miss))
    say("  NOTE: ", sum(miss), " entries in '", what, "' treated as missing: ",
        paste(head(unique(x[miss]), 6), collapse = ", "))
  y
}

## 4.5 Bring each table into the standard format and combine ---------------------
# Numbers in text columns (e.g. treatment 20.1 next to "AB") are written in
# a standard way, so "20.1", "20.10" and 20.1 are the same treatment.
clean_label <- function(x) {
  x <- trimws(as.character(x))
  num <- suppressWarnings(as.numeric(x))
  x[!is.na(num)] <- as.character(num[!is.na(num)])
  x
}
standardise <- function(df, src) {
  cols  <- c(COL_EXPERIMENT, COL_MATRIX, COL_SPECIES, COL_TREATMENT, COL_TIME, COL_REPLICATE, COL_PARENT)
  found <- unname(sapply(cols, function(cl) find_col(df, cl, required = FALSE)))
  if (any(is.na(found))) {
    say("  ", src, ": skipped (columns not found: ", paste(cols[is.na(found)], collapse = ", "), ")")
    return(NULL)
  }
  cb <- find_col(df, COL_BTP, required = FALSE)
  if (is.na(cb) && (any(BIO_VERSIONS) || BTP_UPTAKE_FROM_MEDIUM))
    stop("Column '", COL_BTP, "' not found in ", src, ". Available: ",
         paste0("'", names(df), "'", collapse = ", "),
         ". Adjust COL_BTP or set USE_BIOTRANSFORMATION <- FALSE.")
  cu <- find_col(df, COL_UNIT, required = FALSE)
  cc <- find_col(df, COL_COMPOUND, required = FALSE)
  say("  ", src, ": ", nrow(df), " rows")
  data.frame(
    code      = clean_label(df[[found[1]]]),
    matrix    = trimws(as.character(df[[found[2]]])),
    species   = trimws(as.character(df[[found[3]]])),
    treatment = clean_label(df[[found[4]]]),
    time      = to_num(df[[found[5]]], paste(src, found[5])),
    replicate = trimws(as.character(df[[found[6]]])),
    compound  = if (is.na(cc)) NA_character_ else trimws(as.character(df[[cc]])),
    parent    = to_num(df[[found[7]]], paste(src, found[7])),
    btp       = if (is.na(cb)) NA_real_ else to_num(df[[cb]], paste(src, cb)),
    unit      = if (is.na(cu)) NA_character_ else trimws(as.character(df[[cu]])),
    row_order = seq_len(nrow(df)),                 # file order (for step changes)
    stringsAsFactors = FALSE)
}

rule(); say("TK model with biotransformation (BYOM TK-metabolism) - TK Workshop, Exercise 2"); rule()
say("Run name : ", RUN_NAME, "   (", format(Sys.time(), "%d-%b-%Y %H:%M"), ")")
say("Data file: ", DATA_FILE)
say("Output   : ", normalizePath(OUTPUT_DIR))

TABS <- read_tables(DATA_FILE)
DAT  <- do.call(rbind, Map(standardise, TABS, names(TABS)))
if (is.null(DAT) || !nrow(DAT)) {
  for (cl in c(COL_EXPERIMENT, COL_MATRIX, COL_SPECIES, COL_TREATMENT, COL_TIME, COL_REPLICATE, COL_PARENT))
    find_col(TABS[[1]], cl)                     # gives an informative error
  stop("No usable data found in the file.")
}
# Keep only the rows of the compound chosen in SECTION 1.4 (COMPOUND)
if (any(!is.na(DAT$compound))) {
  comps <- unique(DAT$compound[!is.na(DAT$compound)])
  say("  Compounds in the file: ", paste(comps, collapse = ", "))
  sel <- tolower(DAT$compound) %in% tolower(COMPOUND)
  if (!any(sel))
    stop("COMPOUND '", COMPOUND, "' (SECTION 1.4) is not in the Compound column. Available: ",
         paste(comps, collapse = ", "))
  DAT <- DAT[sel, ]
  say("  Selected compound: ", COMPOUND, " (", nrow(DAT), " rows)")
}
DAT$row_order <- seq_len(nrow(DAT))
DAT$replicate[is.na(DAT$replicate) | DAT$replicate %in% c("", "NA")] <- "1"
if (any(is.na(DAT$time))) {
  say("  NOTE: ", sum(is.na(DAT$time)), " rows without a time value are removed.")
  DAT <- DAT[!is.na(DAT$time), ]
}
if (any(DAT$time < 0)) {       # samples before t = 0 (e.g. leaves soaking before the animals were added)
  neg <- DAT$time < 0 & tolower(DAT$matrix) == tolower(ORGANISM_MATRIX)
  say("  NOTE: ", sum(DAT$time < 0), " samples before t = 0 are only used to interpolate the exposure at t = 0",
      if (any(neg)) paste0(" (", sum(neg), " tissue samples before t = 0 removed)"), ".")
  DAT <- DAT[!neg, ]
}
if (any(DAT$parent < 0, na.rm = TRUE) || any(DAT$btp < 0, na.rm = TRUE))
  stop("Negative concentrations found in the data. Please check.")

## 4.6 Matrices and unit conversion (ug -> umol) ---------------------------------
is_exp <- tolower(DAT$matrix) == tolower(EXPOSURE_MATRIX)
is_org <- tolower(DAT$matrix) == tolower(ORGANISM_MATRIX)
if (!any(is_exp)) stop("EXPOSURE_MATRIX '", EXPOSURE_MATRIX, "' not found. Matrices in the file: ",
                       paste(unique(DAT$matrix), collapse = ", "))
if (!any(is_org)) stop("ORGANISM_MATRIX '", ORGANISM_MATRIX, "' not found. Matrices in the file: ",
                       paste(unique(DAT$matrix), collapse = ", "))
ign <- setdiff(unique(DAT$matrix[!is_exp & !is_org]), NA)
if (length(ign)) say("  Matrices not used: ", paste(ign, collapse = ", "))
DAT <- DAT[is_exp | is_org, ]
DAT$is_exp <- tolower(DAT$matrix) == tolower(EXPOSURE_MATRIX)

units_found <- function(sel) paste(setdiff(unique(DAT$unit[sel]), NA), collapse = ", ")
say("  Units in the file: medium ", units_found(DAT$is_exp), " | tissue ", units_found(!DAT$is_exp))
F_PAR  <- if (USE_MOLAR_UNITS) 1 / PARENT_MOLAR_MASS else 1
F_BTP  <- if (USE_MOLAR_UNITS) 1 / BTP_MOLAR_MASS    else 1
UNIT     <- if (USE_MOLAR_UNITS) "umol/kg ww" else DATA_UNIT_TISSUE   # tissue
UNIT_MED <- if (USE_MOLAR_UNITS) "umol/L" else DATA_UNIT_MEDIUM       # medium
DAT$parent_m <- DAT$parent * F_PAR       # model units
DAT$btp_m    <- DAT$btp    * F_BTP

## 4.7 Scenarios = experiment code x treatment (x species) --------------------------
# The species is added to the name only where one code and treatment has
# several species (e.g. AE5 "TR": A. aquaticus and G. pulex).
key3 <- paste(DAT$code, DAT$treatment, DAT$species, sep = " | ")
key2 <- paste(DAT$code, DAT$treatment, sep = " | ")
n_sp <- tapply(DAT$species, key2, function(s) length(unique(s)))
DAT$scenario <- ifelse(n_sp[key2] > 1, key3, key2)

CATALOG <- do.call(rbind, lapply(split(DAT, DAT$scenario), function(d) data.frame(
  scenario = d$scenario[1], code = d$code[1], treatment = d$treatment[1], species = d$species[1],
  n_medium = sum(d$is_exp), n_tissue = sum(!d$is_exp),
  tissue_times = paste(signif(sort(unique(d$time[!d$is_exp])), 4), collapse = ", "),
  stringsAsFactors = FALSE)))
CATALOG <- CATALOG[order(CATALOG$code, CATALOG$species, CATALOG$treatment), ]
rownames(CATALOG) <- NULL
TRT_TAG <- ""
say("")
say("Available scenarios (\"code | treatment\"; medium / tissue rows):")
for (i in seq_len(nrow(CATALOG)))
  say(sprintf("  %-34s %-14s %3d / %3d   tissue at t = %s", paste0("\"", CATALOG$scenario[i], "\""),
              CATALOG$species[i], CATALOG$n_medium[i], CATALOG$n_tissue[i], CATALOG$tissue_times[i]))

## 4.8 Models: separate scenarios and pooled scenario groups ----------------------
# expand_spec: "AE2" -> all AE2 scenarios; "AE2 | 15.6" -> one scenario;
# "AE5 | AB" -> all species with code AE5 and treatment AB
expand_spec <- function(spec, allow_multi_species = TRUE) {
  parts <- trimws(strsplit(spec, "|", fixed = TRUE)[[1]])
  parts[1:min(2, length(parts))] <- clean_label(parts[1:min(2, length(parts))])
  hit <- CATALOG[tolower(CATALOG$code) == tolower(parts[1]), ]
  if (length(parts) >= 2) hit <- hit[tolower(hit$treatment) == tolower(parts[2]), ]
  if (length(parts) >= 3) hit <- hit[norm_name(hit$species) == norm_name(parts[3]), ]
  if (nrow(hit) == 0) {
    opts <- CATALOG$scenario[tolower(CATALOG$code) == tolower(parts[1])]
    stop("Scenario '", spec, "' not found in the data.",
         if (length(opts)) paste0("\n  Available for ", parts[1], ": ", paste0("\"", opts, "\"", collapse = ", "))
         else paste0("\n  Experiment codes in the file: ", paste(unique(CATALOG$code), collapse = ", ")))
  }
  if (!allow_multi_species && length(parts) == 2 && nrow(hit) > 1)
    stop("Pooled model: scenario '", spec, "' exists for several species. Add the species, e.g.: ",
         paste0("\"", hit$scenario, "\"", collapse = " or "))
  hit$scenario[hit$n_tissue > 0]
}
pick_scenarios <- function(title) {
  if (!interactive()) stop("\"choose\" needs an interactive R session.")
  sel <- utils::select.list(CATALOG$scenario[CATALOG$n_tissue > 0], multiple = TRUE,
                            graphics = TRUE, title = title)
  if (!length(sel)) stop("No scenario selected.")
  sel
}
RUNS <- list()
# separate models: one model per scenario
if (identical(SEPARATE_MODELS, "choose")) {
  sep <- pick_scenarios("Scenarios: each one a SEPARATE model")
} else if (identical(SEPARATE_MODELS, "all")) {
  sep <- CATALOG$scenario[CATALOG$n_tissue > 0]
} else {
  sep <- unique(unlist(lapply(SEPARATE_MODELS, expand_spec)))
}
for (sc in sep) RUNS[[sc]] <- sc
# pooled models: one model per list entry
if (identical(POOLED_MODELS, "choose")) {
  sel <- pick_scenarios("Scenarios to POOL in one model")
  POOLED_MODELS <- list(sel)
  names(POOLED_MODELS) <- paste0(paste(unique(sub(" \\|.*$", "", sel)), collapse = "_"), "_pooled")
}
if (length(POOLED_MODELS)) {
  if (!is.list(POOLED_MODELS)) POOLED_MODELS <- list(POOLED_MODELS)
  nm <- names(POOLED_MODELS)
  if (is.null(nm)) nm <- rep("", length(POOLED_MODELS))
  for (i in seq_along(POOLED_MODELS)) {
    v <- unique(unlist(lapply(POOLED_MODELS[[i]], expand_spec, allow_multi_species = FALSE)))
    key <- if (nzchar(nm[i])) nm[i] else paste(v, collapse = " + ")
    if (key %in% names(RUNS)) key <- paste0(key, "_pooled")
    RUNS[[key]] <- v
  }
}
if (!length(RUNS)) stop("No models selected: set SEPARATE_MODELS and/or POOLED_MODELS (SECTION 1.3).")
# every model in every requested version (with/without biotransformation)
JOBS <- expand.grid(run = names(RUNS), bio = BIO_VERSIONS, stringsAsFactors = FALSE)
JOBS <- JOBS[order(match(JOBS$run, names(RUNS)), !JOBS$bio), ]

say("")
say("Models to calibrate (", length(RUNS), " model(s) x ", length(BIO_VERSIONS), " version(s) = ",
    nrow(JOBS), " fits):")
for (i in seq_along(RUNS)) say(sprintf("  %2d. %-28s %s%s", i, names(RUNS)[i], paste(RUNS[[i]], collapse = "  +  "),
                                       if (length(RUNS[[i]]) > 1) "   (pooled)" else ""))
say("  Biotransformation: ", paste(ifelse(BIO_VERSIONS, "ON (parent + BTP fitted)", "OFF (parent only)"),
                                   collapse = "  and  "))
say("  Units            : ", if (USE_MOLAR_UNITS)
  sprintf("converted to umol (M %s = %.2f, M %s = %.2f g/mol)",
          PARENT_NAME, PARENT_MOLAR_MASS, BTP_NAME, BTP_MOLAR_MASS) else "as in the file")

## 4.9 Helpers used per scenario (SECTION 9) -------------------------------------
STEP_DT <- 1e-6      # [d] repeated time points in the medium = step change
# Exposure profile: mean of the replicates per time point. If a time point
# occurs twice for the same replicate (e.g. 1.15 then 0 at t = 2.00125), the
# second value starts STEP_DT later (step change, in file order).
make_profile <- function(d, col) {
  d <- d[!is.na(d[[col]]), c("time", "replicate", "row_order", col)]
  if (!nrow(d)) return(NULL)
  d <- d[order(d$row_order), ]
  d$k <- ave(seq_len(nrow(d)), d$time, d$replicate, FUN = seq_along)  # 1st, 2nd, ... value
  a <- aggregate(d[[col]], by = list(time = d$time, k = d$k), FUN = mean)
  names(a) <- c("time", "k", "conc")
  a <- a[order(a$time, a$k), ]
  a$time <- a$time + (a$k - 1) * STEP_DT
  if (any(a$time < 0)) {                     # samples before t = 0: only the value AT t = 0 is used
    c0 <- if (nrow(a) > 1) approx(a$time, a$conc, xout = 0, rule = 2)$y else a$conc[1]
    a  <- a[a$time >= 0, ]
    if (!nrow(a) || a$time[1] > 0) a <- rbind(data.frame(time = 0, k = 1, conc = c0), a)
  }
  a[, c("time", "conc")]
}
# Tissue observations of one state (uses ORG of the current run)
make_obs <- function(state) {
  col <- paste0(state, "_m")
  o <- ORG[!is.na(ORG[[col]]), c("scenario", "replicate", "time", col)]
  names(o)[4] <- "value"
  o
}
# Mean tissue concentration at t = 0 (start value of the model)
mean_t0 <- function(o) {
  v <- o$value[o$time == 0]
  if (length(v)) mean(v) else 0
}
# End of exposure: last medium time before the parent drops below 5% of its
# maximum for good; otherwise the last tissue sampling time
end_of_exposure <- function(sc, pp, t_tissue) {
  if (!is.null(END_OF_EXPOSURE) && sc %in% names(END_OF_EXPOSURE)) return(unname(END_OF_EXPOSURE[sc]))
  # 1) the last sudden drop of the medium concentration (transfer to clean medium)
  if (nrow(pp) > 1) {
    drop <- which(diff(pp$time) <= 2 * STEP_DT & pp$conc[-1] < 0.2 * pp$conc[-nrow(pp)])
    if (length(drop)) return(pp$time[max(drop)])
  }
  # 2) otherwise: the last time point before the medium stays below 5% of its maximum
  low <- pp$conc < 0.05 * max(pp$conc)
  n <- length(low)
  first_low <- NA
  for (i in seq_len(n)) if (all(low[i:n])) { first_low <- i; break }
  if (!is.na(first_low) && first_low > 1) return(pp$time[first_low - 1])
  max(t_tissue)
}
# All inputs of one scenario; returns list(skip = reason) if it cannot be used
setup_scenario <- function(sc) {
  d <- DAT[DAT$scenario == sc, ]
  e <- d[d$is_exp, ]
  o <- d[!d$is_exp, ]
  pp <- make_profile(e, "parent_m")
  pm <- make_profile(e, "btp_m")
  if (is.null(pp)) return(list(skip = paste("no", EXPOSURE_MATRIX, "data")))
  if (all(pp$conc == 0)) return(list(skip = "no parent in the medium (e.g. control)"))
  if (!nrow(o)) return(list(skip = paste("no", ORGANISM_MATRIX, "data")))
  if (pp$time[1] > 0)
    say("    NOTE: no medium value at t = 0; the first value (t = ", pp$time[1], " d) is used from t = 0 on.")
  if (BTP_UPTAKE_FROM_MEDIUM && is.null(pm))
    stop("BTP_UPTAKE_FROM_MEDIUM = TRUE, but no BTP data in the medium for '", sc, "'.")
  list(skip = NULL, expo = list(parent = pp, btp = pm), org = o, exp = e)
}



###############################################################################
## SECTION 5 : The model
###############################################################################

## 5.1 THE MODEL AS ODEs  (equivalent of BYOM derivatives.m) ------------------
# This function documents the model and is used when SOLVER = "ode" and for
# the solver self-check. If you change the model, change it here and run with
# SOLVER = "ode" (the analytical solution below is only valid for THIS model).
tk_derivatives <- function(t, y, parms, cw_p, cw_m) {
  Cp <- y[1]                       # internal parent concentration
  Cm <- y[2]                       # internal BTP concentration
  c_p <- cw_p(t)                   # parent in medium at time t
  c_m <- cw_m(t)                   # BTP in medium at time t
  with(as.list(parms), {
    dCp <- ku * c_p - ke * Cp - km * Cp          # uptake - elimination - biotransformation
    dCm <- km * Cp - kem * Cm + kum * c_m        # formation - elimination (+ uptake)
    list(c(dCp, dCm))
  })
}

## 5.2 Pre-computation of the time grid for one scenario ---------------------------
# The exposure is piecewise linear between "breakpoints" (all measurement times
# of parent and BTP in the medium, plus t = 0). For each interval we store
# the concentration at its start (a) and its slope (b). x0 = start values.
forcing_at <- function(prof, t) {
  if (is.null(prof)) return(rep(0, length(t)))
  if (nrow(prof) == 1) return(rep(prof$conc, length(t)))
  approx(prof$time, prof$conc, xout = t, rule = 2)$y   # rule 2: constant outside
}
make_grid <- function(expo, t_out, x0) {
  bp <- sort(unique(c(0, expo$parent$time, if (!is.null(expo$btp)) expo$btp$time)))
  a  <- forcing_at(expo$parent, bp)
  am <- forcing_at(expo$btp, bp)
  nb <- length(bp)
  dt <- diff(bp)
  b  <- if (nb > 1) c(diff(a)  / dt, 0) else 0   # slope; 0 after the last point
  bm <- if (nb > 1) c(diff(am) / dt, 0) else 0
  idx <- findInterval(t_out, bp)                  # interval of each output time
  list(bp = bp, dt = dt, a = a, b = b, am = am, bm = bm,
       idx = idx, s = t_out - bp[idx], t_out = t_out, x0 = x0)
}

## 5.3 Exact (analytical) solution  ------------------------------------------------
# Helper functions, numerically stable for k -> 0:
#   phi(k,s) = integral_0^s exp(-k(s-u)) du     = (1 - exp(-k s)) / k
#   psi(k,s) = integral_0^s exp(-k(s-u)) u du   = (k s - 1 + exp(-k s)) / k^2
phi <- function(k, s) {
  x <- k * s
  ifelse(abs(x) < 1e-6, s * (1 - x / 2 + x * x / 6), -expm1(-x) / k)
}
psi <- function(k, s) {
  x <- k * s
  ifelse(abs(x) < 1e-4, s * s * (0.5 - x / 6 + x * x / 24), (x + expm1(-x)) / (k * k))
}
# Solution over a time span s within one interval, starting from Cp0, Cm0 with
# medium concentration Cw_p = a + b*s (and Cw_m = am + bm*s):
#   Cp(s) = Cp0 e^(-k1 s) + ku [a phi(k1,s) + b psi(k1,s)]          k1 = ke + km
#   Cm(s) = Cm0 e^(-k2 s) + km * convolution(e^(-k2 .), Cp)           k2 = kem
#           + kum [am phi(k2,s) + bm psi(k2,s)]
tk_step <- function(s, Cp0, Cm0, a, b, am, bm, ku, k1, km, k2, kum) {
  Cp <- Cp0 * exp(-k1 * s) + ku * (a * phi(k1, s) + b * psi(k1, s))
  A  <- ku * a / k1 - ku * b / (k1 * k1)          # Cp = A + B s + D exp(-k1 s)
  B  <- ku * b / k1
  D  <- Cp0 - A
  kmin <- min(k1, k2)
  e12  <- exp(-kmin * s) * phi(abs(k2 - k1), s)    # (e^-k1s - e^-k2s)/(k2 - k1), stable
  Cm <- Cm0 * exp(-k2 * s) + km * (A * phi(k2, s) + B * psi(k2, s) + D * e12) +
        kum * (am * phi(k2, s) + bm * psi(k2, s))
  list(Cp = Cp, Cm = Cm)
}
solve_analytical <- function(p, G) {
  ku <- p[["ku"]]; km <- p[["km"]]; k1 <- p[["ke"]] + km; k2 <- p[["kem"]]; kum <- p[["kum"]]
  nb <- length(G$bp)
  Cp <- numeric(nb); Cm <- numeric(nb)
  Cp[1] <- G$x0[[1]]; Cm[1] <- G$x0[[2]]
  if (nb > 1) for (j in 1:(nb - 1)) {                 # states at the breakpoints
    r <- tk_step(G$dt[j], Cp[j], Cm[j], G$a[j], G$b[j], G$am[j], G$bm[j], ku, k1, km, k2, kum)
    Cp[j + 1] <- r$Cp; Cm[j + 1] <- r$Cm
  }
  i <- G$idx                                           # all output times at once
  r <- tk_step(G$s, Cp[i], Cm[i], G$a[i], G$b[i], G$am[i], G$bm[i], ku, k1, km, k2, kum)
  cbind(parent = r$Cp, btp = r$Cm)
}

## 5.4 Numerical solution with deSolve (optional) ------------------------------------
solve_ode <- function(p, G) {
  mk <- function(v) if (length(G$bp) == 1) function(t) v[1] else
    approxfun(G$bp, v, rule = 2)              # same piecewise-linear exposure
  cw_p <- mk(G$a); cw_m <- mk(G$am)
  nb  <- length(G$bp)
  y   <- c(Cp = G$x0[[1]], Cm = G$x0[[2]])
  res <- matrix(NA_real_, length(G$t_out), 2, dimnames = list(NULL, c("parent", "btp")))
  # The solver is restarted at every exposure breakpoint, so it never steps
  # over a (sudden) change in the exposure concentration.
  for (j in seq_len(nb)) {
    t0  <- G$bp[j]
    t1  <- if (j < nb) G$bp[j + 1] else max(c(G$t_out, t0))
    sel <- which(G$idx == j)
    tt  <- sort(unique(c(t0, G$t_out[sel], t1)))
    if (length(tt) > 1) {
      out <- deSolve::ode(y = y, times = tt, func = tk_derivatives, parms = p,
                          method = "lsoda", rtol = 1e-10, atol = 1e-12,
                          cw_p = cw_p, cw_m = cw_m)
      if (length(sel)) res[sel, ] <- out[match(G$t_out[sel], out[, "time"]), c("Cp", "Cm")]
      y <- out[nrow(out), c("Cp", "Cm")]
    } else if (length(sel)) {
      res[sel, ] <- matrix(y, nrow = length(sel), ncol = 2, byrow = TRUE)
    }
  }
  res
}
solve_tk <- function(p, G) {
  if (SOLVER == "ode") solve_ode(p, G) else solve_analytical(p, G)
}

## 5.5 Model prediction at the data time points of the current run ------------------
# (GRID_DATA and SCEN_IDS are set per model run in SECTION 9)
predict_data <- function(p) do.call(rbind, lapply(SCEN_IDS, function(sc) solve_tk(p, GRID_DATA[[sc]])))


###############################################################################
## SECTION 6 : Likelihood  (equivalent of BYOM transfer.m)
###############################################################################

LAM <- DATA_TRANSFORM
transform_data <- function(x) if (LAM == 0) log(pmax(x, 1e-10)) else x^LAM
back_transform <- function(y) if (LAM == 0) exp(y) else pmax(y, 0)^(1 / LAM)

calc_mll <- function(p, details = FALSE) {
  PRED <- predict_data(p)
  if (any(!is.finite(PRED))) return(Inf)
  mll <- 0
  gof <- c(); sig <- c()
  for (st in FIT_STATES) {
    o   <- OBS[[st]]
    D   <- transform_data(o$value)
    M   <- transform_data(pmax(PRED[o$row, st], 0))
    res <- D - M
    n   <- length(res)
    ssq <- sum(res^2)
    mll <- mll + (n / 2) * log(ssq) + n / 2         # BYOM: -(n/2)log(SSQ) - N*SSQ/(2*SSQ)
    if (details) {
      gof[st] <- 1 - ssq / sum((D - mean(D))^2)     # R-square (transformed, as BYOM)
      sig[st] <- sqrt(ssq / n)                      # residual SD (transformed scale, ML)
    }
  }
  if (!is.finite(mll)) mll <- Inf
  if (details) return(list(mll = mll, r2 = gof, sigma = sig, pred = PRED))
  mll
}

## 6.1 Model version: parameter vector <-> fitted vector ---------------------------
# One curve in the tissue panel of the fit plot: column of the model output,
# tissue data shown with it (NA = none), legend label, colour, line type,
# 95% CI band (TRUE/FALSE), factor from model units to the units of the data file.
curve_spec <- function(col, obs, label, colour, lty = 1, band = TRUE, f = 1 / F_PAR)
  list(col = col, obs = obs, label = label, colour = colour, lty = lty, band = band, f = f)
# setup_version() sets everything that depends on the biotransformation switch
# (fitted parameters, fitted states, labels). Called for every fit in SECTION 9.
setup_version <- function(bio) {
  USE_BIOTRANSFORMATION  <<- bio
  BTP_UPTAKE_FROM_MEDIUM <<- BTP_UPTAKE_SETTING && bio   # BTP uptake needs the BTP model
  PAR        <<- check_settings()
  FIT_STATES <<- if (bio) c("parent", "btp") else "parent"
  FIT_NAMES  <<- PAR$name[PAR$fit == 1]
  N_FIT      <<- length(FIT_NAMES)
  IS_LOG     <<- PAR[FIT_NAMES, "scale"] == "log"
  P_FIXED    <<- setNames(PAR$start, PAR$name)             # natural scale, all parameters
  b <- cbind(PAR[FIT_NAMES, "min"], PAR[FIT_NAMES, "max"])  # bounds on fit scale
  b[IS_LOG, ] <- log10(b[IS_LOG, , drop = FALSE])
  rownames(b) <- FIT_NAMES
  BNDS       <<- b
  L_COL      <<- N_FIT + 1
  MODEL_LABEL <<- if (bio) "TK model with biotransformation" else "one-compartment TK model (without biotransformation)"
  CURVE_SPEC  <<- c(list(curve_spec("parent", "parent", if (bio) paste0(PARENT_NAME, " (parent)") else PARENT_NAME, CLR_PARENT)),
                    if (bio) list(curve_spec("btp", "btp", paste0(BTP_NAME, " (BTP)"), CLR_SECOND, f = 1 / F_BTP)))
  DESCR <<- c(ku  = "uptake rate constant (parent)",
              ke  = if (bio) "elimination rate constant (parent)"
                    else "total elimination rate constant (parent, incl. biotransformation)",
              km  = "biotransformation rate constant (parent -> BTP)",
              kem = "elimination rate constant (BTP)",
              kum = "uptake rate constant of BTP from medium")
  DER_NOTE <<- c(
    BCF_parent      = if (bio) "ku/(ke+km): kinetic BCF of the parent [L/kg ww]"
                      else "ku/ke: kinetic BCF of the parent [L/kg ww]",
    t_half_parent_d = "ln2/(ke+km): elimination half-life of the parent [d]",
    t95_parent_d    = "ln20/(ke+km): time to 95% of steady state (parent) [d]",
    BCF_total_molar = "ku/(ke+km)*(1+km/kem): steady-state (parent+BTP, molar) / parent in medium [L/kg ww]",
    fraction_biotransformed = "km/(ke+km): fraction of the parent eliminated by biotransformation [-]",
    t_half_BTP_d    = "ln2/kem: elimination half-life of the BTP [d]",
    BTP_to_parent_ratio_ss  = "km/kem: steady-state molar ratio BTP/parent in the organism [-]")
  invisible(NULL)
}
setup_version(BIO_VERSIONS[1])

to_fit_scale <- function(v) { v[IS_LOG] <- log10(v[IS_LOG]); v }
to_natural   <- function(pf) {
  p <- P_FIXED
  v <- pf; v[IS_LOG] <- 10^v[IS_LOG]
  p[FIT_NAMES] <- v
  p
}
# Objective function: minus log-likelihood for the vector of fitted parameters
# (on fit scale). Outside the min-max bounds: Inf (as BYOM).
OBJ <- function(pf) {
  if (any(!is.finite(pf)) || any(pf < BNDS[, 1] | pf > BNDS[, 2])) return(Inf)
  calc_mll(to_natural(pf))
}



###############################################################################
## SECTION 7 : Optimisation algorithms
###############################################################################

## 7.1 Nelder-Mead simplex (clone of Matlab's fminsearch, as used by BYOM) ---------
eps_of <- function(x) .Machine$double.eps * max(abs(x), 1e-300)
fminsearch_nm <- function(fn, x0, tol_x = 1e-4, tol_f = 1e-4,
                          max_fev = 200 * length(x0), max_iter = 200 * length(x0),
                          on_eval = NULL) {
  n <- length(x0)
  if (n == 0) return(list(par = x0, value = fn(x0), fevals = 1, converged = TRUE))
  fcount <- 0
  f <- function(x) { fcount <<- fcount + 1; if (!is.null(on_eval)) on_eval(fcount); fn(x) }
  V  <- matrix(x0, nrow = n, ncol = n + 1)
  fv <- numeric(n + 1)
  fv[1] <- f(x0)
  for (j in seq_len(n)) {                     # initial simplex: +5% (or 0.00025)
    y <- x0
    y[j] <- if (y[j] != 0) 1.05 * y[j] else 0.00025
    V[, j + 1] <- y
    fv[j + 1]  <- f(y)
  }
  o <- order(fv); fv <- fv[o]; V <- V[, o, drop = FALSE]
  iter <- 1; converged <- FALSE
  while (fcount < max_fev && iter < max_iter) {
    dfv <- max(abs(fv[1] - fv[-1]))
    dx  <- max(abs(V[, -1, drop = FALSE] - V[, 1]))
    if (is.finite(dfv) && dfv <= max(tol_f, 10 * eps_of(fv[1])) &&
        dx <= max(tol_x, 10 * eps_of(max(V[, 1])))) { converged <- TRUE; break }
    xbar <- rowMeans(V[, 1:n, drop = FALSE])
    xr   <- 2 * xbar - V[, n + 1]                       # reflection
    fxr  <- f(xr)
    if (fxr < fv[1]) {
      xe  <- 3 * xbar - 2 * V[, n + 1]                  # expansion
      fxe <- f(xe)
      if (fxe < fxr) { V[, n + 1] <- xe; fv[n + 1] <- fxe } else { V[, n + 1] <- xr; fv[n + 1] <- fxr }
    } else if (fxr < fv[n]) {
      V[, n + 1] <- xr; fv[n + 1] <- fxr
    } else {
      shrink <- FALSE
      if (fxr < fv[n + 1]) {
        xc <- 1.5 * xbar - 0.5 * V[, n + 1]             # outside contraction
        fxc <- f(xc)
        if (fxc <= fxr) { V[, n + 1] <- xc; fv[n + 1] <- fxc } else shrink <- TRUE
      } else {
        xcc <- 0.5 * xbar + 0.5 * V[, n + 1]            # inside contraction
        fxcc <- f(xcc)
        if (fxcc < fv[n + 1]) { V[, n + 1] <- xcc; fv[n + 1] <- fxcc } else shrink <- TRUE
      }
      if (shrink) for (j in 2:(n + 1)) {                # shrink towards the best
        V[, j] <- V[, 1] + 0.5 * (V[, j] - V[, 1])
        fv[j]  <- f(V[, j])
      }
    }
    o <- order(fv); fv <- fv[o]; V <- V[, o, drop = FALSE]
    iter <- iter + 1
  }
  list(par = V[, 1], value = fv[1], fevals = fcount, converged = converged)
}
# BYOM setup_simplex: rough (TolX = TolFun = 1e-2, 30*n evaluations) or normal
run_simplex <- function(x0, rough, fn = OBJ, on_eval = NULL) {
  if (rough) fminsearch_nm(fn, x0, 1e-2, 1e-2, 30 * length(x0), on_eval = on_eval)
  else       fminsearch_nm(fn, x0, on_eval = on_eval)
}
# Simplex with parameter i_fix kept at its value in p (for profile likelihoods)
simplex_fixed <- function(p, i_fix, rough) {
  sub <- function(x) { q <- p; q[-i_fix] <- x; OBJ(q) }
  r <- run_simplex(p[-i_fix], rough, sub)
  q <- p; q[-i_fix] <- r$par
  list(par = q, value = r$value)
}

## 7.2 Settings of the parameter-space explorer (BYOM setup_settings.m) ----------
ps_settings <- function(rough) {
  S <- list()
  S$crit_table <- pmin(11.07, c(3.8415, 5.9915, 7.8147, 9.4877, 11.07, 12.592,
                                14.0671, 15.5073, 16.9190, 18.3070))  # chi2 95%, df 1..10
  S$crit_prop  <- c(-2 * 0.3764, 0.3764) + S$crit_table[1]
  S$crit_add   <- c(15, 7, 3, 2, 1.5, 1, 0.5, 0.25)
  S$n_tr       <- c(NA, 60, 40, 30, 20, 10, 10, 10)
  S$f_d        <- c(NA, 1, 0.65, 0.42, 0.27, 0.18, 0.12, 0.08)
  S$n_ok       <- c(5, 50, 200, 200, 400, 800, 800, 800, 800, 800)
  S$n_conf_all <- rbind(c(500, 500), c(2500, 1500), c(5000, 2000), c(7000, 3000),
                        c(10000, 5000), c(15000, 7500), c(15000, 7500), c(15000, 7500),
                        c(15000, 7500), c(15000, 7500))
  S$tries <- 10
  if (rough >= 1) {
    S$n_conf_all <- rbind(c(100, 100), c(250, 150), c(500, 200), c(700, 300),
                          c(1000, 500), c(1500, 750), c(1500, 750), c(1500, 750),
                          c(1500, 750), c(1500, 750))
    S$tries <- 8
  }
  S$crit_add_extra <- 3
  S$f_d_extra      <- 1
  S$d_extra        <- 0.1
  S$gap_extra      <- if (rough >= 1) 0.5 else 0.25
  S$real_better    <- 0.05
  S$n_max          <- if (rough > 1) 4 else 12
  S
}

## 7.3 Random mutations around accepted sets (BYOM rand_mutations.m) --------------
ps_mutations <- function(coll_ok, n_tr_i, d_grid_i, label) {
  n_cont <- nrow(coll_ok)
  out <- matrix(Inf, n_cont * n_tr_i, N_FIT + 1)
  pb <- pb_new(n_cont * n_tr_i, label); k <- 0
  for (i_ok in seq_len(n_cont)) {
    p_try <- matrix(0, n_tr_i, N_FIT)
    for (i_p in seq_len(N_FIT)) {
      v <- coll_ok[i_ok, i_p] + (runif(n_tr_i) * 2 - 1) * d_grid_i[i_p]
      p_try[, i_p] <- pmin(pmax(v, BNDS[i_p, 1]), BNDS[i_p, 2])
    }
    for (i_t in seq_len(n_tr_i)) {
      k <- k + 1
      out[k, ] <- c(p_try[i_t, ], OBJ(p_try[i_t, ]))
      pb_set(pb, k)
    }
  }
  pb_close(pb)
  sort_coll(out[is.finite(out[, N_FIT + 1]), , drop = FALSE])
}

## 7.4 Profile likelihoods (BYOM calc_proflik_ps.m) ---------------------------------
ps_profiles <- function(coll_all, S, label) {
  L <- N_FIT + 1
  gridsz <- 50
  mll <- coll_all[1, L]; mll_rem <- mll
  chicrit_single <- 0.5 * S$crit_table[1]
  chicrit_joint  <- max(0.5 * S$crit_table[N_FIT], chicrit_single + 1.5)
  ind_final <- sum(coll_all[, L] < mll + chicrit_joint)
  cloud <- coll_all[1:ind_final, 1:N_FIT, drop = FALSE]
  edges <- cbind(apply(cloud, 2, min), apply(cloud, 2, max))
  ok_list  <- list()
  par_best <- NULL
  coll_prof <- vector("list", N_FIT)
  pb <- pb_new(N_FIT * (2 * gridsz - 1), label); k <- 0
  better <- function(v, what) {      # report a better optimum found while profiling
    if (mll_rem - v > S$real_better && mll - v > S$real_better)
      say("    Better optimum found when ", what, ": ", fmt(v, 6), " (best was ", fmt(mll, 6), ")")
  }
  for (i_p in seq_len(N_FIT)) {
    par_rng <- seq(edges[i_p, 1], edges[i_p, 2], length.out = gridsz)
    gridsp  <- 0.5 * diff(edges[i_p, ]) / (gridsz - 1)
    ind_best <- max(which(par_rng <= coll_all[1, i_p]))
    par_rng  <- par_rng + (coll_all[1, i_p] - par_rng[ind_best])   # include best value
    parprof <- matrix(NA_real_, gridsz, L)
    for (i_g in seq_len(gridsz)) {              # forward pass
      ind_tst <- which(coll_all[, i_p] > par_rng[i_g] - gridsp &
                       coll_all[, i_p] < par_rng[i_g] + gridsp)[1]
      if (is.na(ind_tst)) {
        ind_tst <- which.min(abs(coll_all[, i_p] - par_rng[i_g])); mll_compare <- Inf
      } else mll_compare <- coll_all[ind_tst, L]
      p1 <- coll_all[ind_tst, 1:N_FIT]; p1[i_p] <- par_rng[i_g]
      r1 <- simplex_fixed(p1, i_p, rough = TRUE)
      if (i_g > 1) {
        p2 <- parprof[i_g - 1, 1:N_FIT]; p2[i_p] <- par_rng[i_g]
        r2 <- simplex_fixed(p2, i_p, rough = TRUE)
      } else r2 <- list(value = Inf)
      st <- if (r2$value < r1$value) r2$par else r1$par
      r  <- simplex_fixed(st, i_p, rough = FALSE)
      parprof[i_g, ] <- c(r$par, r$value)
      if (r$value < mll + chicrit_single + 1) {
        if (mll_compare - r$value > S$gap_extra) {
          ok_list[[length(ok_list) + 1]] <- coll_all[ind_tst, ]
          ok_list[[length(ok_list) + 1]] <- parprof[i_g, ]
        }
        if (r$value < mll) {
          par_best <- parprof[i_g, ]; better(r$value, paste("profiling", FIT_NAMES[i_p])); mll <- r$value
        }
      }
      k <- k + 1; pb_set(pb, k)
    }
    for (i_g in (gridsz - 1):1) {               # reverse pass
      mll_tst1 <- parprof[i_g, L]
      p2 <- parprof[i_g + 1, 1:N_FIT]; p2[i_p] <- parprof[i_g, i_p]
      r2 <- simplex_fixed(p2, i_p, rough = FALSE)
      if (r2$value < mll_tst1) {
        parprof[i_g, ] <- c(r2$par, r2$value)
        if (mll_tst1 - r2$value > S$gap_extra) ok_list[[length(ok_list) + 1]] <- parprof[i_g, ]
        if (r2$value < mll) {
          par_best <- parprof[i_g, ]; better(r2$value, paste("reverse profiling", FIT_NAMES[i_p])); mll <- r2$value
        }
      }
      k <- k + 1; pb_set(pb, k)
    }
    n_ext <- 0                                    # extend profile to lower values
    while (gridsp > 0 && n_ext < 500 && parprof[1, i_p] > BNDS[i_p, 1] &&
           parprof[1, L] < mll + chicrit_single + 1) {
      n_ext <- n_ext + 1
      p <- parprof[1, 1:N_FIT]; p[i_p] <- max(p[i_p] - gridsp, BNDS[i_p, 1])
      r <- simplex_fixed(p, i_p, rough = TRUE)
      r <- simplex_fixed(r$par, i_p, rough = FALSE)
      parprof <- rbind(c(r$par, r$value), parprof)
      ok_list[[length(ok_list) + 1]] <- parprof[1, ]
      if (r$value < mll) { par_best <- parprof[1, ]; better(r$value, paste("extending", FIT_NAMES[i_p], "down")); mll <- r$value }
    }
    if (n_ext > 0) say("    Profile for ", FIT_NAMES[i_p], " extended to lower values (", n_ext, " points)")
    n_ext <- 0                                    # extend profile to higher values
    while (gridsp > 0 && n_ext < 500 && parprof[nrow(parprof), i_p] < BNDS[i_p, 2] &&
           parprof[nrow(parprof), L] < mll + chicrit_single + 1) {
      n_ext <- n_ext + 1
      p <- parprof[nrow(parprof), 1:N_FIT]; p[i_p] <- min(p[i_p] + gridsp, BNDS[i_p, 2])
      r <- simplex_fixed(p, i_p, rough = TRUE)
      r <- simplex_fixed(r$par, i_p, rough = FALSE)
      parprof <- rbind(parprof, c(r$par, r$value))
      ok_list[[length(ok_list) + 1]] <- parprof[nrow(parprof), ]
      if (r$value < mll) { par_best <- parprof[nrow(parprof), ]; better(r$value, paste("extending", FIT_NAMES[i_p], "up")); mll <- r$value }
    }
    if (n_ext > 0) say("    Profile for ", FIT_NAMES[i_p], " extended to higher values (", n_ext, " points)")
    coll_prof[[i_p]] <- parprof
  }
  pb_close(pb)
  coll_ok <- if (length(ok_list)) do.call(rbind, ok_list) else matrix(0, 0, L)
  phat <- NULL
  if (!is.null(par_best)) {                       # better optimum: polish it
    say("    Finished profiling, running a simplex on the best set found ...")
    r <- run_simplex(par_best[1:N_FIT], rough = FALSE)
    phat <- r$par; mll <- r$value
  }
  flag_short <- FALSE                              # pruned profiles for the CIs
  prof_pruned <- vector("list", N_FIT)
  for (i_p in seq_len(N_FIT)) {
    tmp <- coll_prof[[i_p]][, c(i_p, L), drop = FALSE]
    tmp <- tmp[tmp[, 1] <= BNDS[i_p, 2] & tmp[, 1] >= BNDS[i_p, 1], , drop = FALSE]
    thr <- if (N_FIT > 1) mll + chicrit_joint else mll + chicrit_joint + 0.5
    ii  <- which(tmp[, 2] < thr)
    lo  <- if (length(ii)) max(1, min(ii) - 1) else 1
    hi  <- if (length(ii)) min(nrow(tmp), max(ii) + 1) else nrow(tmp)
    pr  <- tmp[lo:hi, , drop = FALSE]
    pr[is.infinite(pr[, 2]), 2] <- mll + 100
    prof_pruned[[i_p]] <- pr
    if (nrow(pr) < 40) flag_short <- TRUE
  }
  list(coll_prof = coll_prof, coll_ok = coll_ok, phat = phat, mll = mll,
       prof_pruned = prof_pruned, flag_short = flag_short)
}

## 7.5 Check sample against the profiles (BYOM test_proflik_ps.m) ------------------
ps_test_profiles <- function(coll_prof, coll_all, S) {
  L <- N_FIT + 1
  ok_list <- list(); flag <- c(0, 0)
  chicrit_single <- 0.5 * S$crit_table[1]
  mll <- coll_all[1, L]
  for (i_p in seq_len(N_FIT)) {
    parprof  <- coll_prof[[i_p]]
    coll_tst <- coll_all[coll_all[, L] < mll + chicrit_single + 1, , drop = FALSE]
    if (min(coll_tst[, i_p]) < min(parprof[, i_p]) || max(coll_tst[, i_p]) > max(parprof[, i_p])) {
      flag[1] <- flag[1] + 1; flag[2] <- Inf
    }
    parprof <- parprof[parprof[, L] < mll + chicrit_single + 1, , drop = FALSE]
    if (nrow(parprof) >= 3) for (i_g in 2:(nrow(parprof) - 1)) {
      gsp <- 0.5 * diff(parprof[c(i_g - 1, i_g, i_g + 1), i_p])
      ind_tst <- which(coll_all[, i_p] > parprof[i_g, i_p] - gsp[1] &
                       coll_all[, i_p] < parprof[i_g, i_p] + gsp[2])[1]
      if (is.na(ind_tst)) {
        ind_tst <- which.min(abs(coll_all[, i_p] - parprof[i_g, i_p])); mll_compare <- Inf
      } else mll_compare <- coll_all[ind_tst, L]
      mll_prof <- parprof[i_g, L]
      if (mll_prof < mll_compare) {
        if (mll_compare - mll_prof > S$gap_extra) {
          ok_list[[length(ok_list) + 1]] <- coll_all[ind_tst, ]
          ok_list[[length(ok_list) + 1]] <- parprof[i_g, ]
        }
      } else if (mll_compare < mll_prof) {
        min_MLL <- min(parprof[c(i_g - 1, i_g, i_g + 1), L])
        if (mll_compare < min_MLL) { flag[1] <- flag[1] + 1; flag[2] <- max(flag[2], min_MLL - mll_compare) }
      }
    }
  }
  list(flag = flag, coll_ok = if (length(ok_list)) do.call(rbind, ok_list) else matrix(0, 0, L))
}

## 7.6 The parameter-space explorer (BYOM calc_parspace.m) ----------------------------
# Returns the best fit, the 95% CIs and the sample of parameter sets (fit scale).
run_parspace <- function(S, do_profiles) {
  L <- N_FIT + 1
  n_ok   <- S$n_ok[N_FIT]
  n_conf <- S$n_conf_all[N_FIT, ]
  n_tr <- S$n_tr; f_d <- S$f_d
  chicrit_joint  <- 0.5 * S$crit_table[N_FIT]   # joint 95% CI of all parameters
  chicrit_single <- 0.5 * S$crit_table[1]       # single-parameter 95% CI ("inner rim")
  chicrit_rnd    <- chicrit_joint + S$crit_add
  chicrit_max    <- max(1.2 * chicrit_joint, chicrit_single + 1.5)
  status <- function(i_f, i_s, best) say(sprintf("    Status: %d sets within total CI and %d within inner. Best fit: %s",
                                                 i_f, i_s, fmt(best, 7)))

  ## --- Round 1: regular grid (or latin hypercube for > 4 parameters) ---
  n_rnd <- 1
  d_grid <- (BNDS[, 2] - BNDS[, 1]) / (S$tries - 1)
  if (N_FIT < 5 || (N_FIT == 5 && PS_ROUGHNESS == 0)) {
    grid <- as.matrix(expand.grid(lapply(seq_len(N_FIT), function(i)
      seq(BNDS[i, 1], BNDS[i, 2], length.out = S$tries))))
  } else {
    n_lhs <- if (PS_ROUGHNESS == 0) (if (N_FIT < 6) 30000 else if (N_FIT < 8) 60000 else 100000) else 10000
    n_ok  <- min(n_ok, 500)
    grid <- sapply(seq_len(N_FIT), function(i)
      BNDS[i, 1] + (sample(n_lhs) - runif(n_lhs)) / n_lhs * (BNDS[i, 2] - BNDS[i, 1]))
    say("  Using a latin-hypercube sample in the first round.")
  }
  colnames(grid) <- NULL
  say(sprintf("  Starting round 1 with initial grid of %d parameter sets", nrow(grid)))
  coll_all <- cbind(grid, Inf)
  pb <- pb_new(nrow(grid), "Round 1: initial grid")
  for (i in seq_len(nrow(grid))) { coll_all[i, L] <- OBJ(grid[i, ]); pb_set(pb, i) }
  pb_close(pb)
  coll_all <- sort_coll(coll_all[is.finite(coll_all[, L]), , drop = FALSE])
  if (!nrow(coll_all)) stop("All parameter sets in the initial grid gave an infinite MLL. Check data/ranges.")
  mll <- coll_all[1, L]
  ind_cont <- sum(coll_all[, L] < mll + chicrit_rnd[1])
  ind_cont <- min(max(ind_cont, n_ok), nrow(coll_all))
  if (ind_cont > n_conf[1]) {
    ind_cont2 <- which(coll_all[, L] - mll > chicrit_rnd[2])[1]
    if (is.na(ind_cont2)) ind_cont2 <- nrow(coll_all)
    ind_cont <- max(n_conf[1], ind_cont2)
  }
  coll_ok <- coll_all[1:ind_cont, , drop = FALSE]
  say("    Status: best fit so far is (minloglik) ", fmt(mll, 7))

  n_rnd <- 2
  n_tr_i <- n_tr[n_rnd]; f_d_i <- f_d[n_rnd]; chicrit_i <- chicrit_rnd[n_rnd]
  if (ind_cont > 0.5 * n_conf[1]) {
    n_tr_i <- floor(n_tr_i / 2)
    if (ind_cont >= n_conf[1]) n_tr_i <- floor(n_tr_i / 2)
  }
  n_tr_i <- min(n_tr_i, max(2, floor(10 * n_conf[1] / ind_cont)))
  coll_all <- coll_all[coll_all[, L] < mll + chicrit_max, , drop = FALSE]
  flag_stop <- FALSE; flag_inner <- 0

  ## --- Rounds 2 ... : random mutations around the best sets ---
  while (!flag_stop) {
    say(sprintf("  Starting round %d, refining a selection of %d parameter sets, with %d tries each",
                n_rnd, nrow(coll_ok), n_tr_i))
    coll_tries <- ps_mutations(coll_ok, n_tr_i, f_d_i * d_grid, sprintf("Round %d: mutations", n_rnd))
    coll_all <- sort_coll(rbind(coll_all, coll_tries))
    r <- run_simplex(coll_all[1, 1:N_FIT], rough = TRUE)
    mll <- r$value
    coll_all <- rbind(c(r$par, mll), coll_all)
    ind_final  <- sum(coll_all[, L] < mll + chicrit_joint)
    ind_single <- sum(coll_all[, L] < mll + chicrit_single)
    ind_cont_t <- sum(coll_tries[, L] < mll + chicrit_i)
    ind_cont_a <- sum(coll_all[, L] < mll + chicrit_i)
    ind_inner  <- sum(coll_all[, L] < mll + chicrit_single + 0.2)
    ind_cont_t2 <- sum(coll_tries[, L] < mll + chicrit_max)
    ind_cont_a2 <- sum(coll_all[, L] < mll + chicrit_max)
    status(ind_final, ind_single, coll_all[1, L])

    if (ind_final >= n_conf[1]) {                      # enough in the joint CI
      if (ind_single >= n_conf[2]) {
        flag_stop <- TRUE                              # ... and in the inner rim: done
      } else {
        say("    Next round will focus on inner rim (outer rim has enough points)")
        coll_ok <- coll_all[1:ind_inner, , drop = FALSE]
        flag_inner <- 1
        if (ind_inner < n_ok) {
          coll_ok <- rows_ab(coll_all, 1, n_ok)
        } else if (ind_inner > 0.5 * n_conf[2]) {
          coll_ok <- rbind(rows_ab(coll_all, ind_single - n_ok, ind_single + n_ok), coll_all[1, ])
        }
      }
    } else {                                           # not enough in the joint CI yet
      if (ind_cont_t > n_ok) {
        if (ind_cont_t > 2 * n_conf[1]) {
          if (ind_final > 0.5 * n_conf[1]) {
            coll_ok <- rows_ab(coll_tries, ind_single - n_ok, min(ind_cont_t, ind_single + n_ok))
          } else {
            ind_cont_t <- if (coll_tries[2 * n_conf[1], L] > mll + chicrit_max) 2 * n_conf[1] else ind_cont_t2
            coll_ok <- coll_tries[1:ind_cont_t, , drop = FALSE]
          }
        } else coll_ok <- coll_tries[1:ind_cont_t, , drop = FALSE]
        coll_ok <- rbind(coll_ok, coll_all[1, ])
      } else {
        if (ind_cont_a > n_ok) {
          if (ind_cont_a > 2 * n_conf[1])
            ind_cont_a <- if (coll_all[2 * n_conf[1], L] > mll + chicrit_max) 2 * n_conf[1] else ind_cont_a2
          coll_ok <- coll_all[1:ind_cont_a, , drop = FALSE]
        } else coll_ok <- rows_ab(coll_all, 1, n_ok)
      }
    }
    if (n_rnd == S$n_max && !flag_stop) {
      flag_stop <- TRUE
      say("  We have now done ", n_rnd, " rounds without reaching the stopping criterion ... we stop here!")
    }
    if (!flag_stop) {                                  # settings for the next round
      n_rnd <- n_rnd + 1
      j <- min(n_rnd, length(n_tr))
      n_tr_i <- n_tr[j]; f_d_i <- f_d[j]; chicrit_i <- chicrit_rnd[j]
      crit_ntry <- c(ind_final / n_conf[1], ind_single / n_conf[2])
      if (crit_ntry[1 + flag_inner] > 0.75 || nrow(coll_ok) > 2000) {
        n_tr_i <- floor(n_tr_i / 2)
      } else if (n_rnd > 3) {
        if (nrow(coll_ok) < 0.5 * n_conf[1 + flag_inner] && nrow(coll_ok) < 1000) {
          n_tr_i <- 2 * n_tr_i
          if (n_rnd > 4 && crit_ntry[1 + flag_inner] < 0.25) n_tr_i <- 2 * n_tr_i
        } else if (n_rnd > 7) n_tr_i <- 2 * n_tr_i
      }
      n_tr_i <- min(n_tr_i, max(2, floor(10 * n_conf[1 + flag_inner] / nrow(coll_ok))))
    }
    coll_all <- coll_all[coll_all[, L] < mll + chicrit_max, , drop = FALSE]
  }

  ## --- Final simplex on the best set ---
  say("  Finished sampling, running a simplex optimisation ...")
  r <- run_simplex(coll_all[1, 1:N_FIT], rough = FALSE)
  mll <- r$value
  coll_all <- rbind(c(r$par, mll), coll_all)
  ind_final  <- sum(coll_all[, L] < mll + chicrit_joint)
  ind_single <- sum(coll_all[, L] < mll + chicrit_single)
  status(ind_final, ind_single, mll)

  ci <- matrix(NA_real_, N_FIT, 4, dimnames = list(FIT_NAMES, c("lo", "hi", "joint_lo", "joint_hi")))
  flags <- matrix(0, N_FIT, 3, dimnames = list(FIT_NAMES, c("lo_bound", "hi_bound", "broken")))
  coll_prof <- NULL; prof_pruned <- NULL

  if (!do_profiles) {
    ## --- no profiling: CIs from the edges of the sample (approximate) ---
    i1 <- max(1, sum(coll_all[, L] < mll + 0.5 * S$crit_prop[1]))
    i2 <- sum(coll_all[, L] < mll + 0.5 * S$crit_prop[2])
    if (i2 > i1) {
      ci[, 1] <- apply(coll_all[i1:i2, 1:N_FIT, drop = FALSE], 2, min)
      ci[, 2] <- apply(coll_all[i1:i2, 1:N_FIT, drop = FALSE], 2, max)
    }
    ci[, 3] <- apply(coll_all[1:ind_final, 1:N_FIT, drop = FALSE], 2, min)
    ci[, 4] <- apply(coll_all[1:ind_final, 1:N_FIT, drop = FALSE], 2, max)
    say("  Profiling skipped: CIs are the edges of the inner cloud (may underestimate the true CIs).")
  } else {
    ## --- profile likelihoods + extra sampling where profile and sample disagree ---
    n_rnd <- n_rnd + 1
    say(sprintf("  Starting round %d, creating the profile likelihoods for each parameter", n_rnd))
    PR <- ps_profiles(coll_all, S, sprintf("Round %d: profile likelihoods", n_rnd))
    coll_prof <- PR$coll_prof; coll_ok <- PR$coll_ok; prof_pruned <- PR$prof_pruned
    flag_short <- PR$flag_short
    if (!is.null(PR$phat)) {
      mll <- PR$mll
      coll_all <- rbind(c(PR$phat, mll), coll_all)
      if (coll_all[2, L] - coll_all[1, L] > S$real_better) coll_ok <- rbind(coll_all[1, ], coll_ok)
      ind_final  <- sum(coll_all[, L] < mll + chicrit_joint)
      ind_single <- sum(coll_all[, L] < mll + chicrit_single)
    }
    status(ind_final, ind_single, coll_all[1, L])
    if (nrow(coll_ok) > 0 || ind_single < n_conf[2]) {
      if (nrow(coll_ok) == 0 || (ind_single < n_conf[2] && nrow(coll_ok) < 10))
        coll_ok <- rbind(coll_ok, rows_ab(coll_all, ind_single - n_ok, ind_single + n_ok))
      tmp <- rbind(coll_ok[, 1:N_FIT, drop = FALSE], coll_all[1:ind_final, 1:N_FIT, drop = FALSE])
      # NB: BYOM 6.2 used min() for both edges here (giving zero jump size);
      # the range (max - min) of the cloud is clearly what is intended.
      d_grid <- S$d_extra * (apply(tmp, 2, max) - apply(tmp, 2, min))
      f_d_i  <- S$f_d_extra
      coll_ok <- sort_coll(coll_ok)
      coll_ok <- coll_ok[coll_ok[, L] < mll + S$crit_add_extra, , drop = FALSE]
      n_rnd_x <- 0; n_test <- 0
      if (ind_single < n_conf[2]) {
        say("    Profiling has led to a new optimum, which left insufficient sets in the inner rim: extra sampling rounds will be started.")
      } else if (nrow(coll_ok) > 0) {
        say("    Profiling has detected gaps between profile and sample, which requires extra sampling rounds.")
      }
      while (nrow(coll_ok) > 0) {
        if (n_rnd >= 60) {                             # safety net (not in BYOM)
          say("  Stopping extra sampling after ", n_rnd, " rounds. Check the parameter-space plot.")
          break
        }
        n_cont <- nrow(coll_ok)
        n_rnd <- n_rnd + 1; n_rnd_x <- n_rnd_x + 1
        n_tr_i <- min(40, max(2, floor(5 * n_conf[1] / n_cont)))
        say(sprintf("  Starting round %d, (extra %d) refining a selection of %d parameter sets, with %d tries each",
                    n_rnd, n_rnd_x, n_cont, n_tr_i))
        coll_tries <- ps_mutations(coll_ok, n_tr_i, f_d_i * d_grid, sprintf("Round %d: extra sampling", n_rnd))
        coll_all <- sort_coll(rbind(coll_all, coll_tries))
        coll_all <- coll_all[coll_all[, L] < mll + chicrit_max, , drop = FALSE]
        ind_final  <- sum(coll_all[, L] < mll + chicrit_joint)
        ind_single <- sum(coll_all[, L] < mll + chicrit_single)
        status(ind_final, ind_single, coll_all[1, L])
        TP <- ps_test_profiles(coll_prof, coll_all, S)
        coll_ok <- TP$coll_ok
        if (ind_single < n_conf[2]) {
          ind_cont <- sum(coll_tries[, L] < mll + S$crit_add_extra)
          if (ind_cont < 0.5 * n_conf[2]) coll_ok <- rbind(coll_ok, rows_ab(coll_tries, 1, ind_cont))
          else coll_ok <- rows_ab(coll_tries, 1, ind_cont)
          f_d_i <- max(0.1, f_d_i * 0.7)
        } else if (nrow(coll_ok) > 0 && n_rnd_x < 5) {
          f_d_i <- max(0.1, f_d_i * 0.7)
        } else if ((TP$flag[2] > 0.05 || flag_short) && n_test < 2) {
          n_rnd <- n_rnd + 1; n_test <- n_test + 1
          say(sprintf("  Starting round %d, refining the profile likelihoods for each parameter again", n_rnd))
          PR <- ps_profiles(coll_all, S, sprintf("Round %d: profile likelihoods", n_rnd))
          coll_prof <- PR$coll_prof; coll_ok <- PR$coll_ok; prof_pruned <- PR$prof_pruned
          flag_short <- PR$flag_short
          if (!is.null(PR$phat)) {
            mll <- PR$mll
            coll_all <- rbind(c(PR$phat, mll), coll_all)
            ind_final  <- sum(coll_all[, L] < mll + chicrit_joint)
            ind_single <- sum(coll_all[, L] < mll + chicrit_single)
          }
          status(ind_final, ind_single, coll_all[1, L])
          if (nrow(coll_ok) > 0) {
            n_rnd_x <- 0; f_d_i <- S$f_d_extra
            say("    New profiling has detected gaps, which requires additional extra sampling rounds!")
          }
        } else {
          if (nrow(coll_ok) > 0)
            say("  Exiting parameter space explorer, but still ", nrow(coll_ok),
                " sets flagged (too much distance between profile and sample). Check plot.")
          if (TP$flag[2] > 0.05)
            say("  Exiting parameter space explorer, but there still seem to be parameter sets outside of the profile curves. Check plot.")
          coll_ok <- coll_ok[0, , drop = FALSE]
        }
      }
    }
    ## --- CIs from the profile likelihoods ---
    coll_all <- coll_all[coll_all[, L] < mll + chicrit_max, , drop = FALSE]
    ci[, 3] <- apply(coll_all[1:ind_final, 1:N_FIT, drop = FALSE], 2, min)
    ci[, 4] <- apply(coll_all[1:ind_final, 1:N_FIT, drop = FALSE], 2, max)
    for (i_p in seq_len(N_FIT)) {
      coll_all <- sort_coll(rbind(coll_all, coll_prof[[i_p]]))   # profile points join the sample
      ci_res <- ci_from_profile(prof_pruned[[i_p]], mll, chicrit_single, BNDS[i_p, ])
      ci[i_p, 1:2] <- ci_res$ci
      flags[i_p, ] <- ci_res$flags
    }
  }
  best <- coll_all[1, 1:N_FIT]
  coll_all <- coll_all[coll_all[, L] < mll + max(chicrit_joint, 0.5 * S$crit_prop[2]), , drop = FALSE]
  list(best = best, mll = coll_all[1, L], sample = coll_all, ci = ci, flags = flags,
       profiles = coll_prof, n_joint = sum(coll_all[, L] < coll_all[1, L] + chicrit_joint),
       n_inner = sum(coll_all[, L] < coll_all[1, L] + chicrit_single))
}

## 7.7 95% CI from a profile: where it crosses MLL + 1.92 -----------------------------
ci_from_profile <- function(prof, mll, chicrit, bnds) {
  pr <- prof[order(prof[, 1]), , drop = FALSE]
  y  <- pr[, 2] - mll - chicrit
  ok <- is.finite(y); pr <- pr[ok, , drop = FALSE]; y <- y[ok]
  fl <- c(lo_bound = 0, hi_bound = 0, broken = 0)
  ci <- c(NA_real_, NA_real_)
  if (!length(y) || !any(y < 0)) return(list(ci = ci, flags = fl))
  cross <- function(i, j) pr[i, 1] + (0 - y[i]) * (pr[j, 1] - pr[i, 1]) / (y[j] - y[i])
  if (y[1] < 0) { ci[1] <- bnds[1]; fl[1] <- 1 } else {
    i <- which(y < 0)[1]; ci[1] <- cross(i - 1, i)
  }
  n <- length(y)
  if (y[n] < 0) { ci[2] <- bnds[2]; fl[2] <- 1 } else {
    i <- max(which(y < 0)); ci[2] <- cross(i, i + 1)
  }
  if (sum(diff(sign(y)) != 0) + fl[1] + fl[2] > 2) fl[3] <- 1     # broken CI
  list(ci = ci, flags = fl)
}

## 7.8 Profile likelihoods after a simplex fit (FIT_METHOD = "simplex") -------------
# Simple stepwise profiling from the best fit to both sides until the MLL
# exceeds the 95% criterion (+1) or a bound is hit; step size adapts.
simplex_profiles <- function(best, mll, label) {
  chicrit <- 0.5 * 3.8415
  profs <- vector("list", N_FIT)
  pb <- pb_new(N_FIT * 2, label)
  k <- 0; better_found <- NULL
  for (i_p in seq_len(N_FIT)) {
    pts <- rbind(c(best, mll))
    for (dirn in c(-1, 1)) {
      p <- best; v_last <- mll
      step <- 0.02 * (BNDS[i_p, 2] - BNDS[i_p, 1])
      n <- 0
      repeat {
        n <- n + 1
        x_new <- p[i_p] + dirn * step
        x_new <- min(max(x_new, BNDS[i_p, 1]), BNDS[i_p, 2])
        q <- p; q[i_p] <- x_new
        r <- simplex_fixed(q, i_p, rough = TRUE)
        r <- simplex_fixed(r$par, i_p, rough = FALSE)
        dv <- r$value - v_last
        if (dv > 0.5 && step > 1e-4 * (BNDS[i_p, 2] - BNDS[i_p, 1])) { step <- step / 2; next }
        pts <- rbind(pts, c(r$par, r$value))
        p <- r$par; v_last <- r$value
        if (r$value < mll - 0.01 && (is.null(better_found) || r$value < better_found[N_FIT + 1]))
          better_found <- c(r$par, r$value)
        if (r$value > mll + chicrit + 1 || x_new <= BNDS[i_p, 1] || x_new >= BNDS[i_p, 2] || n > 200) break
        if (dv < 0.05) step <- step * 2
      }
      k <- k + 1; pb_set(pb, k)
    }
    profs[[i_p]] <- pts[order(pts[, i_p]), , drop = FALSE]
  }
  pb_close(pb)
  list(profiles = profs, better = better_found)
}


###############################################################################
## SECTION 8 : Functions for results, metrics, plots (used in SECTION 9)
###############################################################################

CHI_PROP <- 0.5 * (3.8415 + 0.3764)   # BYOM: sets with 2*dMLL < chi2(1)+0.3764 -> CI bands
Z95      <- qnorm(0.975)              # 1.96 for the 95% prediction interval
NAMES_ST <- c(parent = PARENT_NAME, btp = BTP_NAME, total = "Total (parent + BTP)")
UNITS <- c(ku = "L/kg ww/d", ke = "1/d", km = "1/d", kem = "1/d", kum = "L/kg ww/d")
back <- function(v, nm) if (PAR[nm, "scale"] == "log") 10^v else v

## 8.1 Derived quantities (incl. BCF) ---------------------------------------------
derived <- function(p) {
  k1 <- p[["ke"]] + p[["km"]]
  d <- c(BCF_parent      = p[["ku"]] / k1,
         t_half_parent_d = log(2) / k1,
         t95_parent_d    = log(20) / k1)
  if (USE_BIOTRANSFORMATION)
    d <- c(d, BCF_total_molar         = p[["ku"]] / k1 * (1 + p[["km"]] / p[["kem"]]),
              fraction_biotransformed = p[["km"]] / k1,
              t_half_BTP_d            = log(2) / p[["kem"]],
              BTP_to_parent_ratio_ss  = p[["km"]] / p[["kem"]])
  d
}
## 8.2 Goodness-of-fit metrics ---------------------------------------------------------
# PI_LO/PI_HI: 95% prediction interval at the observations (transformed scale),
# DET: best-fit details, T_EOE: end of exposure per scenario
calc_metrics <- function() {
  out <- list()
  for (st in FIT_STATES) {
    o    <- OBS[[st]]
    pred <- pmax(DET$pred[o$row, st], 0)
    tD   <- transform_data(o$value)
    inside <- tD >= PI_LO[o$row, st] - 1e-12 & tD <= PI_HI[o$row, st] + 1e-12
    for (sc in c("all", SCEN_IDS)) {
      s <- if (sc == "all") rep(TRUE, nrow(o)) else o$scenario == sc
      if (!any(s)) next
      obs <- o$value[s]; pr <- pred[s]
      r2 <- if (sc == "all") DET$r2[[st]] else {
        Ds <- tD[s]; 1 - sum((Ds - transform_data(pr))^2) / sum((Ds - mean(Ds))^2) }
      nrmse <- if (mean(obs) > 0) 100 * sqrt(mean((pr - obs)^2)) / mean(obs) else NA
      ppc   <- 100 * mean(inside[s])
      # PE at the end of exposure (per scenario; "all" = mean over scenarios,
      # acceptable only if every scenario has |PE| <= PE_MAX)
      pe_sc <- sapply(if (sc == "all") SCEN_IDS else sc, function(s2) {
        os <- o[o$scenario == s2, ]; ps <- pred[o$scenario == s2]
        tt <- os$time[os$time <= T_EOE[[s2]] + 1e-9]
        if (!length(tt)) return(c(NA, NA))
        t_s <- max(tt); sel <- abs(os$time - t_s) < 1e-9
        m_obs <- mean(os$value[sel])
        c(if (m_obs > 0) 100 * (mean(ps[sel]) - m_obs) / m_obs else NA, t_s)
      })
      pe_sc <- matrix(pe_sc, nrow = 2)
      out[[length(out) + 1]] <- data.frame(
        state = NAMES_ST[[st]], scenario = if (sc == "all") "all (pooled)" else sc,
        n = sum(s), R2 = r2, NRMSE_pct = nrmse, PPC_pct = ppc,
        PE_pct = if (all(is.na(pe_sc[1, ]))) NA else mean(pe_sc[1, ], na.rm = TRUE),
        t_PE_d = if (sc == "all") NA else pe_sc[2, 1],
        PPC_ok = ppc >= PPC_MIN, NRMSE_ok = !is.na(nrmse) && nrmse <= NRMSE_MAX,
        PE_ok = all(!is.na(pe_sc[1, ])) && all(abs(pe_sc[1, ]) <= PE_MAX),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

## 8.3 Plot helpers ---------------------------------------------------------------------
show_and_save <- function(draw, name, width, height, n_pages = 1) {
  if (SHOW_PLOTS && interactive())           # a small plot window must not stop the script
    tryCatch(draw(), error = function(e)
      say("  NOTE: the plot could not be drawn in the plot window (", conditionMessage(e),
          "). Enlarge the plot window; the plot is saved as png."))
  if (SAVE_PLOTS) {
    if (!dir.exists(MODEL_DIR)) dir.create(MODEL_DIR, recursive = TRUE)
    grDevices::png(out_file(paste0(name, if (n_pages > 1) "_%02d", ".png")), width = width, height = height,
                   units = "in", res = 150); draw(); grDevices::dev.off()
  }
}
par_text <- function() {
  txt <- paste(sapply(FIT_NAMES, function(nm) sprintf("%s = %s", nm, fmt(P_BEST[[nm]]))), collapse = ",  ")
  sprintf("%s  |  MLL = %s, AIC = %s  |  %s", txt, fmt(MLL, 6), fmt(AIC, 5),
          switch(FIT$method, parspace = "parameter-space explorer", simplex = "simplex", none = "start values (not fitted)"))
}
PER_PAGE <- 4                        # scenarios (rows) per page in the fit plot
pretty_unit <- function(u) sub("^u", "\u00b5", u)                         # "ug/L" -> micro-g/L
shade <- function(colour) grDevices::adjustcolor(colour, alpha.f = 0.22)   # transparent colour (CI band)
# Empty panel with horizontal grid lines, a panel letter and a title
new_panel <- function(xlim, ylim, xlab, ylab, letter, heading) {
  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, yaxs = "i", xlab = "", ylab = "", bty = "l")
  grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
  box(bty = "l")
  title(xlab = xlab, line = 2.4)
  title(ylab = ylab, line = 3.6)
  mtext(letter, side = 3, line = 2, at = grconvertX(0.02, from = "nfc", to = "user"), adj = 0, font = 2, cex = 1.1)
  mtext(heading, side = 3, line = 0.5, adj = 0, font = 2, cex = 0.9)
}

## 8.4 Fit plot: exposure (left) and tissue (right), one row per scenario ---------------
# All concentrations are shown in the units of the data file (e.g. ug/L and
# ug/kg ww), also when the model is fitted in molar units.
#   line          : model with the best-fit parameters
#   shaded band   : 95% confidence interval of the model curve
#   points        : measured concentrations (all replicates)
#   dotted line   : end of exposure
draw_fits <- function() {
  pages <- split(SCEN_IDS, ceiling(seq_along(SCEN_IDS) / PER_PAGE))
  n_col <- 2
  op <- par(no.readonly = TRUE); on.exit(par(op))
  for (pg in pages) {
    layout(matrix(seq_len(n_col * length(pg)), nrow = length(pg), byrow = TRUE),
           widths = c(rep(1, n_col - 1), 1.7))
    par(mar = c(4, 5.2, 3.6, 1), oma = c(2.4, 0, 1.8, 0), mgp = c(2.4, 0.7, 0), las = 1)
    for (sc in pg) {
      C  <- CURVES[[sc]]
      ex <- EXPO[[sc]]
      k  <- (match(sc, pg) - 1) * n_col              # panel letters a, b, c, ...
      xl <- range(C$t)

      ## --- exposure medium: measured concentrations and the interpolation used by the model ---
      e  <- EXP[EXP$scenario == sc & EXP$time >= 0, ]
      cw <- forcing_at(ex$parent, C$t) / F_PAR
      new_panel(xl, c(0, 1.08 * max(c(cw, e$parent, 1e-12), na.rm = TRUE)), "Time [d]",
                paste0("Medium concentration [", pretty_unit(DATA_UNIT_MEDIUM), "]"), letters[k + 1],
                if (length(SCEN_IDS) > 1) paste0("Exposure medium (", sc, ")") else "Exposure medium")
      lines(C$t, cw, lwd = 2.5, col = CLR_PARENT)
      points(e$time, e$parent, pch = 21, bg = CLR_PARENT, col = "white", cex = 1.5, xpd = NA)
      if (BTP_UPTAKE_FROM_MEDIUM && !is.null(ex$btp)) {     # BTP in the medium (only shown if it is taken up)
        lines(C$t, forcing_at(ex$btp, C$t) / F_BTP, lwd = 2, lty = 2, col = CLR_SECOND)
        points(e$time, e$btp, pch = 21, bg = CLR_SECOND, col = "white", cex = 1.2, xpd = NA)
      }
      ## --- tissue: model curves (with 95% CI band) and measured concentrations ---
      sp    <- CATALOG$species[match(sc, CATALOG$scenario)]
      right <- C$t > 0.5 * max(C$t)                  # right half of the panel: room for the legend?
      ymax  <- 1e-12
      ymax_right <- 0
      for (s in CURVE_SPEC) {
        y <- C$best[, s$col] * s$f
        if (s$band && !is.null(C$hi)) y <- pmax(y, C$hi[, s$col] * s$f)
        ymax       <- max(ymax, y, na.rm = TRUE)
        ymax_right <- max(ymax_right, y[right], na.rm = TRUE)
        if (!is.na(s$obs)) {
          o <- OBS[[s$obs]][OBS[[s$obs]]$scenario == sc, ]
          ymax       <- max(ymax, o$value * s$f)
          ymax_right <- max(ymax_right, o$value[o$time > 0.5 * max(C$t)] * s$f)
        }
      }
      top <- if (ymax_right > 0.68 * ymax) ymax * (1.2 + 0.12 * length(CURVE_SPEC)) else 1.08 * ymax
      new_panel(xl, c(0, top), "Time [d]",
                paste0("Tissue concentration [", pretty_unit(DATA_UNIT_TISSUE), "]"), letters[k + n_col],
                bquote(bold(.(paste0(ORGANISM_MATRIX, " ("))) * bolditalic(.(sp)) * bold(")")))
      abline(v = T_EOE[[sc]], lty = 3, col = "grey30")
      for (s in CURVE_SPEC) if (s$band && !is.null(C$lo))
        polygon(c(C$t, rev(C$t)), c(C$lo[, s$col], rev(C$hi[, s$col])) * s$f, col = shade(s$colour), border = NA)
      for (s in CURVE_SPEC)
        lines(C$t, C$best[, s$col] * s$f, lwd = if (s$band) 2.5 else 2, lty = s$lty, col = s$colour)
      for (s in CURVE_SPEC) if (!is.na(s$obs)) {
        o <- OBS[[s$obs]][OBS[[s$obs]]$scenario == sc, ]
        points(o$time, o$value * s$f, pch = 21, bg = s$colour, col = "white", cex = 1.5, xpd = NA)
      }
      legend("topright", inset = 0.02, legend = sapply(CURVE_SPEC, function(s) s$label),
             col = sapply(CURVE_SPEC, function(s) s$colour), lty = sapply(CURVE_SPEC, function(s) s$lty),
             lwd = 2.5, pch = sapply(CURVE_SPEC, function(s) if (is.na(s$obs)) NA else 16),
             pt.cex = 1.3, bg = "white", box.col = "grey40", cex = 0.9)
    }
    mtext(paste0(RUN_NAME, " - ", GRP, ": ", MODEL_LABEL), outer = TRUE, line = 0.5, font = 2)
    mtext(par_text(), side = 1, outer = TRUE, line = 1, cex = 0.7)
  }
}

## 8.5 Parameter space: sample + profile likelihoods (as the BYOM/openGUTS plot) --------
# How to read this plot:
#   Diagonal panels : profile likelihood of one parameter (red line). x = value of
#                     the parameter, y = log-likelihood ratio (how much worse than
#                     the best fit). The 95% CI of the parameter is the range where
#                     the red line is below the solid horizontal line (1.92).
#                     A narrow, V-shaped profile = well-defined parameter.
#   Other panels    : the tested parameter sets for two parameters at a time.
#                     A round cloud = the two parameters are independent; a tilted,
#                     stretched cloud = they are correlated.
#   green           : parameter sets within the 95% CI of the single parameters
#   cyan            : parameter sets within the joint 95% confidence region
#   yellow          : best fit
#   blue line       : search bound (min/max in SECTION 1.6). If the cloud or the CI
#                     touches it: widen the range, or the data cannot define this parameter.
draw_parspace <- function() {
  smp   <- FIT$sample
  d     <- smp[, L_COL] - MLL                    # log-likelihood ratio to the best fit
  crit  <- 0.5 * 3.8415                          # 95% criterion for a single parameter (chi2, df = 1)
  inner <- d < crit
  ord   <- order(-d)                             # worst first, so the best sets are drawn on top
  ymax  <- max(1.3 * crit, min(max(d), 3 * crit))
  best  <- FIT$best
  lim <- lapply(seq_len(N_FIT), function(j) {    # plotted range of each parameter (on fit scale)
    v <- c(smp[, j], best[j])
    if (!is.null(FIT$profiles)) {
      pr <- FIT$profiles[[j]]
      v  <- c(v, pr[pr[, L_COL] - MLL < ymax, j])
    }
    r <- range(v)
    if (diff(r) <= 0) r <- r + c(-0.5, 0.5)
    r + c(-0.04, 0.04) * diff(r)
  })
  # Axis with tick labels in real parameter values (also for parameters fitted
  # on log scale). axisTicks() gives well-spaced ticks (1-2-5 pattern) on log
  # axes; ticks outside the plotted range are dropped, and a small gap.axis
  # lets neighbouring labels sit closer before any of them is left out.
  par_axis <- function(side, j, labels) {
    r <- lim[[j]]
    if (IS_LOG[j]) {
      tk <- grDevices::axisTicks(r, log = TRUE, nint = 5)
      if (length(tk) < 3)                        # very narrow range: linear ticks
        tk <- pretty(10^r, n = 5)
      tk <- tk[tk > 0]
      at <- log10(tk)
    } else {
      tk <- pretty(r, n = 5)
      at <- tk
    }
    keep <- at >= r[1] & at <= r[2]
    tk <- tk[keep]; at <- at[keep]
    if (!length(tk)) return(invisible(NULL))
    lab <- if (labels) sapply(tk, function(v) format(v, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) else FALSE
    ax <- function(...) if (getRversion() >= "4.0.0") axis(..., gap.axis = 0.25) else axis(...)
    ax(side, at = at, labels = lab, tcl = -0.3, col = "grey30")
  }
  panel <- function(xl, yl) {                    # empty white panel on the grey background
    plot.new()
    plot.window(xlim = xl, ylim = yl, xaxs = "i", yaxs = "i")
    rect(xl[1], yl[1], xl[2], yl[2], col = "white", border = NA)
  }
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = c(N_FIT, N_FIT), mar = c(0.8, 0.8, 0.4, 0.4), oma = c(4.5, 5.5, 3.4, 1), mgp = c(2.2, 0.6, 0),
      las = 1, bg = "grey92")
  for (i in seq_len(N_FIT)) for (j in seq_len(N_FIT)) {
    if (i == j) {                                # diagonal: profile likelihood of parameter j
      panel(lim[[j]], c(0, ymax))
      points(smp[, j], d, pch = 16, cex = 0.4, col = ifelse(inner, CLR_PS_INNER, CLR_PS_JOINT))
      abline(h = crit)                                           # 95% criterion
      abline(h = 0.5 * (3.8415 + c(-2, 1) * 0.3764), lty = 3)    # band around the criterion (as BYOM)
      if (!is.null(FIT$profiles)) {
        pr <- FIT$profiles[[j]]
        pr <- pr[order(pr[, j]), , drop = FALSE]
        lines(pr[, j], pr[, L_COL] - MLL, col = "red", lwd = 2)
      }
      abline(v = BNDS[j, ], col = "blue", lwd = 2)               # search bounds (only visible if reached)
      points(best[j], 0, pch = 21, bg = "yellow", cex = 1.8, xpd = NA)
      box(col = "grey30")
      par_axis(1, j, labels = (i == N_FIT))
      axis(2, at = pretty(c(0, ymax), n = 4), labels = (j == 1), tcl = -0.3, col = "grey30")
      if (j == 1) mtext("loglik ratio", side = 2, line = 3.2, las = 0, cex = 0.85)
    } else if (i > j) {                          # below the diagonal: parameter j (x) against parameter i (y)
      panel(lim[[j]], lim[[i]])
      points(smp[ord, j], smp[ord, i], pch = 21, cex = 0.9, lwd = 0.6,
             col = ifelse(inner[ord], "black", CLR_PS_JOINT), bg = ifelse(inner[ord], CLR_PS_INNER, NA))
      abline(v = BNDS[j, ], h = BNDS[i, ], col = "blue", lwd = 2)
      points(best[j], best[i], pch = 21, bg = "yellow", cex = 1.8)
      box(col = "grey30")
      par_axis(1, j, labels = (i == N_FIT))
      par_axis(2, i, labels = (j == 1))
      if (j == 1) mtext(FIT_NAMES[i], side = 2, line = 3.2, las = 0, cex = 0.85)
    } else {
      plot.new()                                 # above the diagonal: empty
    }
    if (i == N_FIT) mtext(FIT_NAMES[j], side = 1, line = 2.4, cex = 0.85)
  }
  mtext(paste0(RUN_NAME, " - ", GRP, ": parameter space"), outer = TRUE, line = 1.9, font = 2, cex = 0.95)
  mtext("green: within 95% CI of single parameters | cyan: joint 95% region | red: profile likelihood | yellow: best fit",
        outer = TRUE, line = 0.5, cex = 0.6)
}


## 8.6 Compartment plot: one panel per model curve, all scenarios overlaid --------------
# Written only for models with more than one curve (e.g. structure + lipid).
# All scenarios of the run share the panels, so pooled scenarios (e.g.
# different exposure concentrations) can be compared per compartment:
#   line          : model with the best-fit parameters
#   shaded band   : 95% confidence interval of the model curve
#   points        : measured concentrations (only the whole organism is measured)
#   dotted line   : end of exposure
SC_COLOURS <- c("#A858C8", "#4EA6F0", "#6E9A1F", "#D98E2B", "#C0392B", "#16A085",
                "#7F8C8D", "#2C3E50")                     # one colour per scenario
draw_compartments <- function() {
  n  <- length(CURVE_SPEC)
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = c(1, n), mar = c(4, 5.2, 3.6, 1), oma = c(0, 0, 2.2, 0),
      mgp = c(2.4, 0.7, 0), las = 1)
  cols <- setNames(rep_len(if (length(SCEN_IDS) > 1) SC_COLOURS else CLR_PARENT,
                           length(SCEN_IDS)), SCEN_IDS)
  xl <- range(sapply(CURVES, function(C) range(C$t)))
  for (k in seq_len(n)) {
    s <- CURVE_SPEC[[k]]
    ymax <- 1e-12                                # common y range of this panel
    for (sc in SCEN_IDS) {
      C <- CURVES[[sc]]
      ymax <- max(ymax, C$best[, s$col] * s$f,
                  if (!is.null(C$hi)) C$hi[, s$col] * s$f, na.rm = TRUE)
      if (!is.na(s$obs))
        ymax <- max(ymax, OBS[[s$obs]]$value[OBS[[s$obs]]$scenario == sc] * s$f)
    }
    new_panel(xl, c(0, 1.08 * ymax), "Time [d]",
              paste0("Concentration [", pretty_unit(DATA_UNIT_TISSUE), "]"),
              LETTERS[k], s$label)
    abline(v = unique(unlist(T_EOE[SCEN_IDS])), lty = 3, col = "grey30")
    for (sc in SCEN_IDS) {
      C <- CURVES[[sc]]
      if (!is.null(C$lo))
        polygon(c(C$t, rev(C$t)), c(C$lo[, s$col], rev(C$hi[, s$col])) * s$f,
                col = shade(cols[[sc]]), border = NA)
      lines(C$t, C$best[, s$col] * s$f, lwd = 2.5, col = cols[[sc]])
      if (!is.na(s$obs)) {
        o <- OBS[[s$obs]][OBS[[s$obs]]$scenario == sc, ]
        points(o$time, o$value * s$f, pch = 21, bg = shade(cols[[sc]]), col = cols[[sc]],
               cex = 1.3, xpd = NA)
      }
    }
    if (k == 1 && length(SCEN_IDS) > 1)
      legend("topright", inset = 0.02, legend = SCEN_IDS, col = cols, lwd = 2.5,
             bg = "white", box.col = "grey40", cex = 0.9)
  }
  mtext(paste0(RUN_NAME, " - ", GRP, ": compartments  (", MODEL_LABEL, ")"),
        outer = TRUE, line = 0.6, font = 2)
}
###############################################################################
## SECTION 9 : Fit every model run (pooled scenarios share one parameter set)
###############################################################################

SUMMARY  <- list()                # one row per model run
ALL_TAB  <- list()                # parameters of every model (Excel sheet 1, SECTION 10)
FIT_DATA <- list()                # observations + model curves (Excel sheet 3, SECTION 10)
SKIPPED <- character(0)           # scenarios / runs that could not be fitted
SOLVER_CHECKED <- FALSE

HEADER_END <- length(.LOG$lines)  # general part of the log (written into every model log)

for (JOB in seq_len(nrow(JOBS))) {

  GRP <- JOBS$run[JOB]
  setup_version(JOBS$bio[JOB])
  VERSION <- if (USE_BIOTRANSFORMATION) "bio" else "nobio"
  JOB_ID  <- if (length(BIO_VERSIONS) > 1) paste0(GRP, " [", VERSION, "]") else GRP
  TRT_TAG <- gsub("[^A-Za-z0-9_-]+", "_", if (length(BIO_VERSIONS) > 1) paste(GRP, VERSION) else GRP)
  LOG_FROM <- length(.LOG$lines) + 1
  say("")
  rule("#")
  say(sprintf("MODEL %d of %d: %s   (%s)   biotransformation %s", JOB, nrow(JOBS), GRP,
              paste(RUNS[[GRP]], collapse = " + "), if (USE_BIOTRANSFORMATION) "ON" else "OFF"))
  rule("#")

  ## 9.1 Data of the scenario(s) of this run --------------------------------------------
  SET <- list()
  for (sc in RUNS[[GRP]]) {
    say("  Scenario ", sc)
    st_sc <- setup_scenario(sc)
    if (!is.null(st_sc$skip)) {
      say("    SKIPPED: ", st_sc$skip)
      SKIPPED[sc] <- st_sc$skip; next
    }
    SET[[sc]] <- st_sc
  }
  if (!length(SET)) { SKIPPED[JOB_ID] <- "no usable scenario"; next }
  SCEN_IDS <- names(SET)
  MODEL_DIR <- file.path(OUTPUT_DIR, paste0(TRT_TAG, "_", RUN_STAMP))   # one output folder per model
  EXPO <- lapply(SET, function(x) x$expo)
  ORG  <- do.call(rbind, lapply(SET, function(x) x$org))
  EXP  <- do.call(rbind, lapply(SET, function(x) x$exp))

  OBS <- lapply(c(parent = "parent", btp = "btp"), make_obs)
  too_few <- FIT_STATES[sapply(FIT_STATES, function(st) nrow(OBS[[st]]) < 3)]
  if (length(too_few)) {
    say("  SKIPPED: fewer than 3 tissue data for: ", paste(NAMES_ST[too_few], collapse = ", "))
    SKIPPED[JOB_ID] <- "too few tissue data"; next
  }
  if (length(unique(OBS$parent$time[OBS$parent$time > 0])) < 2) {
    say("  SKIPPED: tissue data at fewer than 2 time points after t = 0 - rate constants cannot be estimated.")
    SKIPPED[JOB_ID] <- "tissue data at < 2 time points"; next
  }
  for (st in FIT_STATES) if (all(OBS[[st]]$value == 0))
    say("  WARNING: all tissue data for ", NAMES_ST[st], " are 0.")
  if (DATA_TRANSFORM == 0 && any(unlist(lapply(OBS[FIT_STATES], function(o) o$value == 0))))
    say("  WARNING: log-transformation with zeros in the data (zeros are set to 1e-10, as BYOM).")

  ## 9.2 Grids at the data time points, start values, end of exposure ---------------------
  GRID_DATA <- list(); OFFSET <- c(); T_EOE <- list(); n_rows <- 0
  for (sc in SCEN_IDS) {
    o_p <- OBS$parent[OBS$parent$scenario == sc, ]
    o_m <- OBS$btp[OBS$btp$scenario == sc, ]
    x0 <- if (X0_SOURCE == "data") c(parent = mean_t0(o_p), btp = mean_t0(o_m)) else
      c(parent = X0_PARENT * F_PAR, btp = X0_BTP * F_BTP)
    if (!USE_BIOTRANSFORMATION) x0[["btp"]] <- 0
    t_sc <- sort(unique(c(o_p$time, o_m$time)))
    GRID_DATA[[sc]] <- make_grid(EXPO[[sc]], t_sc, x0)
    OFFSET[sc] <- n_rows
    n_rows <- n_rows + length(t_sc)
    T_EOE[[sc]] <- end_of_exposure(sc, EXPO[[sc]]$parent, t_sc)
    say(sprintf("    medium %d time points | tissue %d %s, %d %s values | start: %s / %s %s | end of exposure: %s d",
                length(unique(round(EXPO[[sc]]$parent$time, 4))), nrow(o_p), PARENT_NAME, nrow(o_m), BTP_NAME,
                fmt(x0[["parent"]]), fmt(x0[["btp"]]), UNIT, fmt(T_EOE[[sc]], 5)))
    if (X0_SOURCE == "data" && !any(o_p$time == 0))
      say("    NOTE: no tissue data at t = 0 for ", sc, "; the model starts at 0.")
  }
  for (st in names(OBS)) {            # row of each observation in the stacked prediction
    o <- OBS[[st]]
    OBS[[st]]$row <- if (nrow(o)) unname(OFFSET[o$scenario]) +
      unname(mapply(function(sc, t) match(t, GRID_DATA[[sc]]$t_out), o$scenario, o$time)) else integer(0)
  }

  ## 9.3 Solver self-check (analytical vs. deSolve), once, if deSolve is available ---------
  if (!SOLVER_CHECKED && requireNamespace("deSolve", quietly = TRUE)) {
    p_chk <- P_FIXED
    if (USE_BIOTRANSFORMATION && p_chk[["km"]] == 0) p_chk[["km"]] <- 0.1
    d_a <- solve_analytical(p_chk, GRID_DATA[[1]])
    d_o <- solve_ode(p_chk, GRID_DATA[[1]])
    rel <- max(abs(d_a - d_o)) / max(1e-12, max(abs(d_o)))
    say("  Solver check (analytical vs. deSolve at start values): max. relative difference = ",
        formatC(rel, format = "e", digits = 1), if (rel < 1e-4) "  -> OK" else "  -> CHECK!")
    SOLVER_CHECKED <- TRUE
  }


  ## 9.4 Fit --------------------------------------------------------------------------------
  set.seed(RANDOM_SEED)
  say("")
  say("Fitting ", GRP, " (", N_FIT, " parameters: ", paste(FIT_NAMES, collapse = ", "), "; data transformation: ",
      switch(as.character(LAM), "0" = "log", "0.5" = "square root", "1" = "none", paste("power", LAM)), ")")
  FIT <- NULL
  T_FIT <- Sys.time()
  if (!FIT_PARAMETERS || N_FIT == 0) {
    say("  No fitting (FIT_PARAMETERS = FALSE or no parameter with fit = 1): model with start values.")
    pf <- to_fit_scale(P_FIXED[FIT_NAMES])
    FIT <- list(method = "none", best = pf, mll = OBJ(pf), sample = NULL, ci = NULL)
  } else if (FIT_METHOD == "parspace") {
    say("  Method: parameter-space explorer (roughness ", PS_ROUGHNESS,
        ", profiles ", if (PS_PROFILES) "on" else "off", ")")
    if (N_FIT > 4) say("  WARNING: the parameter-space explorer is not thoroughly tested for more than 4 free parameters!")
    res <- run_parspace(ps_settings(PS_ROUGHNESS), PS_PROFILES)
    FIT <- c(list(method = "parspace"), res)
  } else {
    say("  Method: Nelder-Mead simplex (", SIMPLEX_RUNS, " runs) from the start values")
    pf0 <- to_fit_scale(P_FIXED[FIT_NAMES])
    if (!is.finite(OBJ(pf0))) stop("The start values give an infinite MLL. Check start values and data.")
    pb <- pb_new(30 * N_FIT, "Simplex run 1 (rough)")
    r <- run_simplex(pf0, rough = TRUE, on_eval = function(i) pb_set(pb, i)); pb_close(pb)
    say(sprintf("  Simplex run 1 (rough): MLL = %s (%d evaluations)", fmt(r$value, 7), r$fevals))
    for (k in seq_len(max(0, SIMPLEX_RUNS - 1))) {
      pb <- pb_new(200 * N_FIT, sprintf("Simplex run %d", k + 1))
      r <- run_simplex(r$par, rough = FALSE, on_eval = function(i) pb_set(pb, i)); pb_close(pb)
      say(sprintf("  Simplex run %d: MLL = %s (%d evaluations%s)", k + 1, fmt(r$value, 7), r$fevals,
                  if (r$converged) "" else ", NOT converged - increase SIMPLEX_RUNS"))
    }
    FIT <- list(method = "simplex", best = r$par, mll = r$value, sample = NULL, ci = NULL)
    if (SIMPLEX_PROFILE_CI) {
      for (attempt in 1:3) {
        SP <- simplex_profiles(FIT$best, FIT$mll, "Profile likelihoods")
        if (is.null(SP$better)) break
        say("  Profiling found a better optimum (MLL ", fmt(SP$better[N_FIT + 1], 7), "); re-optimising ...")
        r <- run_simplex(SP$better[1:N_FIT], rough = FALSE)
        FIT$best <- r$par; FIT$mll <- r$value
      }
      ci <- matrix(NA_real_, N_FIT, 4, dimnames = list(FIT_NAMES, c("lo", "hi", "joint_lo", "joint_hi")))
      flags <- matrix(0, N_FIT, 3, dimnames = list(FIT_NAMES, c("lo_bound", "hi_bound", "broken")))
      for (i_p in seq_len(N_FIT)) {
        cr <- ci_from_profile(SP$profiles[[i_p]][, c(i_p, N_FIT + 1)], FIT$mll, 0.5 * 3.8415, BNDS[i_p, ])
        ci[i_p, 1:2] <- cr$ci; flags[i_p, ] <- cr$flags
      }
      FIT$ci <- ci; FIT$flags <- flags; FIT$profiles <- SP$profiles
      # parameter sets on the profiles inside the 95% region -> approximate CI bands
      smp <- sort_coll(do.call(rbind, SP$profiles))
      FIT$sample <- smp[smp[, N_FIT + 1] < FIT$mll + CHI_PROP, , drop = FALSE]
    }
  }
  FIT$fit_names <- FIT_NAMES
  FIT$time_needed <- elapsed(T_FIT)

  ## 9.5 Results: parameters, goodness of fit ---------------------------------------------
  P_BEST <- to_natural(FIT$best)
  DET    <- calc_mll(P_BEST, details = TRUE)
  MLL    <- DET$mll
  AIC    <- 2 * MLL + 2 * N_FIT
  RES <- data.frame(parameter = PAR$name, estimate = unname(P_BEST[PAR$name]),
                    CI95_lower = NA_real_, CI95_upper = NA_real_,
                    jointCI_lower = NA_real_, jointCI_upper = NA_real_,
                    fitted = PAR$fit, scale = PAR$scale, min = PAR$min, max = PAR$max,
                    remark = "", unit = UNITS[PAR$name], description = DESCR[PAR$name],
                    stringsAsFactors = FALSE)
  if (!is.null(FIT$ci)) for (nm in FIT_NAMES) {
    i <- match(nm, RES$parameter)
    RES[i, c("CI95_lower", "CI95_upper", "jointCI_lower", "jointCI_upper")] <- back(FIT$ci[nm, ], nm)
    rmk <- c(if (FIT$flags[nm, 1] == 1) "CI hits lower bound",
             if (FIT$flags[nm, 2] == 1) "CI hits upper bound",
             if (FIT$flags[nm, 3] == 1) "broken CI")
    RES$remark[i] <- paste(rmk, collapse = "; ")
  }
  if (!USE_BIOTRANSFORMATION) RES$remark[RES$parameter %in% c("km", "kem")] <- "not used (USE_BIOTRANSFORMATION = FALSE)"
  if (!BTP_UPTAKE_FROM_MEDIUM) RES$remark[RES$parameter == "kum"] <- "not used (BTP_UPTAKE_FROM_MEDIUM = FALSE)"
  rownames(RES) <- RES$parameter

  say("")
  rule(); say("Results ", GRP); rule()
  say("Method        : ", FIT$method, if (FIT$method == "parspace")
    sprintf("  (sample: %d sets in joint CI, %d in inner CI)", FIT$n_joint, FIT$n_inner) else "")
  say("Units         : concentrations in ", UNIT, ", time in days")
  say("Minus log-likelihood: ", fmt(MLL, 7), "   (AIC = ", fmt(AIC, 7), ")")
  for (st in FIT_STATES)
    say(sprintf("R-square %-14s: %6.3f  (n = %d, transformed data)", NAMES_ST[st],
                DET$r2[st], nrow(OBS[[st]])))
  say("")
  say("Best estimates and 95% CIs (profile likelihood)")
  say(strrep("-", 78))
  for (nm in PAR$name) {
    r <- RES[nm, ]
    if (r$fitted == 1) {
      say(sprintf("%-4s best: %10s (%10s - %-10s) fit: 1 (%s) %s", nm, fmt(r$estimate),
                  fmt(r$CI95_lower), fmt(r$CI95_upper), r$scale, if (nzchar(r$remark)) paste0(" * ", r$remark) else ""))
    } else if (!grepl("not used", r$remark)) {
      say(sprintf("%-4s fixed: %10s", nm, fmt(r$estimate)))
    }
  }
  say(strrep("-", 78))
  if (FIT$method == "simplex" && !SIMPLEX_PROFILE_CI) say("(No CIs: set SIMPLEX_PROFILE_CI = TRUE or use FIT_METHOD = \"parspace\")")
  if (any(grepl("bound", RES$remark))) say("* A CI runs into a search bound: consider widening min/max (SECTION 1.6).")


  ## 9.6 Derived quantities incl. BCF, with CIs --------------------------------------------
  DER <- data.frame(quantity = names(derived(P_BEST)), estimate = derived(P_BEST),
                    CI_lower = NA_real_, CI_upper = NA_real_, stringsAsFactors = FALSE)
  CI_SET <- NULL
  if (!is.null(FIT$sample) && nrow(FIT$sample) > 1) {
    CI_SET <- FIT$sample[FIT$sample[, L_COL] < MLL + CHI_PROP, , drop = FALSE]
    if (nrow(CI_SET) > CI_MAX_SETS)
      CI_SET <- CI_SET[c(1, sample(2:nrow(CI_SET), CI_MAX_SETS - 1)), , drop = FALSE]
    dmat <- t(apply(CI_SET[, 1:N_FIT, drop = FALSE], 1, function(pf) derived(to_natural(pf))))
    DER$CI_lower <- apply(dmat, 2, min); DER$CI_upper <- apply(dmat, 2, max)
  }
  DER$meaning <- DER_NOTE[DER$quantity]
  rownames(DER) <- DER$quantity
  say("")
  say("Derived quantities (CIs: range over the parameter sets within the 95% region)")
  for (i in seq_len(nrow(DER)))
    say(sprintf("  %-24s %10s (%10s - %-10s)", DER$quantity[i], fmt(DER$estimate[i]),
                fmt(DER$CI_lower[i]), fmt(DER$CI_upper[i])))
  if (BTP_UPTAKE_FROM_MEDIUM) say("  (Note: with BTP uptake from the medium the steady-state ratio km/kem is incomplete.)")

  ## 9.7 Model curves, 95% CI bands and 95% prediction intervals ------------------------
  # CI band  : envelope of the model curves of all parameter sets in the 95% region
  #            (BYOM calc_conf type 3); shown as shaded band in the fit plot
  # PI band  : same envelope, each curve widened by +-1.96 residual SD on the
  #            transformed scale (parameter uncertainty + residual error);
  #            only used for the PPC metric
  SIG <- DET$sigma
  CURVES <- list(); G_CURVE <- list()
  for (sc in SCEN_IDS) {
    pp <- EXPO[[sc]]$parent
    t_max <- max(c(GRID_DATA[[sc]]$t_out, pp$time))
    t_c <- sort(unique(c(seq(0, t_max, length.out = N_TIME_POINTS_CURVE), pp$time,
                         GRID_DATA[[sc]]$t_out)))
    G_CURVE[[sc]] <- make_grid(EXPO[[sc]], t_c, GRID_DATA[[sc]]$x0)
    X <- solve_tk(P_BEST, G_CURVE[[sc]])
    CURVES[[sc]] <- list(t = t_c, best = cbind(X, total = X[, 1] + X[, 2]), lo = NULL, hi = NULL)
  }
  pi_update <- function(lo, hi, X) {           # X: model on natural scale (parent, btp)
    for (st in FIT_STATES) {
      tm <- transform_data(pmax(X[, st], 0))
      lo[, st] <- pmin(lo[, st], tm - Z95 * SIG[[st]])
      hi[, st] <- pmax(hi[, st], tm + Z95 * SIG[[st]])
    }
    list(lo = lo, hi = hi)
  }
  init_pi <- function(n) { m <- matrix(Inf, n, 2, dimnames = list(NULL, c("parent", "btp"))); list(lo = m, hi = -m) }
  PI_D <- pi_update(init_pi(nrow(DET$pred))$lo, init_pi(nrow(DET$pred))$hi, DET$pred)   # at the data
  PI_C <- lapply(SCEN_IDS, function(sc) { z <- init_pi(length(CURVES[[sc]]$t)); pi_update(z$lo, z$hi, CURVES[[sc]]$best) })
  names(PI_C) <- SCEN_IDS
  if (!is.null(CI_SET) && nrow(CI_SET) > 1) {
    say("")
    say("Calculating 95% CI and prediction bands from ", nrow(CI_SET), " parameter sets ...")
    for (sc in SCEN_IDS) { CURVES[[sc]]$lo <- CURVES[[sc]]$best; CURVES[[sc]]$hi <- CURVES[[sc]]$best }
    pb <- pb_new(nrow(CI_SET), "CI and prediction bands")
    for (k in seq_len(nrow(CI_SET))) {
      p <- to_natural(CI_SET[k, 1:N_FIT])
      for (sc in SCEN_IDS) {
        X <- solve_tk(p, G_CURVE[[sc]]); X <- cbind(X, total = X[, 1] + X[, 2])
        CURVES[[sc]]$lo <- pmin(CURVES[[sc]]$lo, X)
        CURVES[[sc]]$hi <- pmax(CURVES[[sc]]$hi, X)
        PI_C[[sc]] <- pi_update(PI_C[[sc]]$lo, PI_C[[sc]]$hi, X)
      }
      PI_D <- pi_update(PI_D$lo, PI_D$hi, predict_data(p))
      pb_set(pb, k)
    }
    pb_close(pb)
  } else {
    say("")
    say("Prediction interval from the best fit only (no parameter sample available).")
  }
  PI_LO <- PI_D$lo; PI_HI <- PI_D$hi                 # transformed scale, at the data
  for (sc in SCEN_IDS) {                              # natural scale, for plots/tables
    CURVES[[sc]]$pi_lo <- back_transform(PI_C[[sc]]$lo)
    CURVES[[sc]]$pi_hi <- back_transform(PI_C[[sc]]$hi)
    CURVES[[sc]]$pi_lo[!is.finite(PI_C[[sc]]$lo)] <- NA
    CURVES[[sc]]$pi_hi[!is.finite(PI_C[[sc]]$hi)] <- NA
  }

  ## 9.8 Goodness-of-fit metrics ------------------------------------------------------
  MET <- calc_metrics()
  say("")
  say(sprintf("Goodness of fit (PPC acceptable >= %g%%, NRMSE acceptable <= %g%%, |PE| acceptable <= %g%%; PE > 0 = overprediction)",
              PPC_MIN, NRMSE_MAX, PE_MAX))
  say(sprintf("  %-22s %-26s %4s %7s %9s %9s %9s", "state", "scenario", "n", "R2", "NRMSE%", "PPC%", "PE%"))
  for (i in seq_len(nrow(MET)))
    say(sprintf("  %-22s %-26s %4d %7.3f %8s%s %8s%s %8s%s", MET$state[i], MET$scenario[i], MET$n[i], MET$R2[i],
                fmt(MET$NRMSE_pct[i], 3), if (isTRUE(MET$NRMSE_ok[i])) " " else "!",
                fmt(MET$PPC_pct[i], 3), if (isTRUE(MET$PPC_ok[i])) " " else "!",
                fmt(MET$PE_pct[i], 3), if (isTRUE(MET$PE_ok[i])) " " else "!"))
  say("  (! = outside the acceptance criterion; PE at the end of exposure, 'all' = mean over scenarios, ok only if every scenario is ok)")

  ## 9.9 Plots ------------------------------------------------------------------------------
  show_and_save(draw_fits, "medium_tissue_concentrations", width = 10.5,
                height = 4.1 * min(length(SCEN_IDS), PER_PAGE) + 0.9,
                n_pages = ceiling(length(SCEN_IDS) / PER_PAGE))
  if (length(CURVE_SPEC) > 1)                 # extra file: each compartment in its own panel
    show_and_save(draw_compartments, "compartment_concentrations",
                  width = 4.3 * length(CURVE_SPEC), height = 4.6)
  if (!is.null(FIT$sample) && nrow(FIT$sample) > 1 && N_FIT >= 1)
    show_and_save(draw_parspace, "parameter_space", width = max(6.5, 2.4 * N_FIT + 1.5), height = max(6.5, 2.4 * N_FIT + 1.5))

  say("Time required for ", GRP, ": ", FIT$time_needed)

  ## 9.10 Row for the summary over all model runs -----------------------------------------
  row <- data.frame(model_run = GRP, version = VERSION, scenarios = paste(SCEN_IDS, collapse = " + "),
                    biotransformation = USE_BIOTRANSFORMATION, method = FIT$method,
                    MLL = MLL, AIC = AIC, stringsAsFactors = FALSE)
  for (nm in FIT_NAMES) {
    row[[nm]] <- RES[nm, "estimate"]
    row[[paste0(nm, "_CI95_lower")]] <- RES[nm, "CI95_lower"]
    row[[paste0(nm, "_CI95_upper")]] <- RES[nm, "CI95_upper"]
  }
  row$BCF_parent <- DER["BCF_parent", "estimate"]
  row$BCF_parent_CI_lower <- DER["BCF_parent", "CI_lower"]
  row$BCF_parent_CI_upper <- DER["BCF_parent", "CI_upper"]
  for (st in FIT_STATES) {
    m <- MET[MET$state == NAMES_ST[[st]] & MET$scenario == "all (pooled)", ]
    row[[paste0("R2_", st)]]        <- m$R2
    row[[paste0("PPC_pct_", st)]]   <- m$PPC_pct
    row[[paste0("NRMSE_pct_", st)]] <- m$NRMSE_pct
    row[[paste0("PE_pct_", st)]]    <- m$PE_pct
  }
  row$all_criteria_met <- all(MET$PPC_ok[MET$scenario == "all (pooled)"]) &&
                          all(MET$NRMSE_ok[MET$scenario == "all (pooled)"]) &&
                          all(MET$PE_ok[MET$scenario == "all (pooled)"])
  row$remarks <- paste(paste0(RES$parameter, ": ", RES$remark)[nzchar(RES$remark) & RES$fitted == 1], collapse = "; ")
  SUMMARY[[JOB_ID]] <- row

  ## 9.11 Collect the results for the Excel file of the whole run (SECTION 10) ----------
  # Sheet 1: rate constants and derived quantities (BCF, half-life, ...) of
  # every model, each with the best estimate and the 95% confidence interval.
  used <- !grepl("not used", RES$remark)
  TAB <- rbind(
    data.frame(parameter = RES$parameter[used], description = RES$description[used],
               estimate = RES$estimate[used], CI95_lower = RES$CI95_lower[used], CI95_upper = RES$CI95_upper[used],
               unit = RES$unit[used],
               note = ifelse(RES$fitted[used] == 1, RES$remark[used], "fixed (not fitted)"),
               stringsAsFactors = FALSE),
    data.frame(parameter = DER$quantity, description = trimws(sub("\\[.*$", "", DER$meaning)),
               estimate = DER$estimate, CI95_lower = DER$CI_lower, CI95_upper = DER$CI_upper,
               unit = sub("^.*\\[(.*)\\]$", "\\1", DER$meaning),
               note = "derived from the rate constants", stringsAsFactors = FALSE))
  TAB <- TAB[, c(setdiff(names(TAB), "description"), "description")]   # description last
  ALL_TAB[[JOB_ID]] <- cbind(model = JOB_ID, TAB)
  # Sheet 3: observations and model curves (as in the fit plot, in the units of
  # the data file), so the fits can be re-plotted with your own script.
  for (sc in SCEN_IDS) {
    C <- CURVES[[sc]]
    for (s in CURVE_SPEC) {
      FIT_DATA[[length(FIT_DATA) + 1]] <- data.frame(
        model = JOB_ID, scenario = sc, series = s$label, type = "model curve",
        time = C$t, value = C$best[, s$col] * s$f,
        CI95_lower = if (s$band && !is.null(C$lo)) C$lo[, s$col] * s$f else NA_real_,
        CI95_upper = if (s$band && !is.null(C$hi)) C$hi[, s$col] * s$f else NA_real_,
        stringsAsFactors = FALSE)
      if (!is.na(s$obs)) {
        o <- OBS[[s$obs]][OBS[[s$obs]]$scenario == sc, ]
        if (nrow(o)) FIT_DATA[[length(FIT_DATA) + 1]] <- data.frame(
          model = JOB_ID, scenario = sc, series = s$label, type = "observation",
          time = o$time, value = o$value * s$f,
          CI95_lower = NA_real_, CI95_upper = NA_real_, stringsAsFactors = FALSE)
      }
    }
  }

  ## 9.12 Log of this model -------------------------------------------------------------
  writeLines(c(.LOG$lines[seq_len(HEADER_END)], .LOG$lines[LOG_FROM:length(.LOG$lines)]), out_file("log.txt"))
}


###############################################################################
## SECTION 10 : Summary over all model runs
###############################################################################

TRT_TAG <- ""
say("")
rule("#"); say("SUMMARY OVER ALL MODEL RUNS"); rule("#")
if (length(SUMMARY)) {
  all_cols <- unique(unlist(lapply(SUMMARY, names)))
  SUM_TAB <- do.call(rbind, lapply(SUMMARY, function(r) { r[setdiff(all_cols, names(r))] <- NA; r[all_cols] }))
  for (i in seq_len(nrow(SUM_TAB))) {
    say(sprintf("%s   (%s)", SUM_TAB$model_run[i], SUM_TAB$scenarios[i]))
    say(sprintf("    MLL = %s, AIC = %s", fmt(SUM_TAB$MLL[i], 6), fmt(SUM_TAB$AIC[i], 5)))
    say("    version: ", if (SUM_TAB$biotransformation[i]) "with biotransformation" else "without biotransformation")
    for (nm in intersect(c("ku", "ke", "km", "kem", "kum"), names(SUM_TAB))) if (!is.na(SUM_TAB[[nm]][i]))
      say(sprintf("    %-10s %10s (%10s - %-10s)", nm, fmt(SUM_TAB[[nm]][i]),
                  fmt(SUM_TAB[[paste0(nm, "_CI95_lower")]][i]), fmt(SUM_TAB[[paste0(nm, "_CI95_upper")]][i])))
    say(sprintf("    %-10s %10s (%10s - %-10s)", "BCF", fmt(SUM_TAB$BCF_parent[i]),
                fmt(SUM_TAB$BCF_parent_CI_lower[i]), fmt(SUM_TAB$BCF_parent_CI_upper[i])))
    for (st in c("parent", "btp")) if (!is.null(SUM_TAB[[paste0("R2_", st)]]) && !is.na(SUM_TAB[[paste0("R2_", st)]][i]))
      say(sprintf("    %-10s R2 = %s, PPC = %s%%, NRMSE = %s%%, PE = %s%%", NAMES_ST[[st]],
                  fmt(SUM_TAB[[paste0("R2_", st)]][i], 3), fmt(SUM_TAB[[paste0("PPC_pct_", st)]][i], 3),
                  fmt(SUM_TAB[[paste0("NRMSE_pct_", st)]][i], 3), fmt(SUM_TAB[[paste0("PE_pct_", st)]][i], 3)))
  }
  ## 10.1 ONE Excel file for the whole run --------------------------------------------
  qual_cols <- intersect(c("model_run", "version", "model", "biotransformation", "scenarios", "method", "MLL", "AIC"),
                         names(SUM_TAB))
  QUAL <- SUM_TAB[, c(qual_cols, grep("^(R2|PPC_pct|NRMSE_pct|PE_pct)_", names(SUM_TAB), value = TRUE)), drop = FALSE]
  PTAB <- do.call(rbind, ALL_TAB)
  rownames(QUAL) <- rownames(PTAB) <- NULL
  XLSX_FILE <- file.path(OUTPUT_DIR, paste0(RUN_NAME, "_", RUN_STAMP, "_TK_results.xlsx"))
  write_excel(list(TK_parameters = PTAB, quality_of_fit = QUAL,
                   data_and_curves = do.call(rbind, FIT_DATA)), XLSX_FILE)
  say("")
  say("Excel results of ALL models in one file: ", XLSX_FILE)
  say("  'TK_parameters'  : fitted and derived parameters with 95% CIs")
  say("  'quality_of_fit' : MLL, AIC, R2, PPC, NRMSE, PE per model")
  say("  'data_and_curves': observations and model curves (best fit + 95% CI) in the")
  say("                     units of the data file, as input for your own plots")
} else say("Nothing could be fitted.")
if (length(SKIPPED))
  say("Skipped: ", paste(paste0(names(SKIPPED), " (", SKIPPED, ")"), collapse = ", "))
rule()
say("Total time: ", elapsed(T_START))
say("Results in: ", normalizePath(OUTPUT_DIR))
rule()

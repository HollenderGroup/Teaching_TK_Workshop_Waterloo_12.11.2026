# Model mathematics for dummies

*How the TK workshop scripts work, in plain language.*

The script does four things in order: it describes the animal with an equation, predicts tissue concentrations from that equation, finds the parameter values that best match your data, and works out how certain those values are.

**1. The model: a bucket with a tap and a drain**

Think of the organism as a bucket. The chemical flows in from the water and leaks out:

change per day = **ku × water concentration** − **ke × tissue concentration**

The inflow (ku) depends on how much is in the water. The outflow (ke) depends on how much is already in the animal: the fuller the bucket, the faster it leaks. Eventually inflow equals outflow, and the level stops rising. That is steady state, with BCF = ku/ke. The other scripts just add more buckets (lipid, receptor, metabolite) or more taps (food).

**2. Predicting the tissue curve**

Exposure concentrations are only measured at a few time points, so the script joins them with straight lines. On each straight segment, the bucket equation has an exact textbook solution: an exponential approach towards the current water level. The script steps from one segment to the next, using the end of one as the start of the next. This is fast and exact. The one exception is receptor binding: its saturation term (free receptors run out) has no exact formula, so a numerical solver takes many small steps instead.

**3. Finding the best parameters**

For a given ku and ke, the script compares the predicted curve with every measured tissue value. It first takes the square root of both, so high and low concentrations count more evenly. It then adds up the squared differences. A small sum means a good fit. This sum is turned into a “minus log-likelihood” (MLL), which is just a rescaled version of it: lower is better. The task is to find the ku and ke with the lowest MLL.

**4. The parameter-space explorer: how sure are we?**

Instead of searching only for the single best point, the script tests thousands of random combinations of ku and ke within your min–max ranges. It then samples more densely around the good ones. Each combination gets its MLL. That gives the two kinds of points in the plot:

- **Green**: combinations that are almost as good as the best one, i.e. within 1.92 MLL units of it. 1.92 is the statistical cut-off for 95% confidence on one parameter.

- **Cyan**: a looser cut-off for the 95% region of all parameters together.

The **red line** is the profile likelihood. For each value of ku, the script finds the best possible ke and records how much worse that fit is than the overall best. Where the red line crosses the 1.92 line, you read off the 95% confidence interval. A narrow V means the data pin the parameter down. A flat line means they don’t. A tilted cloud (like ku vs ke) means the two trade off against each other: a higher uptake rate can be compensated by faster elimination.

**5. The shaded bands and the comparison numbers**

The **shaded band** in the fit plot is made by drawing the model curve for every green parameter set and shading the range between the highest and lowest curve.

**AIC** = 2 × MLL + 2 × number of fitted parameters. It rewards a good fit but penalizes extra parameters. That makes it fair to compare a one-compartment model with a two-compartment model on the same data: the extra compartment has to improve the fit enough to pay for its extra parameters.

**Modifying the script manually:**

To add a pathway, you edit six places in TK_1_one_compartment.R. The example below adds uptake from sediment at a fixed concentration:

**dC/dt = ku·Cw(t) + ks·Csed − ke·C**

The line numbers are approximate.

**1. Parameter table (SECTION 1.6, ~line 235).** Add one row per new parameter:

```r
ku 100 1 1 5000 log
ke 1 1 0.01 100 log
ks 0.1 1 0.0001 100 log
Csed 0.05 0 0.001 10 log
```

Csed has fit = 0, so it stays fixed at its start value. Enter it in model units: with USE_MOLAR_UNITS = TRUE that is µmol/kg. The script only converts the data, not this table.

**2. Allowed names (SECTION 3, check_settings, ~line 410).**

```r
need <- c("ku", "ke", "ks", "Csed")
```

Without this, the script stops with “PARAMETERS must contain exactly these names”.

**3. Solver (SECTION 1.7, ~line 277).**

```r
SOLVER <- "ode"
```

The exact solution in 5.3 only works for the original equation, so switch to the numerical solver. This needs the deSolve package. The solver check in 9.3 will then print -\> CHECK!, which is expected: it compares your new model with the old exact solution.

**4. The equation (SECTION 5.1, ~line 790).** This is the only place the model is defined:

```r
dC <- ku * Cw + ks * Csed - ke * C # water + sediment - elimination
```

Inside with(as.list(parms), ...) every parameter is available by its name.

**5. Descriptions (SECTION 6.1, DESCR, ~line 934).** These text labels go into the Excel sheet:

```r
DESCR <<- c(ku = "uptake rate constant",
ke = "elimination rate constant",
ks = "uptake rate constant from sediment",
Csed = "sediment concentration")
```

**6. Units and derived values (SECTION 8, ~line 1552).**

```r
UNITS <- c(ku = "L/kg ww/d", ke = "1/d", ks = "kg dw/kg ww/d", Csed = "umol/kg dw")
```

If you want a new derived value, such as a sediment accumulation factor, add it in two places. In derived(), add BSAF = p\[\["ks"\]\] / p\[\["ke"\]\]. In DER_NOTE, add a matching line: BSAF = "ks/ke: ... \[kg dw/kg ww\]". The unit goes at the end in square brackets, because the Excel output reads it from there.

**If the new source changes over time and is measured** (for example food or sediment concentrations from the data file), the fixed Csed parameter isn’t enough. The script also has to read those rows and interpolate them, like it does for the medium. TK_4_dietary_uptake.R is exactly that case. To see what changed there, compare these parts with script 1:

- FOOD_MATRIX in SECTION 1.3

- the FOOD rows in SECTION 4

- make_grid in 5.2, which adds af and bf

- solve_ode in 5.4, which passes cf = cf

To check that your edit worked, run with FIT_METHOD \<- "simplex" first, which is quick. The new parameter should appear in the console table with an estimate. Once the fit plot looks sensible, switch back to "parspace".

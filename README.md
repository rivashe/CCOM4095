# MicroEven

**A Shiny app for computing microbiome evenness metrics from MicrobiomeAnalyst outputs**

> **Author:** Eliud R. Rivas Hernandez — Colorado State University / University of Puerto Rico, Arecibo  
> **Advisors:** Dr. Tiffany Weir (CSU) · Prof. Emilio Perez Arnau (UPRA)

---

## What is MicroEven?

MicrobiomeAnalyst exports Shannon diversity and Observed Features (richness) from its alpha-diversity module, but it does **not** compute evenness. MicroEven fills that gap. Upload your two MicrobiomeAnalyst CSV files, and the app will:

- Compute **Pielou's evenness** (J′ = H′ / ln S) and **Shannon equitability** (exp(H′) / S)
- Display per-sample and group-level summary tables
- Generate boxplots, violin plots, and scatter plots
- Run the appropriate statistical test for your study design (linear mixed model, Wilcoxon, Kruskal-Wallis, or ANOVA/t-test)
- Export results and figures as CSV, PNG, and PDF

---

## Quickstart — launch the app in one step

You do not need to download any files. Run this in an R session:

```r
install.packages("shiny")   # skip if shiny is already installed
library(shiny)
runGitHub(repo = "CCOM4095", username = "rivashe", ref = "main")
```

Shiny will download the app, install any missing packages, and open it in your browser automatically.

---

## Local installation (optional)

If you prefer to run the app from a local copy:

**1. Install dependencies**

```r
install.packages(c(
  "shiny", "dplyr", "tidyr", "readr", "ggplot2", "DT", "bslib",
  "lme4", "lmerTest", "emmeans", "broom.mixed", "FSA", "car"
))
```

**2. Download and run**

```r
library(shiny)
runGitHub(repo = "CCOM4095", username = "rivashe", ref = "main")
```

Or, if you have cloned the repository:

```r
shiny::runApp("app.R")
```

You can also open `app.R` in RStudio and click **Run App**.

---

## Input files

MicroEven requires **two CSV files** exported from the MicrobiomeAnalyst alpha-diversity module: one for **Observed Features** and one for **Shannon diversity**. These are the files MicrobiomeAnalyst gives you when you download alpha-diversity results — no reformatting is needed.

### Required columns

Both files must contain these three columns (exact names, case-sensitive):

| Column | What it contains |
|---|---|
| `sample_id` | A unique identifier for each sample. Used to join the two files. |
| `variable` | The metric name. Must contain `"Observed"` in the Observed Features file and `"Shannon"` in the Shannon file. |
| `value` | The numeric metric value for that sample. |

### Optional metadata columns

Any additional columns in your files (treatment group, week, time point, subject ID, cage, visit, etc.) are automatically detected and preserved. The app will offer them as grouping, faceting, and subject-ID options in the sidebar.

### What the CSV should look like

**Observed Features file (`alpha_diversity_observed.csv`)**

```
sample_id,variable,value,group,week,subject_id
S01,Observed,142,Control,0,P01
S02,Observed,158,Control,0,P02
S03,Observed,134,FLAX,0,P03
S04,Observed,167,FLAX,0,P04
S05,Observed,139,Control,6,P01
S06,Observed,161,Control,6,P02
S07,Observed,128,FLAX,6,P03
S08,Observed,175,FLAX,6,P04
```

**Shannon file (`alpha_diversity_shannon.csv`)**

```
sample_id,variable,value,group,week,subject_id
S01,Shannon,4.21,Control,0,P01
S02,Shannon,4.53,Control,0,P02
S03,Shannon,4.07,FLAX,0,P03
S04,Shannon,4.68,FLAX,0,P04
S05,Shannon,4.19,Control,6,P01
S06,Shannon,4.59,Control,6,P02
S07,Shannon,4.01,FLAX,6,P03
S08,Shannon,4.75,FLAX,6,P04
```

> **Tip:** If MicrobiomeAnalyst added an unnamed row-number column (it shows up as a blank header or `...1`), MicroEven will detect and remove it automatically.

### Common errors and fixes

| Error message | What it means | How to fix |
|---|---|---|
| `File is missing required column(s): sample_id` | The column names don't match exactly | Check capitalization — must be `sample_id`, `variable`, `value` (all lowercase) |
| `Expected metric 'Observed' not found in file` | You uploaded the Shannon file where the Observed file is expected, or vice versa | Swap the files in the upload boxes |
| `Expected metric 'Shannon' not found in file` | Same issue, opposite file | Upload the correct file to each box |
| Samples joined: 0 | The `sample_id` values don't match between the two files | Make sure both files use the same sample identifiers |

---

## How to use the app

1. **Upload files** — Use the sidebar to upload your Observed Features CSV and Shannon CSV.
2. **Set grouping** — Choose a metadata column (e.g., `group`, `treatment`) as the grouping variable for plots and statistics.
3. **Optional: set faceting** — Choose a second variable (e.g., `week`) to split the plot into panels.
4. **Optional: set subject ID** — If the same subjects were measured at multiple time points, select the subject column. The app will then use a linear mixed model instead of assuming independence.
5. **Explore the tabs** — Per-sample table, group summary, plot, and statistics are each on their own tab.
6. **Download** — Export tables as CSV and figures as PNG or PDF from the sidebar.

---

## Statistical tests

The app picks the right test automatically based on your inputs, or you can override it manually.

| Situation | Test used |
|---|---|
| Repeated measures (subject ID set) | Linear mixed model with random intercept per subject (lmerTest, Satterthwaite df); pairwise contrasts via emmeans with Tukey adjustment |
| Two independent groups | Wilcoxon rank-sum (non-parametric, safe at small n) |
| Three or more independent groups | Kruskal-Wallis with Dunn's post-hoc (Benjamini-Hochberg FDR) |
| Parametric (manual override) | Welch t-test (2 groups) or one-way ANOVA with Tukey HSD (3+ groups) |

Residual diagnostics (Shapiro-Wilk test + QQ plot) are shown for LMM and parametric tests so you can check assumptions.

---

## Metrics computed

**Pielou's evenness (J′)**  
J′ = H′ / ln(S), where H′ is the Shannon index and S is the number of observed features. Ranges from 0 (dominated by one taxon) to 1 (perfectly even). This is the standard evenness measure reported in most microbiome papers.

**Shannon equitability**  
exp(H′) / S. A bounded approximation of evenness derived from Shannon diversity. Note: this is *not* the same as Simpson evenness, which requires the Simpson index. The app labels it clearly to avoid misinterpretation.

---

## Citation

If MicroEven contributed to a thesis, poster, or publication, please cite it as:

> Rivas Hernandez, E. R. (2026). *MicroEven: a Shiny application for computing microbiome evenness metrics and group-level statistical comparisons.* Colorado State University.

The app's **About / Help** tab contains BibTeX and full references for all methods and packages.

---

## Contact

**Eliud R. Rivas Hernandez** · [eliud.rivas@upr.edu](mailto:eliud.rivas@upr.edu)

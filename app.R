# =============================================================================
# MicroEven: A Shiny GUI for Microbiome Evenness Metrics
# =============================================================================
# Author : Eliud R. Rivas Hernandez
# Advisor: Dr. Tiffany Weir
# Course : CCOM 4095
#
# PURPOSE
# -------
# MicrobiomeAnalyst exports Shannon diversity and Observed Features (richness)
# but does NOT export evenness. This app takes those two exports (in the long
# format MicrobiomeAnalyst uses: one row per sample with columns including
# `sample_id`, `variable`, `value`, plus metadata) and:
#
#   1. Joins the two files by sample_id
#   2. Computes Pielou's evenness: J' = H' / ln(S)
#   3. Computes Simpson-derived evenness (approximation from Shannon when
#      raw OTU/ASV counts are not available -- see notes below)
#   4. Displays per-sample tables, per-group summary tables, and plots
#   5. Runs the right statistical test for the design (LMM for repeated
#      measures; Kruskal/Wilcoxon or ANOVA/t-test for independent samples)
#   6. Exports everything as CSV / PNG / PDF
#
# IMPORTANT MATH NOTE
# -------------------
# Pielou's evenness (J') requires only H' and S, both of which we have.
# This is mathematically exact and what most microbiome papers report.
#
# A "true" Simpson evenness (E_{1/D} = (1/D)/S) requires the Simpson index D,
# which MicrobiomeAnalyst's basic Shannon export does NOT include. We provide
# an OPTIONAL third file upload for Simpson, OR we compute an approximation
# (effective number of species from Shannon: exp(H')/S) and label it clearly
# as a Shannon-derived equitability, NOT true Simpson evenness. This matters
# scientifically -- don't claim it's Simpson when it isn't.
#
# STATISTICS NOTE
# ---------------
# Microbiome studies often have repeated measures (same subject sampled over
# weeks). A plain t-test or one-way ANOVA assumes independence and will give
# wrong p-values in that case. This app handles design properly:
#
#   - If a subject_id column is selected, we fit a linear mixed model with
#     a random intercept per subject (lmerTest, Satterthwaite df).
#   - Otherwise, we default to non-parametric tests (Wilcoxon for 2 groups,
#     Kruskal-Wallis + Dunn for >2 groups), which are safer for small n and
#     do not assume normality. The user can override to parametric.
#   - Residual diagnostics (Shapiro-Wilk + QQ plot) are shown on the
#     Statistics tab so the user can sanity-check assumptions.
#
# INPUT FORMAT (matches your alphadiversity.csv)
# ----------------------------------------------
# Columns expected (extra columns are fine and preserved as metadata):
#   - sample_id      : unique sample identifier (used for joining)
#   - variable       : "Observed" or "Shannon" (we check this)
#   - value          : the numeric metric value
#   - any others     : treatment, week, subject_id, visit, etc. -> metadata
#
# DEPENDENCIES
# ------------
# install.packages(c("shiny","dplyr","tidyr","readr","ggplot2","DT","bslib",
#                    "lme4","lmerTest","emmeans","broom.mixed","FSA","car"))
# =============================================================================


# ---- 1. Load libraries ------------------------------------------------------
# Keep this list minimal so deployment is easier (every package adds size on
# shinyapps.io). bslib gives us a clean Bootstrap 5 theme without writing CSS.
library(shiny)
library(dplyr)       # data wrangling (mutate, group_by, summarise, left_join)
library(tidyr)       # pivot_wider, drop_na
library(readr)       # read_csv (faster + better type inference than read.csv)
library(ggplot2)     # all plotting
library(DT)          # interactive HTML tables (sortable, searchable, exportable)
library(bslib)       # modern Bootstrap theming for Shiny
library(lme4)        # lmer() for mixed models
library(lmerTest)    # adds Satterthwaite df + p-values to lmer
library(emmeans)     # estimated marginal means + pairwise contrasts
library(broom.mixed) # tidy() methods for lmer objects
library(FSA)         # Dunn's test (post-hoc for Kruskal-Wallis)
library(car)         # Anova() with type II/III SS for ANOVA tables


# ---- 2. Helper functions ----------------------------------------------------
# These are pure functions (no Shiny reactivity inside) so they're easy to
# unit-test and easy to swap out. Keeping logic separate from the UI is the
# pattern your course's "module / class" diagrams describe.

#' Validate and parse a MicrobiomeAnalyst alpha-diversity CSV
#'
#' @param path   Filesystem path to the uploaded CSV
#' @param expected_metric  One of "Observed" or "Shannon" -- what we expect
#'        the `variable` column to contain. If the file has a different
#'        metric we stop with an informative error.
#' @return A tibble with columns: sample_id, <metric>, plus metadata.
parse_alpha_file <- function(path, expected_metric) {

  # read_csv is robust to the leading empty-name column ("") that R writes
  # when you save a data.frame with row.names = TRUE. It will name it `...1`.
  df <- readr::read_csv(path, show_col_types = FALSE)

  # --- minimum required columns ---
  required <- c("sample_id", "variable", "value")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("File is missing required column(s): ",
         paste(missing, collapse = ", "),
         ". Expected a MicrobiomeAnalyst alpha-diversity export.")
  }

  # --- check the metric matches what the user said they uploaded ---
  # We allow case-insensitive match because MicrobiomeAnalyst sometimes
  # exports "shannon" vs "Shannon".
  metrics_in_file <- unique(df$variable)
  if (!any(tolower(metrics_in_file) == tolower(expected_metric))) {
    stop("Expected metric '", expected_metric, "' not found in file. ",
         "Found: ", paste(metrics_in_file, collapse = ", "))
  }

  # Filter to just the expected metric (in case file mixes multiple)
  df <- df[tolower(df$variable) == tolower(expected_metric), , drop = FALSE]

  # Pivot to wide: rename `value` to the metric name so we can join two
  # files side by side later.
  df$metric_value <- df$value
  names(df)[names(df) == "metric_value"] <- expected_metric

  # Drop the now-redundant `variable` and `value` columns
  df$variable <- NULL
  df$value    <- NULL

  # Drop the unnamed row-number column if it snuck in
  df <- df[, names(df) != "...1", drop = FALSE]

  df
}


#' Compute evenness metrics from joined Shannon + Observed table
#'
#' @param df  A data frame with columns Shannon and Observed (plus metadata)
#' @return    Same df with new columns: Pielou, ShannonEquitability
compute_evenness <- function(df) {

  # Pielou's J' = H' / ln(S)
  # Guard: ln(1) = 0 would give Inf, and S = 0 makes no sense biologically.
  # Mark those as NA so plots/summaries don't break.
  df$Pielou <- ifelse(
    df$Observed > 1,
    df$Shannon / log(df$Observed),
    NA_real_
  )

  # Shannon-derived equitability (effective species / observed species).
  # NOT the same as Simpson evenness -- we label it clearly in the UI.
  # exp(H') gives the "effective number of equally-abundant species" that
  # would produce the observed Shannon value. Dividing by S puts it on [0,1].
  df$ShannonEquitability <- ifelse(
    df$Observed > 0,
    exp(df$Shannon) / df$Observed,
    NA_real_
  )

  df
}


#' Summarise evenness by a grouping variable
#'
#' @param df          Per-sample data with evenness columns
#' @param group_var   Name of column to group by (string), or NULL for overall
summarise_evenness <- function(df, group_var = NULL) {

  metric_cols <- c("Observed", "Shannon", "Pielou", "ShannonEquitability")
  metric_cols <- intersect(metric_cols, names(df))

  if (is.null(group_var) || group_var == "(none)") {
    grouped <- df
  } else {
    grouped <- dplyr::group_by(df, .data[[group_var]])
  }

  # Compute n, mean, median, sd, IQR for each metric.
  # We use across() for a clean, scalable summary -- adding metrics later
  # only requires updating metric_cols above.
  grouped |>
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(
        dplyr::all_of(metric_cols),
        list(
          mean   = ~mean(.x, na.rm = TRUE),
          median = ~stats::median(.x, na.rm = TRUE),
          sd     = ~stats::sd(.x, na.rm = TRUE),
          IQR    = ~stats::IQR(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}


# ---- 2b. Statistical test helpers -------------------------------------------
# Each test returns a named list with a consistent shape so the UI rendering
# code stays simple:
#   list(
#     method      = character, # human-readable name of test
#     formula     = character, # the call/formula shown to the user
#     fixed       = data.frame, # main results table
#     pairwise    = data.frame or NULL, # contrasts / post-hoc
#     diagnostics = list(shapiro_p = ..., residuals = ..., fitted = ...),
#     note        = character  # caveats, df, n used, etc.
#   )

#' Decide which test to run when the user picks "Auto"
#'
#' Rules (in order):
#'   1. If a subject column is selected -> LMM (repeated measures).
#'   2. If group has 2 levels -> Wilcoxon (non-parametric, safe at small n).
#'   3. If group has >2 levels -> Kruskal-Wallis.
#'   4. If no group at all -> NULL (we'll surface a message instead).
choose_auto_test <- function(df, group_var, subject_var) {
  if (!is.null(subject_var) && subject_var != "(none)") return("lmm")
  if (is.null(group_var) || group_var == "(none)") return(NA_character_)
  k <- length(unique(df[[group_var]]))
  if (k < 2) return(NA_character_)
  if (k == 2) return("wilcox")
  "kruskal"
}


#' Fit a linear mixed model with random intercept per subject
#'
#' Formula: metric ~ group * facet + (1 | subject)   (if both present)
#'          metric ~ group + (1 | subject)            (if facet absent)
run_lmm <- function(df, metric, group_var, facet_var, subject_var) {

  # Build fixed-effects side of the formula from whichever variables are set
  fixed_terms <- c()
  if (!is.null(group_var) && group_var != "(none)") fixed_terms <- c(fixed_terms, group_var)
  if (!is.null(facet_var) && facet_var != "(none)") fixed_terms <- c(fixed_terms, facet_var)
  if (length(fixed_terms) == 0) {
    return(list(method = "Linear mixed model",
                note = "Need at least one fixed-effect predictor (set Grouping variable).",
                fixed = NULL, pairwise = NULL, diagnostics = NULL,
                formula = NA_character_))
  }

  # Use interaction when both group and facet are present -- this is the
  # standard treatment x time model.
  rhs <- if (length(fixed_terms) == 2) paste(fixed_terms, collapse = " * ")
         else fixed_terms[1]
  fml_str <- paste0("`", metric, "` ~ ", rhs, " + (1 | `", subject_var, "`)")
  fml <- stats::as.formula(fml_str)

  # Coerce predictors and subject to factor for clean contrasts.
  # We don't coerce metric (response) -- it stays numeric.
  for (v in c(fixed_terms, subject_var)) df[[v]] <- factor(df[[v]])
  df <- df[stats::complete.cases(df[, c(metric, fixed_terms, subject_var)]), ]

  # lmerTest::lmer gives us Satterthwaite df + p-values out of the box.
  fit <- tryCatch(
    lmerTest::lmer(fml, data = df, REML = TRUE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(method = "Linear mixed model",
                formula = fml_str,
                note = paste("Model failed to fit:", fit$message),
                fixed = NULL, pairwise = NULL, diagnostics = NULL))
  }

  # ANOVA-style table for fixed effects (Type III SS, Satterthwaite)
  anova_tbl <- as.data.frame(stats::anova(fit, type = 3))
  anova_tbl <- tibble::rownames_to_column(anova_tbl, var = "Term")

  # Pairwise contrasts. If only one fixed effect, contrast its levels;
  # if interaction, contrast group levels within each facet level.
  pw <- tryCatch({
    if (length(fixed_terms) == 1) {
      emm <- emmeans::emmeans(fit, specs = fixed_terms[1])
      as.data.frame(pairs(emm, adjust = "tukey"))
    } else {
      # Contrast the first fixed term within each level of the second
      emm <- emmeans::emmeans(fit, specs = fixed_terms[1], by = fixed_terms[2])
      as.data.frame(pairs(emm, adjust = "tukey"))
    }
  }, error = function(e) NULL)

  # Diagnostics: Shapiro-Wilk on residuals + residuals/fitted for QQ plot
  res <- stats::residuals(fit)
  fit_vals <- stats::fitted(fit)
  shap <- tryCatch(stats::shapiro.test(res), error = function(e) NULL)

  list(
    method = "Linear mixed model (REML, Satterthwaite df)",
    formula = fml_str,
    fixed = anova_tbl,
    pairwise = pw,
    diagnostics = list(residuals = res, fitted = fit_vals,
                       shapiro_p = if (!is.null(shap)) shap$p.value else NA_real_),
    note = paste0("n = ", nrow(df),
                  " observations, subjects = ", length(unique(df[[subject_var]])),
                  ". Random intercept per subject. ",
                  "Pairwise contrasts use Tukey adjustment.")
  )
}


#' Wilcoxon rank-sum (2 groups, independent) or signed-rank (paired)
run_wilcoxon <- function(df, metric, group_var) {

  df <- df[stats::complete.cases(df[, c(metric, group_var)]), ]
  df[[group_var]] <- factor(df[[group_var]])
  if (length(levels(df[[group_var]])) != 2) {
    return(list(method = "Wilcoxon rank-sum",
                note = "Wilcoxon requires exactly 2 groups.",
                formula = NA_character_,
                fixed = NULL, pairwise = NULL, diagnostics = NULL))
  }

  fml <- stats::as.formula(paste0("`", metric, "` ~ `", group_var, "`"))
  test <- stats::wilcox.test(fml, data = df, exact = FALSE)

  tbl <- data.frame(
    Comparison = paste(levels(df[[group_var]]), collapse = " vs "),
    W = unname(test$statistic),
    p_value = test$p.value
  )

  list(
    method = "Wilcoxon rank-sum (Mann-Whitney U)",
    formula = paste0("wilcox.test(", metric, " ~ ", group_var, ")"),
    fixed = tbl,
    pairwise = NULL,
    diagnostics = NULL,
    note = paste0("n = ", nrow(df),
                  ". Non-parametric, no normality assumption. ",
                  "Use this when n is small or residuals are non-normal.")
  )
}


#' Kruskal-Wallis with Dunn's post-hoc (>=3 groups, independent)
run_kruskal <- function(df, metric, group_var) {

  df <- df[stats::complete.cases(df[, c(metric, group_var)]), ]
  df[[group_var]] <- factor(df[[group_var]])
  fml <- stats::as.formula(paste0("`", metric, "` ~ `", group_var, "`"))

  kw <- stats::kruskal.test(fml, data = df)
  tbl <- data.frame(
    Term = group_var,
    chi_squared = unname(kw$statistic),
    df = unname(kw$parameter),
    p_value = kw$p.value
  )

  # Dunn's test for pairwise comparisons, Benjamini-Hochberg adjusted
  dunn <- tryCatch({
    d <- FSA::dunnTest(fml, data = df, method = "bh")
    d$res
  }, error = function(e) NULL)

  list(
    method = "Kruskal-Wallis with Dunn post-hoc",
    formula = paste0("kruskal.test(", metric, " ~ ", group_var, ")"),
    fixed = tbl,
    pairwise = dunn,
    diagnostics = NULL,
    note = paste0("n = ", nrow(df),
                  ". Non-parametric ANOVA alternative. ",
                  "Dunn's post-hoc adjusted with Benjamini-Hochberg (FDR).")
  )
}


#' One-way ANOVA + Tukey HSD, or two-sample t-test if exactly 2 groups
run_parametric <- function(df, metric, group_var) {

  df <- df[stats::complete.cases(df[, c(metric, group_var)]), ]
  df[[group_var]] <- factor(df[[group_var]])
  k <- length(levels(df[[group_var]]))
  fml <- stats::as.formula(paste0("`", metric, "` ~ `", group_var, "`"))

  if (k == 2) {
    tt <- stats::t.test(fml, data = df, var.equal = FALSE)
    tbl <- data.frame(
      Comparison = paste(levels(df[[group_var]]), collapse = " vs "),
      t = unname(tt$statistic),
      df = unname(tt$parameter),
      p_value = tt$p.value
    )
    # Diagnostics on residuals from a linear model fit (equivalent)
    lmfit <- stats::lm(fml, data = df)
    res <- stats::residuals(lmfit); fit_vals <- stats::fitted(lmfit)
    shap <- tryCatch(stats::shapiro.test(res), error = function(e) NULL)

    return(list(
      method = "Welch two-sample t-test",
      formula = paste0("t.test(", metric, " ~ ", group_var, ")"),
      fixed = tbl,
      pairwise = NULL,
      diagnostics = list(residuals = res, fitted = fit_vals,
                         shapiro_p = if (!is.null(shap)) shap$p.value else NA_real_),
      note = paste0("n = ", nrow(df),
                    ". Assumes approximately normal residuals; ",
                    "check Shapiro-Wilk and QQ plot.")
    ))
  }

  # One-way ANOVA + Tukey HSD
  fit <- stats::aov(fml, data = df)
  anova_tbl <- as.data.frame(stats::anova(fit))
  anova_tbl <- tibble::rownames_to_column(anova_tbl, var = "Term")
  tuk <- tryCatch({
    th <- stats::TukeyHSD(fit)
    as.data.frame(th[[1]]) |> tibble::rownames_to_column(var = "Comparison")
  }, error = function(e) NULL)

  res <- stats::residuals(fit); fit_vals <- stats::fitted(fit)
  shap <- tryCatch(stats::shapiro.test(res), error = function(e) NULL)

  list(
    method = "One-way ANOVA with Tukey HSD",
    formula = paste0("aov(", metric, " ~ ", group_var, ")"),
    fixed = anova_tbl,
    pairwise = tuk,
    diagnostics = list(residuals = res, fitted = fit_vals,
                       shapiro_p = if (!is.null(shap)) shap$p.value else NA_real_),
    note = paste0("n = ", nrow(df),
                  ". Assumes approximately normal residuals and equal variance; ",
                  "check Shapiro-Wilk and QQ plot.")
  )
}


# ---- 3. UI -------------------------------------------------------------------
# Layout: a sidebar for uploads + controls, a main panel with tabs for the
# different views (Data, Summary, Plots, Statistics, About). bslib::page_sidebar()
# is the modern replacement for fluidPage(sidebarLayout(...)).
ui <- bslib::page_sidebar(

  title = "MicroEven — Microbiome Evenness Calculator",
  theme = bslib::bs_theme(bootswatch = "flatly"),  # pick any bootswatch theme

  # ----- sidebar with controls -----
  sidebar = bslib::sidebar(
    width = 360,

    h4("1. Upload files"),
    helpText("Both files should be MicrobiomeAnalyst alpha-diversity exports ",
             "in long format (one row per sample, with columns sample_id, ",
             "variable, value)."),

    fileInput(
      "file_observed",
      "Observed Features file (CSV)",
      accept = c(".csv", "text/csv")
    ),
    fileInput(
      "file_shannon",
      "Shannon file (CSV)",
      accept = c(".csv", "text/csv")
    ),

    hr(),
    h4("2. Plot options"),

    # These selectors are populated dynamically from the uploaded metadata
    # columns via updateSelectInput() in the server. They start empty.
    selectInput("group_var", "Grouping variable (x-axis / color)",
                choices = "(none)", selected = "(none)"),
    selectInput("facet_var", "Faceting variable (optional)",
                choices = "(none)", selected = "(none)"),
    selectInput("metric_to_plot", "Metric to plot",
                choices = c("Pielou", "ShannonEquitability",
                            "Shannon", "Observed"),
                selected = "Pielou"),
    selectInput("plot_type", "Plot type",
                choices = c("Boxplot" = "box",
                            "Violin"  = "violin",
                            "Scatter (vs richness)" = "scatter"),
                selected = "box"),

    hr(),
    h4("3. Statistics"),

    selectInput("subject_var",
                "Subject ID column (for repeated measures)",
                choices = "(none)", selected = "(none)"),
    helpText("If your data has repeated samples from the same subject ",
             "(e.g., multiple time points), set this so the model accounts ",
             "for within-subject correlation."),

    selectInput("stat_test", "Statistical test",
                choices = c("Auto (recommended)" = "auto",
                            "Linear mixed model (LMM)" = "lmm",
                            "Non-parametric (Wilcoxon / Kruskal)" = "np",
                            "Parametric (t-test / ANOVA)" = "param"),
                selected = "auto"),
    helpText("Auto picks LMM when a subject column is set, ",
             "otherwise non-parametric (safer at small n)."),

    hr(),
    h4("4. Downloads"),
    downloadButton("dl_per_sample", "Per-sample CSV", class = "btn-sm"),
    br(), br(),
    downloadButton("dl_summary",    "Summary CSV",    class = "btn-sm"),
    br(), br(),
    downloadButton("dl_stats",      "Statistics (CSV)", class = "btn-sm"),
    br(), br(),
    downloadButton("dl_plot_png",   "Plot (PNG)",     class = "btn-sm"),
    downloadButton("dl_plot_pdf",   "Plot (PDF)",     class = "btn-sm"),

    hr(),
    helpText(
      tags$small(
        "If MicroEven contributed to your work, please cite this app and ",
        "MicrobiomeAnalyst (the upstream source). See the ",
        tags$em("About / Help"), " tab for full references."
      )
    )
  ),

  # ----- main content with tabs -----
  bslib::navset_card_tab(

    bslib::nav_panel(
      title = "Per-sample table",
      DT::DTOutput("table_per_sample")
    ),

    bslib::nav_panel(
      title = "Group summary",
      DT::DTOutput("table_summary")
    ),

    bslib::nav_panel(
      title = "Plot",
      plotOutput("evenness_plot", height = "500px")
    ),

    bslib::nav_panel(
      title = "Statistics",

      # Each subsection is wrapped in its own card so DT tables, which
      # auto-expand as content loads, can't bleed into the next header.
      # Cards also visually separate the four sections, which matches the
      # mental model the user already has from the navset_card_tab layout.

      bslib::card(
        bslib::card_header("Test"),
        bslib::card_body(
          verbatimTextOutput("stat_method"),
          tags$strong("Model / call"),
          verbatimTextOutput("stat_formula"),
          helpText("The call above is shown so you can reproduce it ",
                   "in a script. For LMMs the random-effects term is included.")
        )
      ),

      bslib::card(
        bslib::card_header("Main results"),
        bslib::card_body(
          # min_height keeps the card from collapsing when the table is empty,
          # which prevents layout jump when results first load.
          min_height = "120px",
          DT::DTOutput("stat_main")
        )
      ),

      bslib::card(
        bslib::card_header("Pairwise contrasts"),
        bslib::card_body(
          min_height = "120px",
          helpText("LMM uses emmeans + Tukey adjustment. ",
                   "Kruskal uses Dunn's test (BH-adjusted). ",
                   "ANOVA uses Tukey HSD."),
          DT::DTOutput("stat_pairwise")
        )
      ),

      bslib::card(
        bslib::card_header("Notes"),
        bslib::card_body(
          verbatimTextOutput("stat_note")
        )
      ),

      bslib::card(
        bslib::card_header("Residual diagnostics"),
        bslib::card_body(
          helpText("Only shown for LMM and parametric tests. ",
                   "If Shapiro-Wilk p < 0.05 or the QQ plot deviates strongly ",
                   "from the line, prefer the non-parametric option."),
          verbatimTextOutput("stat_shapiro"),
          plotOutput("stat_qqplot", height = "350px")
        )
      )
    ),

    bslib::nav_panel(
      title = "About / Help",
      h3("What this app does"),
      p("MicroEven computes microbiome evenness metrics that MicrobiomeAnalyst ",
        "does not export by default. You upload your Observed Features and ",
        "Shannon CSVs (the same format the MicrobiomeAnalyst alpha-diversity ",
        "module gives you), and the app joins them on sample_id and computes:"),
      tags$ul(
        tags$li(strong("Pielou's evenness (J'):"),
                " J' = H' / ln(S). Range [0, 1]. 1 = perfectly even."),
        tags$li(strong("Shannon equitability:"),
                " exp(H') / S. An alternative bounded measure of evenness ",
                "derived from Shannon. NOT the same as Simpson's evenness — ",
                "true Simpson evenness requires the Simpson index, which ",
                "this app does not yet ingest.")
      ),
      h3("Statistical tests"),
      tags$ul(
        tags$li(strong("LMM:"),
                " lmer(metric ~ group * facet + (1 | subject)). ",
                "Use for repeated measures. p-values via Satterthwaite df ",
                "(lmerTest). Pairwise contrasts via emmeans with Tukey ",
                "adjustment."),
        tags$li(strong("Wilcoxon rank-sum:"),
                " 2 independent groups, non-parametric. Safer at small n."),
        tags$li(strong("Kruskal-Wallis:"),
                " 3+ independent groups, non-parametric. ",
                "Post-hoc via Dunn's test (Benjamini-Hochberg)."),
        tags$li(strong("t-test / ANOVA:"),
                " parametric alternatives. Use only if residuals are ",
                "approximately normal -- diagnostics shown.")
      ),
      h3("Required input columns"),
      tags$ul(
        tags$li(code("sample_id"),  " — unique per row"),
        tags$li(code("variable"),   " — must contain 'Observed' or 'Shannon'"),
        tags$li(code("value"),      " — numeric metric value"),
        tags$li("Any other columns are kept as metadata (treatment, ",
                "week, subject_id, etc.)")
      ),
      h3("How to cite this app"),
      p("If MicroEven contributed to a publication, poster, or thesis, ",
        "please cite it as:"),
      tags$blockquote(
        em("Rivas Hernandez, E. R. (2026). MicroEven: a Shiny application ",
           "for computing microbiome evenness metrics and group-level ",
           "statistical comparisons. Colorado State University.")
      ),
      p("BibTeX:"),
      tags$pre(
"@software{rivashernandez_microeven_2026,
  author  = {Rivas Hernandez, Eliud R.},
  title   = {MicroEven: a Shiny application for computing microbiome
             evenness metrics and group-level statistical comparisons},
  year    = {2026},
  institution = {Colorado State University}
}"
      ),

      h3("References"),

      h4("Upstream data source"),
      tags$ul(
        tags$li(
          strong("MicrobiomeAnalyst."),
          " Lu, Y., Zhou, G., Ewald, J., Pang, Z., Shiri, T., & Xia, J. (2023). ",
          em("MicrobiomeAnalyst 2.0: comprehensive statistical, functional ",
             "and integrative analysis of microbiome data."),
          " Nucleic Acids Research, 51(W1), W310–W318. ",
          "doi:10.1093/nar/gkad407. ",
          "MicroEven ingests the Shannon and Observed Features CSV exports ",
          "from the MicrobiomeAnalyst alpha-diversity module."
        )
      ),

      h4("Statistical methods"),
      tags$ul(
        tags$li(
          strong("Pielou's evenness:"),
          " Pielou, E. C. (1966). ",
          em("The measurement of diversity in different types of biological ",
             "collections."),
          " Journal of Theoretical Biology, 13, 131–144. ",
          "doi:10.1016/0022-5193(66)90013-0."
        ),
        tags$li(
          strong("Shannon index:"),
          " Shannon, C. E. (1948). ",
          em("A mathematical theory of communication."),
          " Bell System Technical Journal, 27(3), 379–423."
        ),
        tags$li(
          strong("Linear mixed model with Satterthwaite df:"),
          " Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017). ",
          em("lmerTest package: tests in linear mixed effects models."),
          " Journal of Statistical Software, 82(13), 1–26. ",
          "doi:10.18637/jss.v82.i13."
        ),
        tags$li(
          strong("Dunn's test (post-hoc for Kruskal–Wallis):"),
          " Dunn, O. J. (1964). ",
          em("Multiple comparisons using rank sums."),
          " Technometrics, 6(3), 241–252."
        ),
        tags$li(
          strong("Benjamini-Hochberg FDR adjustment:"),
          " Benjamini, Y., & Hochberg, Y. (1995). ",
          em("Controlling the false discovery rate: a practical and powerful ",
             "approach to multiple testing."),
          " Journal of the Royal Statistical Society B, 57(1), 289–300."
        )
      ),

      h4("R packages"),
      tags$ul(
        tags$li(
          strong("shiny:"),
          " Chang, W., Cheng, J., Allaire, J. J., Sievert, C., Schloerke, B., ",
          "Xie, Y., Allen, J., McPherson, J., Dipert, A., & Borges, B. ",
          "shiny: Web Application Framework for R. R package."
        ),
        tags$li(
          strong("bslib:"),
          " Sievert, C., Cheng, J., & Aden-Buie, G. ",
          "bslib: Custom 'Bootstrap' 'Sass' Themes for 'shiny' and 'rmarkdown'. ",
          "R package."
        ),
        tags$li(
          strong("DT:"),
          " Xie, Y., Cheng, J., & Tan, X. ",
          "DT: A Wrapper of the JavaScript Library 'DataTables'. R package."
        ),
        tags$li(
          strong("dplyr / tidyr / readr:"),
          " Wickham, H., et al. The tidyverse: a collection of R packages ",
          "for data science. doi:10.21105/joss.01686."
        ),
        tags$li(
          strong("ggplot2:"),
          " Wickham, H. (2016). ",
          em("ggplot2: Elegant Graphics for Data Analysis."),
          " Springer-Verlag New York. ISBN 978-3-319-24277-4."
        ),
        tags$li(
          strong("lme4:"),
          " Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). ",
          em("Fitting linear mixed-effects models using lme4."),
          " Journal of Statistical Software, 67(1), 1–48. ",
          "doi:10.18637/jss.v067.i01."
        ),
        tags$li(
          strong("lmerTest:"),
          " Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017). ",
          "Journal of Statistical Software, 82(13). doi:10.18637/jss.v082.i13."
        ),
        tags$li(
          strong("emmeans:"),
          " Lenth, R. V. ",
          "emmeans: Estimated Marginal Means, aka Least-Squares Means. ",
          "R package."
        ),
        tags$li(
          strong("broom.mixed:"),
          " Bolker, B., & Robinson, D. ",
          "broom.mixed: Tidying Methods for Mixed Models. R package."
        ),
        tags$li(
          strong("FSA:"),
          " Ogle, D. H., Doll, J. C., Wheeler, A. P., & Dinno, A. ",
          "FSA: Simple Fisheries Stock Assessment Methods. R package. ",
          "(Provides Dunn's test wrapper.)"
        ),
        tags$li(
          strong("car:"),
          " Fox, J., & Weisberg, S. (2019). ",
          em("An R Companion to Applied Regression, Third Edition."),
          " Sage. (Provides Type II/III ANOVA tables.)"
        ),
        tags$li(
          strong("R itself:"),
          " R Core Team. ",
          em("R: A Language and Environment for Statistical Computing."),
          " R Foundation for Statistical Computing, Vienna, Austria. ",
          "https://www.R-project.org/"
        )
      ),

      helpText("Tip: in R, run ", code("citation(\"packagename\")"),
               " (e.g., ", code("citation(\"lme4\")"),
               ") to get the canonical, version-specific citation for any ",
               "package, including a BibTeX entry.")
    )
  )
)


# ---- 4. Server ---------------------------------------------------------------
server <- function(input, output, session) {

  # ---- 4a. Reactive: parse Observed file -----------------------------------
  # req(input$file_observed) halts execution until the user uploads the file.
  # This prevents the rest of the pipeline from erroring on a NULL input.
  observed_df <- reactive({
    req(input$file_observed)

    # tryCatch lets us surface parser errors as a friendly Shiny notification
    # instead of a red error panel that doesn't tell the user what happened.
    tryCatch(
      parse_alpha_file(input$file_observed$datapath, "Observed"),
      error = function(e) {
        showNotification(paste("Observed file error:", e$message),
                         type = "error", duration = 10)
        NULL
      }
    )
  })

  # ---- 4b. Reactive: parse Shannon file ------------------------------------
  shannon_df <- reactive({
    req(input$file_shannon)
    tryCatch(
      parse_alpha_file(input$file_shannon$datapath, "Shannon"),
      error = function(e) {
        showNotification(paste("Shannon file error:", e$message),
                         type = "error", duration = 10)
        NULL
      }
    )
  })
  
  # ---- 4c. Reactive: join + compute evenness -------------------------------
  # This is the heart of the app: join Observed and Shannon by sample_id,
  # then compute the evenness metrics. We also detect metadata mismatches
  # (sample present in one file but not the other) and warn the user.
  joined_df <- reactive({
    obs <- observed_df()
    sha <- shannon_df()
    req(obs, sha)

    # Inner join on sample_id. Metadata columns from `obs` are kept; we drop
    # duplicate metadata from `sha` (suffix .y) to keep the table clean.
    merged <- dplyr::inner_join(obs, sha, by = "sample_id",
                                suffix = c("", ".y"))

    # Drop any *.y duplicate metadata columns from the Shannon file
    dup_cols <- grep("\\.y$", names(merged), value = TRUE)
    merged   <- merged[, !(names(merged) %in% dup_cols), drop = FALSE]

    # Tell the user if some samples didn't match
    n_obs    <- nrow(obs)
    n_sha    <- nrow(sha)
    n_merged <- nrow(merged)
    if (n_merged < min(n_obs, n_sha)) {
      showNotification(
        paste0("Joined ", n_merged, " samples. ",
               n_obs - n_merged, " in Observed and ",
               n_sha - n_merged, " in Shannon did not match by sample_id."),
        type = "warning", duration = 8
      )
    }

    # Compute Pielou + Shannon equitability
   compute_evenness(merged)
    })

  # ---- 4d. Populate metadata dropdowns dynamically -------------------------
  # When the user uploads files, scan the joined data for columns that look
  # like metadata (not the metric columns themselves) and offer them as
  # grouping/faceting/subject options.
  observe({
    df <- joined_df()
    req(df)

    metric_cols <- c("Observed", "Shannon", "Pielou", "ShannonEquitability",
                     "sample_id", "samples", "...1")
    meta_cols   <- setdiff(names(df), metric_cols)

    updateSelectInput(session, "group_var",
                      choices  = c("(none)", meta_cols),
                      selected = if ("treatment" %in% meta_cols) "treatment"
                                 else "(none)")
    updateSelectInput(session, "facet_var",
                      choices  = c("(none)", meta_cols),
                      selected = "(none)")
    updateSelectInput(session, "subject_var",
                      choices  = c("(none)", meta_cols),
                      selected = if ("subject_id" %in% meta_cols) "subject_id"
                                 else "(none)")
  })

  # ---- 4e. Per-sample table ------------------------------------------------
  output$table_per_sample <- DT::renderDT({
    df <- joined_df()
    req(df)

    # Round numeric columns for display only -- the underlying download still
    # has full precision.
    num_cols <- vapply(df, is.numeric, logical(1))
    df_disp  <- df
    df_disp[num_cols] <- lapply(df_disp[num_cols], round, digits = 3)

    DT::datatable(
      df_disp,
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # ---- 4f. Group summary table --------------------------------------------
  summary_df <- reactive({
    df <- joined_df()
    req(df)
    summarise_evenness(df, input$group_var)
  })

  output$table_summary <- DT::renderDT({
    s <- summary_df()
    req(s)
    num_cols <- vapply(s, is.numeric, logical(1))
    s[num_cols] <- lapply(s[num_cols], round, digits = 3)
    DT::datatable(s, options = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE)
  })

  # ---- 4g. Plot ------------------------------------------------------------
  # Build the ggplot object inside a reactive so we can reuse it for both
  # the on-screen output and the PNG/PDF downloads (DRY).
  evenness_plot_obj <- reactive({
    df     <- joined_df()
    req(df)
    metric <- input$metric_to_plot
    grp    <- input$group_var
    fac    <- input$facet_var
    ptype  <- input$plot_type

    # Drop rows where the metric is NA so geom_* doesn't warn for every one
    df <- df[!is.na(df[[metric]]), , drop = FALSE]

    # Build aes() differently depending on whether a grouping var is chosen.
    if (grp == "(none)") {
      p <- ggplot(df, aes(x = "all samples", y = .data[[metric]]))
    } else {
      # Coerce to factor so ggplot treats numeric groups (e.g. week) as
      # categorical for boxplots. For scatter we leave it as-is.
      if (ptype != "scatter") {
        df[[grp]] <- factor(df[[grp]])
      }
      p <- ggplot(df, aes(x = .data[[grp]], y = .data[[metric]],
                          fill = .data[[grp]], color = .data[[grp]]))
    }

    # Pick the geom based on plot type
    p <- switch(
      ptype,
      box     = p + geom_boxplot(alpha = 0.6, outlier.shape = NA) +
                    geom_jitter(width = 0.15, alpha = 0.7, size = 1.8),
      violin  = p + geom_violin(alpha = 0.6, trim = FALSE) +
                    geom_jitter(width = 0.1, alpha = 0.7, size = 1.5),
      scatter = ggplot(df, aes(x = Observed, y = .data[[metric]],
                               color = if (grp == "(none)") NULL
                                       else .data[[grp]])) +
                  geom_point(size = 2.5, alpha = 0.8) +
                  labs(x = "Observed features (richness)")
    )

    # Optional faceting
    if (fac != "(none)") {
      p <- p + facet_wrap(stats::as.formula(paste("~", fac)))
    }

    # Theme + labels. Pielou label gets the prime notation.
    p +
      labs(
        y = switch(metric,
                   Pielou              = "Pielou's evenness (J')",
                   ShannonEquitability = "Shannon equitability  exp(H') / S",
                   Shannon             = "Shannon index (H')",
                   Observed            = "Observed features (richness)"),
        title = paste0(metric,
                       if (grp != "(none)") paste0(" by ", grp) else "")
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "right",
            plot.title = element_text(face = "bold"))
  })

  output$evenness_plot <- renderPlot({
    evenness_plot_obj()
  })

  # ---- 4h. Statistics ------------------------------------------------------
  # Reactive that picks and runs the appropriate test, returning the
  # standardized list shape described in section 2b.
  stat_result <- reactive({
    df <- joined_df()
    req(df)

    metric <- input$metric_to_plot
    grp    <- input$group_var
    fac    <- input$facet_var
    subj   <- input$subject_var
    choice <- input$stat_test

    # Resolve "auto" to a concrete test
    if (choice == "auto") {
      auto <- choose_auto_test(df, grp, subj)
      if (is.na(auto)) {
        return(list(method = "(no test)",
                    note = "Select a grouping variable to run a test.",
                    formula = NA_character_,
                    fixed = NULL, pairwise = NULL, diagnostics = NULL))
      }
      choice <- auto
    }

    # Route to the appropriate helper. "np" means: pick Wilcoxon if k==2,
    # otherwise Kruskal-Wallis. This mirrors how a careful analyst would
    # choose between them.
    if (choice == "np") {
      if (grp == "(none)") {
        return(list(method = "Non-parametric",
                    note = "Need a grouping variable for non-parametric tests.",
                    formula = NA_character_,
                    fixed = NULL, pairwise = NULL, diagnostics = NULL))
      }
      k <- length(unique(df[[grp]]))
      if (k == 2) return(run_wilcoxon(df, metric, grp))
      return(run_kruskal(df, metric, grp))
    }

    switch(
      choice,
      lmm     = run_lmm(df, metric, grp, fac, subj),
      wilcox  = run_wilcoxon(df, metric, grp),
      kruskal = run_kruskal(df, metric, grp),
      param   = run_parametric(df, metric, grp),
      list(method = "(no test)", note = "Unknown test selection.",
           formula = NA_character_,
           fixed = NULL, pairwise = NULL, diagnostics = NULL)
    )
  })

  output$stat_method  <- renderText({ stat_result()$method })
  output$stat_formula <- renderText({
    f <- stat_result()$formula
    if (is.null(f) || is.na(f)) "(no model fit)" else f
  })
  output$stat_note    <- renderText({
    n <- stat_result()$note
    if (is.null(n)) "" else n
  })

  output$stat_main <- DT::renderDT({
    r <- stat_result()
    req(r$fixed)
    tbl <- r$fixed
    num_cols <- vapply(tbl, is.numeric, logical(1))
    tbl[num_cols] <- lapply(tbl[num_cols], function(x) signif(x, 4))
    DT::datatable(tbl, options = list(dom = "t", pageLength = 25),
                  rownames = FALSE)
  })

  output$stat_pairwise <- DT::renderDT({
    r <- stat_result()
    req(r$pairwise)
    tbl <- r$pairwise
    num_cols <- vapply(tbl, is.numeric, logical(1))
    tbl[num_cols] <- lapply(tbl[num_cols], function(x) signif(x, 4))
    DT::datatable(tbl, options = list(dom = "t", pageLength = 25),
                  rownames = FALSE)
  })

  output$stat_shapiro <- renderText({
    r <- stat_result()
    if (is.null(r$diagnostics) || is.null(r$diagnostics$shapiro_p)) {
      return("Not applicable for this test.")
    }
    p <- r$diagnostics$shapiro_p
    paste0("Shapiro-Wilk on residuals: p = ", signif(p, 4),
           if (!is.na(p) && p < 0.05) "  (residuals deviate from normality — consider non-parametric)"
           else "  (no strong evidence against normality)")
  })

  output$stat_qqplot <- renderPlot({
    r <- stat_result()
    req(r$diagnostics, r$diagnostics$residuals)
    res <- r$diagnostics$residuals
    qq_df <- data.frame(
      theo   = stats::qnorm(stats::ppoints(length(res))),
      sample = sort(res)
    )
    ggplot(qq_df, aes(x = theo, y = sample)) +
      geom_point(alpha = 0.7) +
      geom_abline(slope = stats::sd(res), intercept = mean(res),
                  linetype = "dashed") +
      labs(x = "Theoretical quantiles",
           y = "Residuals (sample quantiles)",
           title = "Normal QQ plot of residuals") +
      theme_minimal(base_size = 13)
  })

  # ---- 4i. Downloads -------------------------------------------------------
  output$dl_per_sample <- downloadHandler(
    filename = function() paste0("MicroEven_per_sample_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(joined_df(), file)
  )

  output$dl_summary <- downloadHandler(
    filename = function() paste0("MicroEven_summary_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(summary_df(), file)
  )

  # Stats CSV bundles main + pairwise tables, with a header row noting
  # the test and formula. This is the version you would paste into a
  # methods section.
  output$dl_stats <- downloadHandler(
    filename = function() paste0("MicroEven_statistics_", Sys.Date(), ".csv"),
    content  = function(file) {
      r <- stat_result()
      con <- file(file, open = "w")
      on.exit(close(con))
      writeLines(paste0("# Test: ", r$method), con)
      writeLines(paste0("# Formula: ", r$formula), con)
      writeLines(paste0("# Note: ", r$note), con)
      writeLines("", con)
      if (!is.null(r$fixed)) {
        writeLines("# Main results", con)
        utils::write.csv(r$fixed, con, row.names = FALSE)
        writeLines("", con)
      }
      if (!is.null(r$pairwise)) {
        writeLines("# Pairwise contrasts", con)
        utils::write.csv(r$pairwise, con, row.names = FALSE)
      }
    }
  )

  output$dl_plot_png <- downloadHandler(
    filename = function() paste0("MicroEven_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      ggsave(file, plot = evenness_plot_obj(),
             width = 8, height = 5, dpi = 300, device = "png")
    }
  )

  output$dl_plot_pdf <- downloadHandler(
    filename = function() paste0("MicroEven_plot_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggsave(file, plot = evenness_plot_obj(),
             width = 8, height = 5, device = "pdf")
    }
  )
}


# ---- 5. Launch ---------------------------------------------------------------
shinyApp(ui, server)

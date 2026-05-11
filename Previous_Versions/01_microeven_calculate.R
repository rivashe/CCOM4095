# =============================================================================
# MicroEven Step 1 — Calculate Pielou's Evenness
# =============================================================================
# Author : Eliud R. Rivas Hernandez
# Advisor: Dr. Tiffany Weir
# Course : CCOM 4095
#
# PURPOSE
# -------
# STEP 1 of 3. Loads the two MicrobiomeAnalyst alpha-diversity exports
# (Observed Features and Shannon), joins them on sample_id, and computes
# Pielou's evenness (J' = H' / ln S).
#
# Displays:
#   1. Per-sample table
#   2. Group summary table
#
#
#
# MATH NOTE
# ---------
# Pielou's J' = H' / ln(S)
#   H' = Shannon index      (Shannon CSV)
#   S  = Observed features  (Observed CSV)
# Range [0, 1]. J' = 1 means a perfectly even community.
#
# INPUT FORMAT
# ------------
# Required columns: sample_id | variable | value
# Extra columns are kept as metadata (treatment, week, subject_id, etc.)
#
# DEPENDENCIES
# ------------
# install.packages(c("shiny", "dplyr", "readr", "DT", "bslib"))
# =============================================================================

library(shiny)
library(dplyr)
library(readr)
library(DT)
library(bslib)


# ---- Helper functions -------------------------------------------------------

parse_alpha_file <- function(path, expected_metric) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  required <- c("sample_id", "variable", "value")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("File is missing required column(s): ",
         paste(missing, collapse = ", "),
         ". Expected a MicrobiomeAnalyst alpha-diversity export.")
  }

  metrics_in_file <- unique(df$variable)
  if (!any(tolower(metrics_in_file) == tolower(expected_metric))) {
    stop("Expected metric '", expected_metric, "' not found in file. ",
         "Found: ", paste(metrics_in_file, collapse = ", "))
  }

  df <- df[tolower(df$variable) == tolower(expected_metric), , drop = FALSE]
  df[[expected_metric]] <- df$value
  df$variable <- NULL
  df$value    <- NULL
  df <- df[, names(df) != "...1", drop = FALSE]
  df
}

compute_pielou <- function(df) {
  df$Pielou <- ifelse(
    df$Observed > 1,
    df$Shannon / log(df$Observed),
    NA_real_
  )
  df
}

summarise_evenness <- function(df, group_var = NULL) {
  metric_cols <- intersect(c("Observed", "Shannon", "Pielou"), names(df))

  grouped <- if (is.null(group_var) || group_var == "(none)") df
             else dplyr::group_by(df, .data[[group_var]])

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


# ---- UI ---------------------------------------------------------------------

ui <- bslib::page_sidebar(

  title = "MicroEven Step 1 — Calculate Pielou's Evenness",
  theme = bslib::bs_theme(bootswatch = "flatly"),

  sidebar = bslib::sidebar(
    width = 340,

    h4("1. Upload files"),
    helpText("Both files should be MicrobiomeAnalyst alpha-diversity exports ",
             "in long format (one row per sample, with columns sample_id, ",
             "variable, value)."),

    fileInput("file_observed", "Observed Features file (CSV)",
              accept = c(".csv", "text/csv")),
    fileInput("file_shannon",  "Shannon file (CSV)",
              accept = c(".csv", "text/csv")),

    hr(),
    h4("2. Summary options"),
    selectInput("group_var", "Group by (for summary table)",
                choices = "(none)", selected = "(none)"),

    hr(),
    h4("3. Downloads"),
    downloadButton("dl_per_sample", "Per-sample CSV", class = "btn-sm"),
    br(), br(),
    downloadButton("dl_summary",    "Summary CSV",    class = "btn-sm")
  ),

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
      title = "About / Help",
      h3("What this step does"),
      p("Step 1 loads your Observed Features and Shannon exports from ",
        "MicrobiomeAnalyst, joins them on sample_id, and computes ",
        "Pielou's evenness for each sample"),
      h3("Pielou's evenness"),
      tags$ul(
        tags$li(strong("J' = H' / ln(S)."),
                " Range [0, 1]. 1 = perfectly even community."),
        tags$li("Requires only H' and S, both of which come from your ",
                "MicrobiomeAnalyst exports. The calculation is exact.")
      ),
      h3("Required input columns"),
      tags$ul(
        tags$li(code("sample_id"),  " — unique per row"),
        tags$li(code("variable"),   " — must contain 'Observed' or 'Shannon'"),
        tags$li(code("value"),      " — numeric metric value"),
        tags$li("Any other columns are kept as metadata.")
      ),
      h3("Citation"),
      p("Pielou, E.C. (1966). The measurement of diversity in different ",
        "types of biological collections. Journal of Theoretical Biology, ",
        "13, 131–144.")
    )
  )
)


# ---- Server -----------------------------------------------------------------

server <- function(input, output, session) {

  observed_df <- reactive({
    req(input$file_observed)
    tryCatch(
      parse_alpha_file(input$file_observed$datapath, "Observed"),
      error = function(e) {
        showNotification(paste("Observed file error:", e$message),
                         type = "error", duration = 10)
        NULL
      }
    )
  })

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

  joined_df <- reactive({
    obs <- observed_df(); sha <- shannon_df()
    req(obs, sha)

    merged <- dplyr::inner_join(obs, sha, by = "sample_id", suffix = c("", ".y"))
    dup_cols <- grep("\\.y$", names(merged), value = TRUE)
    merged   <- merged[, !(names(merged) %in% dup_cols), drop = FALSE]

    n_obs <- nrow(obs); n_sha <- nrow(sha); n_merged <- nrow(merged)
    if (n_merged < min(n_obs, n_sha)) {
      showNotification(
        paste0("Joined ", n_merged, " samples. ",
               n_obs - n_merged, " in Observed and ",
               n_sha - n_merged, " in Shannon did not match by sample_id."),
        type = "warning", duration = 8
      )
    }

    compute_pielou(merged)
  })

  observe({
    df <- joined_df(); req(df)
    metric_cols <- c("Observed", "Shannon", "Pielou",
                     "sample_id", "samples", "...1")
    meta_cols   <- setdiff(names(df), metric_cols)
    updateSelectInput(session, "group_var",
                      choices  = c("(none)", meta_cols),
                      selected = if ("treatment" %in% meta_cols) "treatment"
                                 else "(none)")
  })

  output$table_per_sample <- DT::renderDT({
    df <- joined_df(); req(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = 3)
    DT::datatable(df, options = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE)
  })

  summary_df <- reactive({
    df <- joined_df(); req(df)
    summarise_evenness(df, input$group_var)
  })

  output$table_summary <- DT::renderDT({
    s <- summary_df(); req(s)
    num_cols <- vapply(s, is.numeric, logical(1))
    s[num_cols] <- lapply(s[num_cols], round, digits = 3)
    DT::datatable(s, options = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE)
  })

  output$dl_per_sample <- downloadHandler(
    filename = function() paste0("step1_per_sample_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(joined_df(), file)
  )

  output$dl_summary <- downloadHandler(
    filename = function() paste0("step1_summary_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(summary_df(), file)
  )
}

shinyApp(ui, server)

# =============================================================================
# MicroEven Step 2 — Tables + Plot
# =============================================================================
# Author : Eliud R. Rivas Hernandez
# Advisor: Dr. Tiffany Weir
# Course : CCOM 4095
#
# PURPOSE
# -------
# STEP 2 of 3. Everything in Step 1, plus an interactive plot panel.
# Adds boxplot / violin / scatter visualisations of Pielou's evenness
# (and other metrics) across groups and facets, with PNG/PDF export.
#
# Displays:
#   1. Per-sample table
#   2. Group summary table
#   3. Plot (boxplot / violin / scatter)
#
# Statistical tests are in Step 3 (the full MicroEven app).
#
# MATH NOTE
# ---------
# Pielou's J' = H' / ln(S).  Range [0, 1].
# Shannon equitability = exp(H') / S.
#   This is NOT the same as Simpson evenness -- labelled clearly in the UI.
#
# INPUT FORMAT
# ------------
# Required columns: sample_id | variable | value
# Extra columns are kept as metadata (treatment, week, subject_id, etc.)
#
# DEPENDENCIES
# ------------
# install.packages(c("shiny", "dplyr", "readr", "ggplot2", "DT", "bslib"))
# =============================================================================

library(shiny)
library(dplyr)
library(readr)
library(ggplot2)
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

compute_evenness <- function(df) {
  df$Pielou <- ifelse(
    df$Observed > 1,
    df$Shannon / log(df$Observed),
    NA_real_
  )
  df$ShannonEquitability <- ifelse(
    df$Observed > 0,
    exp(df$Shannon) / df$Observed,
    NA_real_
  )
  df
}

summarise_evenness <- function(df, group_var = NULL) {
  metric_cols <- c("Observed", "Shannon", "Pielou", "ShannonEquitability")
  metric_cols <- intersect(metric_cols, names(df))

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

  title = "MicroEven Step 2 — Tables + Plot",
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
    h4("2. Plot options"),

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
    h4("3. Downloads"),
    downloadButton("dl_per_sample", "Per-sample CSV", class = "btn-sm"),
    br(), br(),
    downloadButton("dl_summary",    "Summary CSV",    class = "btn-sm"),
    br(), br(),
    downloadButton("dl_plot_png",   "Plot (PNG)",     class = "btn-sm"),
    downloadButton("dl_plot_pdf",   "Plot (PDF)",     class = "btn-sm")
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
      title = "Plot",
      plotOutput("evenness_plot", height = "500px")
    ),

    bslib::nav_panel(
      title = "About / Help",
      h3("What this step does"),
      p("Step 2 builds on Step 1 by adding an interactive plot panel. ",
        "You can switch between boxplot, violin, and scatter geometries, ",
        "choose any metadata column for grouping or faceting, and export ",
        "publication-ready PNG and PDF figures. Statistical tests are in ",
        "Step 3."),
      h3("Metrics"),
      tags$ul(
        tags$li(strong("Pielou's evenness (J'):"),
                " J' = H' / ln(S). Range [0, 1]. 1 = perfectly even."),
        tags$li(strong("Shannon equitability:"),
                " exp(H') / S. An alternative bounded measure of evenness ",
                "derived from Shannon. NOT the same as Simpson's evenness — ",
                "true Simpson evenness requires the Simpson index, which ",
                "this app does not yet ingest.")
      ),
      h3("Required input columns"),
      tags$ul(
        tags$li(code("sample_id"),  " — unique per row"),
        tags$li(code("variable"),   " — must contain 'Observed' or 'Shannon'"),
        tags$li(code("value"),      " — numeric metric value"),
        tags$li("Any other columns are kept as metadata (treatment, ",
                "week, subject_id, etc.)")
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

    compute_evenness(merged)
  })

  observe({
    df <- joined_df(); req(df)
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

  evenness_plot_obj <- reactive({
    df     <- joined_df()
    req(df)
    metric <- input$metric_to_plot
    grp    <- input$group_var
    fac    <- input$facet_var
    ptype  <- input$plot_type

    df <- df[!is.na(df[[metric]]), , drop = FALSE]

    if (grp == "(none)") {
      p <- ggplot(df, aes(x = "all samples", y = .data[[metric]]))
    } else {
      if (ptype != "scatter") {
        df[[grp]] <- factor(df[[grp]])
      }
      p <- ggplot(df, aes(x = .data[[grp]], y = .data[[metric]],
                          fill = .data[[grp]], color = .data[[grp]]))
    }

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

    if (fac != "(none)") {
      p <- p + facet_wrap(stats::as.formula(paste("~", fac)))
    }

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

  output$dl_per_sample <- downloadHandler(
    filename = function() paste0("step2_per_sample_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(joined_df(), file)
  )

  output$dl_summary <- downloadHandler(
    filename = function() paste0("step2_summary_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(summary_df(), file)
  )

  output$dl_plot_png <- downloadHandler(
    filename = function() paste0("step2_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      ggsave(file, plot = evenness_plot_obj(),
             width = 8, height = 5, dpi = 300, device = "png")
    }
  )

  output$dl_plot_pdf <- downloadHandler(
    filename = function() paste0("step2_plot_", Sys.Date(), ".pdf"),
    content  = function(file) {
      ggsave(file, plot = evenness_plot_obj(),
             width = 8, height = 5, device = "pdf")
    }
  )
}

shinyApp(ui, server)

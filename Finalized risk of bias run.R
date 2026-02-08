# This rerun only has risk of bias elements in it so that we can regenerate those results in specific, as they are supposed to be the same as our first test run.

library(readxl)
library(dplyr)
library(metafor)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(stringr) # for filename cleaning

safe_filename <- function(x) {
  x <- gsub("[/\\:*?\"<>|]", "_", x)
  x <- gsub(" ", "_", x)
  return(x)
}

meta_data <- read_excel("Meta-analysis data (3).xlsx")

required_cols <- c(
  "Effect_Size", "Variance", "Risk_of_Bias",
  "Sample_Size", "Analysis_ID", "Outcome_Measure", "Study_Title"
)

missing_cols <- setdiff(required_cols, names(meta_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

meta_data <- meta_data %>%
  mutate(
    Effect_Size      = as.numeric(Effect_Size),
    Variance         = as.numeric(Variance),
    Risk_of_Bias     = as.numeric(Risk_of_Bias),
    Sample_Size      = as.numeric(Sample_Size),
    Analysis_ID      = as.factor(Analysis_ID),
    Outcome_Measure  = as.factor(Outcome_Measure),
    Study_Title         = as.character(Study_Title)
  ) %>%
  filter(
    !is.na(Effect_Size),
    !is.na(Variance),
    Variance > 0,
    !is.na(Outcome_Measure)
  )

analysis_grid <- meta_data %>%
  distinct(Analysis_ID, Outcome_Measure)

for (i in seq_len(nrow(analysis_grid))) {
  
  id <- analysis_grid$Analysis_ID[i]
  om <- analysis_grid$Outcome_Measure[i]
  
  id_chr <- as.character(id)
  om_chr <- as.character(om)
  id_safe <- safe_filename(id_chr)
  om_safe <- safe_filename(om_chr)
  
  cat("\n============================================\n")
  cat("Analysis ID:", id_chr, "\n")
  cat("Outcome Measure:", om_chr, "\n")
  cat("============================================\n")
  
  sub_df <- meta_data %>%
    filter(Analysis_ID == id, Outcome_Measure == om)
  
  if (nrow(sub_df) < 2) {
    cat("Not enough studies for quantitative synthesis.\n")
    next
  }
  
  if ("Study_Title" %in% names(sub_df)) {
    sub_df$Study_Label <- ifelse(
      is.na(sub_df$Study_Title),
      sub_df$Study_Title,
      sub_df$Study_Title
    )
  } else {
    sub_df$Study_Label <- sub_df$Study_Title
  }
  sub_df$Study_Label <- as.character(sub_df$Study_Label)
  
  rob_df <- sub_df %>% filter(!is.na(Risk_of_Bias))
  
    print(summary(rob_res))
  }
    cat("\nRisk of bias meta-regression not run ",
        "(insufficient studies or no variation in Risk_of_Bias).\n")
  
  if (nrow(rob_df) >= 2) {
    rob_plot_df <- rob_df %>%
      mutate(
        Risk_Label = factor(
          Risk_of_Bias,
          levels = 1:4,
          labels = c("Low", "Some concerns", "High", "Critical")
        )
      )
    
    p <- ggplot(
      rob_plot_df,
      aes(
        x = Outcome_Measure,
        y = Study_Label,
        fill = Risk_Label,
        text = paste(
          "Study:", Study_Label,
          "<br>Outcome:", Outcome_Measure,
          "<br>Risk:", Risk_Label
        )
      )
    ) +
      geom_tile(color = "white") +
      scale_fill_manual(
        values = c(
          "Low" = "#4CAF50",
          "Some concerns" = "#FFC107",
          "High" = "#FF9800",
          "Critical" = "#F44336"
        )
      ) +
      labs(
        title = paste("Risk of Bias –", id_chr, "|", om_chr),
        x = "Outcome Measure",
        y = "Study"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    htmlwidgets::saveWidget(
      ggplotly(p, tooltip = "text"),
      paste0("risk_of_bias_", id_safe, "_", om_safe, ".html"),
      selfcontained = TRUE
    )
  }
  
  cat("\nCompleted:", id_chr, "|", om_chr, "\n")
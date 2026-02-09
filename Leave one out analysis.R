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
    Study_Title      = as.character(Study_Title)
  ) %>%
  filter(
    !is.na(Effect_Size),
    !is.na(Variance),
    Variance > 0,
    !is.na(Outcome_Measure)
  )

analysis_grid <- meta_data %>%
  distinct(Analysis_ID, Outcome_Measure)

# Loop through each analysis
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
  
  cat("Number of studies in subset:", nrow(sub_df), "\n")
  
  if (nrow(sub_df) < 2) {
    cat("Not enough studies for quantitative synthesis. Skipping...\n")
    next
  }
  
  sub_df$Study_Label <- as.character(sub_df$Study_Title)
  
  leaveout_results <- list()
  
  for (j in 1:nrow(sub_df)) {
    sub_df_excluded <- sub_df[-j, ]
    
    res <- tryCatch({
      rma(
        yi = Effect_Size,
        vi = Variance,
        data = sub_df_excluded,
        method = "REML"
      )
    }, error = function(e) {
      cat("Error fitting meta-analysis model for exclusion of study:", sub_df$Study_Title[j], "\n")
      cat("Error message:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(res)) {
      next
    }
    
    leaveout_results[[j]] <- list(
      estimate = res$b,
      se = res$se,
      zval = res$zval,
      pval = res$pval,
      ci.lb = res$ci.lb,
      ci.ub = res$ci.ub
    )
  }
  
  cat("Number of leave-one-out results:", length(leaveout_results), "\n")
  
  cat("\nLeave-one-out analysis:\n")
  
  study_titles <- sub_df$Study_Title
  
  for (j in 1:length(leaveout_results)) {
    excluded_study <- study_titles[j]
    cat("\nResults for leaving out study:", excluded_study, "\n")
    
    result <- leaveout_results[[j]]
    cat("Estimate:", result$estimate, "\n")
    cat("Standard Error:", result$se, "\n")
    cat("Z-value:", result$zval, "\n")
    cat("P-value:", result$pval, "\n")
    cat("95% CI:", result$ci.lb, "to", result$ci.ub, "\n")
  }
  
  library(ggplot2)
  library(dplyr)
  library(metafor)
  library(stringr)
  
  leaveout_plot_data <- data.frame(
    Omitted_Study = character(),
    Estimate = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    P_Value = numeric(),
    Analysis_ID = character(),
    Outcome_Measure = character()
  )
  
  for (j in 1:length(leaveout_results)) {
    result <- leaveout_results[[j]]
    
    leaveout_plot_data <- rbind(leaveout_plot_data, data.frame(
      Omitted_Study = study_titles[j], 
      Estimate = result$estimate,
      CI_Lower = result$ci.lb,
      CI_Upper = result$ci.ub,
      P_Value = result$pval,
      Analysis_ID = id_chr,
      Outcome_Measure = om_chr
    ))
  }
  
  cat("\nSummary of 95% CI and p-value for each study:\n")
  for (j in 1:nrow(leaveout_plot_data)) {
    cat(paste0(leaveout_plot_data$Omitted_Study[j], ": ", 
               "Estimate = ", round(leaveout_plot_data$Estimate[j], 3), ", ",
               "95% CI = [", round(leaveout_plot_data$CI_Lower[j], 3), ", ", round(leaveout_plot_data$CI_Upper[j], 3), "], ",
               "p-value = ", round(leaveout_plot_data$P_Value[j], 3), "\n"))
  }
  
  forest_plot <- ggplot(leaveout_plot_data, aes(x = Estimate, y = reorder(Omitted_Study, Estimate))) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
    theme_minimal() +
    labs(
      title = paste("Leave-one-out Analysis for", id_chr, "|", om_chr),
      x = "Effect Size (with 95% CI)",
      y = "Omitted Study"
    ) +
    theme(
      axis.text.y = element_text(size = 8), 
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    ) +
    scale_y_discrete(labels = function(x) str_sub(x, 1, 30)) 
  
  ggsave(paste0("leave1out_forest_plot_", id_safe, "_", om_safe, ".png"), plot = forest_plot, width = 8, height = 6)
  
  print(forest_plot)

  cat("\nCompleted:", id_chr, "|", om_chr, "\n")
}

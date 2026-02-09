library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

meta_data <- read_excel("Meta-analysis data (3).xlsx")

required_cols <- c("Effect_Size","Variance","Analysis_ID",
                   "Outcome_Measure","Study_Title")

missing_cols <- setdiff(required_cols, names(meta_data))
if(length(missing_cols)>0) stop("Missing columns: ", paste(missing_cols,collapse=", "))

meta_data <- meta_data %>%
  mutate(
    Effect_Size = as.numeric(Effect_Size),
    Variance = as.numeric(Variance),
    Analysis_ID = as.character(Analysis_ID),
    Outcome_Measure = as.character(Outcome_Measure)
  ) %>%
  filter(!is.na(Effect_Size), !is.na(Variance), Variance > 0)

analysis_grid <- meta_data %>%
  distinct(Analysis_ID, Outcome_Measure)

compute_lfk <- function(TE,se){
  w <- 1/se^2
  pooled <- sum(w*TE)/sum(w)
  z <- (TE-pooled)/se
  left <- mean(z[z<0])
  right <- mean(z[z>=0])
  (right-left)/sd(z)
}

for(i in seq_len(nrow(analysis_grid))){
  
  id <- analysis_grid$Analysis_ID[i]
  om <- analysis_grid$Outcome_Measure[i]
  
  sub_df <- meta_data %>%
    filter(Analysis_ID==id & Outcome_Measure==om)
  
  cat("Processing:", id, "|", om, " | studies:", nrow(sub_df), "\n")
  
  if(nrow(sub_df) < 2){
    cat("Skipped (needs ≥2 studies)\n")
    next
  }
  
  sub_df$se <- sqrt(sub_df$Variance)
  
  lfk_val <- compute_lfk(sub_df$Effect_Size, sub_df$se)
  
  w <- 1/sub_df$se^2
  pooled <- sum(w*sub_df$Effect_Size)/sum(w)
  
  sub_df$z <- (sub_df$Effect_Size-pooled)/sub_df$se
  sub_df$rank <- rank(abs(sub_df$z))
  
  sub_df <- sub_df %>% arrange(rank)
  
  doi_plot <- ggplot(sub_df, aes(rank,z,group=1))+
    geom_line()+
    geom_point(size=3)+
    geom_hline(yintercept=0,linetype="dashed")+
    theme_minimal()+
    labs(
      title=paste0("DOI Plot (LFK=",round(lfk_val,2),") | Analysis ",id," - ",om),
      x="Rank(|Z|)",
      y="Standardized Effect"
    )
  
  filename <- paste0("doi_plot_", make.names(id), "_", make.names(om), ".png")
  
  ggsave(filename, doi_plot, width=7, height=5)
  print(doi_plot)
}

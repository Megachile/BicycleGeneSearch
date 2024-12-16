#!/usr/bin/env Rscript

# train_bicycle_classifier.R
# Trains a GLM classifier to identify bicycle genes based on gene structure characteristics
# Adapted from Han et al. 2023 (https://doi.org/10.1093/molbev/msad210)

# Load required libraries
suppressPackageStartupMessages({
  library(seqinr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(rtracklayer)
  library(gridExtra)
  library(RFLPtools)
  library(rcompanion)
  library(glmnet)
  library(broom)
  library(reshape2)
  library(scales)
  library(tidyverse)
  library(ggseqlogo)
  library(class)
})

# Set up paths
base_dir <- "/export/martinsons/adam/bicycle_classifier"
data_dir <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "output")
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)





generate_diagnostic_plots <- function(model, fit.summary, optimal_cutoff) {
  # Predict on training data
  predictions <- predict(model, fit.summary, type="response")
  pred_df <- data.frame(
    response = predictions,
    Parent = fit.summary$Parent
  )
  
  # 1. Generate response histogram
  p1 <- ggplot() + 
    geom_histogram(aes(predictions), bins = 100) + 
    scale_y_log10() +
    geom_vline(aes(xintercept=optimal_cutoff), color="red") +
    theme_classic() + 
    theme(legend.position = "none", 
          axis.text=element_text(size=16),
          axis.title=element_text(size=18)) + 
    xlab("Response") + 
    ylab("Number of Transcripts")
  
  # 2. Get predictor correlations
  predictor_data <- fit.summary %>% 
    select(total_exon_length, gene_length, first_exon_length, 
           last_exon_length, exon_mean_length, mode0, mode1, mode2)

  sapply(predictor_data, function(x) sum(is.na(x)))
  
  corr_matrix <- cor(predictor_data)
  print(corr_matrix)
  corr_melted <- melt(corr_matrix)
  colnames(corr_melted) <- c("predictor1", "predictor2", "correlation")
  
  # Correlation heatmap
  p2 <- ggplot(data = corr_melted, 
               aes(x=predictor1, y=predictor2, fill=correlation)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#0072B2", high = "#D55E00", 
                        mid = "white", midpoint = 0, 
                        limit = c(-1,1), space = "Lab") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.text=element_text(size=12),
          axis.title=element_blank(),
          legend.title=element_blank(),
          legend.text = element_text(size=11)) +
    coord_fixed()
  
  # Save plots
  ggsave(file.path(output_dir, "response_distribution.pdf"), p1, width=6, height=5)
  ggsave(file.path(output_dir, "predictor_correlations.pdf"), p2, width=8, height=6)
}





calculate_metrics <- function(predictions, fit.summary, cutoff) {
  test.bicycle <- predictions$Parent[predictions$response > cutoff]
  real.bicycle <- fit.summary$Parent[fit.summary$response == 1] 
  test.nonbicycle <- predictions$Parent[predictions$response <= cutoff]
  real.nonbicycle <- fit.summary$Parent[fit.summary$response == 0]  
  
  TP <- length(intersect(test.bicycle, real.bicycle))
  FP <- length(test.bicycle[!test.bicycle %in% real.bicycle])
  TN <- length(intersect(test.nonbicycle, real.nonbicycle))
  FN <- length(test.nonbicycle[!test.nonbicycle %in% real.nonbicycle])
  
  precision <- TP/(TP + FP)
  recall <- TP/(TP + FN)
  
  return(list(precision = precision, recall = recall))
}

# Function to check if required files exist
check_required_files <- function() {
  required_files <- c(
    file.path(data_dir, "reference/hcor_annotations.gff3"),
    file.path(data_dir, "reference/hcor_gene_cluster2.names"),
    file.path(data_dir, "reference/hcor_bicycle_genes.txt")
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Missing required files:\n", paste(missing_files, collapse="\n"))
  }
}

# Function to extract gene features from GFF
extract_gene_features <- function(gff.exons) {
  # General gene statistics
  exons.summary <- gff.exons %>%
    mutate(size = end - start + 1) %>% 
    group_by(Parent) %>%
    summarise(num_total_exons = n(), 
              total_exon_length = sum(size))
  
  # Gene length calculations
  gene.length.summary <- gff.exons %>%
    group_by(Parent) %>%
    summarise(
      gene_start = dplyr::first(start), 
      gene_end = dplyr::last(end),
      gene_length = abs(gene_end-gene_start)+1
    ) %>%
    select(Parent, gene_length)
  
  # First and last exon statistics
  first.exons.summary <- extract_terminal_exons(gff.exons, "first")
  last.exons.summary <- extract_terminal_exons(gff.exons, "last")
  
  # Internal exon statistics
  internal.exons.summary <- extract_internal_exons(gff.exons)
  
  # Merge all features
  gene.summary <- merge(
    merge(
      merge(
        merge(
          exons.summary, 
          gene.length.summary, 
          by="Parent", 
          all=TRUE
        ),
        first.exons.summary, 
        by="Parent", 
        all=TRUE
      ),
      last.exons.summary, 
      by="Parent", 
      all=TRUE
    ),
    internal.exons.summary, 
    by="Parent", 
    all=TRUE
  )
  
  return(gene.summary)
}

# Function to extract terminal exon features
extract_terminal_exons <- function(gff.exons, type) {
  pos_strand <- gff.exons %>% filter(strand == "+")
  neg_strand <- gff.exons %>% filter(strand == "-")
  
  if (type == "first") {
    pos_summary <- pos_strand %>% 
      group_by(Parent) %>% 
      filter(row_number() == 1) %>%
      mutate(size = end - start + 1) %>%
      summarise(first_exon_length = size)
    
    neg_summary <- neg_strand %>% 
      group_by(Parent) %>%
      filter(row_number() == n()) %>%
      mutate(size = end - start + 1) %>%
      summarise(first_exon_length = size)
    
    return(rbind(pos_summary, neg_summary))
  } else {
    pos_summary <- pos_strand %>% 
      group_by(Parent) %>%
      filter(row_number() != 1 & row_number() == n()) %>%
      mutate(size = end - start + 1) %>%
      summarise(last_exon_length = size)
    
    neg_summary <- neg_strand %>% 
      group_by(Parent) %>%
      filter(row_number() == 1 & row_number() != n()) %>%
      mutate(size = end - start + 1) %>%
      summarise(last_exon_length = size)
    
    return(rbind(pos_summary, neg_summary))
  }
}

# Function to extract internal exon features
extract_internal_exons <- function(gff.exons) {
  internal.exons <- gff.exons %>% 
    group_by(Parent) %>% 
    filter(row_number() != 1 & row_number() != n())
  
  internal.exons.summary <- internal.exons %>%
    mutate(size = end - start + 1) %>%
    group_by(Parent) %>%
    summarise(
      num_internal_exons = n(),
      exon_mean_length = mean(size),
      exon_var = var(size),
      mode0 = sum(phase==0),
      mode1 = sum(phase==1),
      mode2 = sum(phase==2)
    )
  
  return(internal.exons.summary)
}

# Function to label genes
label_genes <- function(gene.summary, bicycle_genes) {
  gene.summary %>% 
    mutate(
      label = case_when(
        gsub("\\..*","", Parent) %in% gsub("\\..*","", bicycle_genes) ~ "bicycle",
        TRUE ~ "non-bicycle"
      )
    )
}

# Function to train GLM model
train_glm_model <- function(fit.summary) {
  glm.full <- glm(
    response ~ total_exon_length + gene_length + first_exon_length + 
                last_exon_length + exon_mean_length + mode0 + mode1 + mode2,
    data = fit.summary,
    family = binomial(link="logit")
  )
  return(glm.full)
}

# Function to find optimal cutoff
find_optimal_cutoff <- function(model, fit.summary) {
  test.predict <- predict(model, fit.summary, type="response")
  test.predict.merged <- cbind.data.frame(
    fit.summary$Parent,
    test.predict
  )
  colnames(test.predict.merged) <- c("Parent", "response")
  
  cutoffs <- seq(0, 0.999, 0.02)
  results <- data.frame(
    cutoff = cutoffs,
    precision = NA,
    recall = NA
  )
  
  for (i in seq_along(cutoffs)) {
    metrics <- calculate_metrics(
      test.predict.merged,
      fit.summary,
      cutoffs[i]
    )
    results$precision[i] <- metrics$precision
    results$recall[i] <- metrics$recall
  }
  
  results$ratio <- results$precision / results$recall
  optimal_cutoff <- results$cutoff[which.min(abs(results$ratio - 1))]
  
  return(optimal_cutoff)
}



format_gff <- function(gff) {
	gff$Parent <- unlist(gff$Parent)
	parents <- unique(gff$Parent)
	gff$Parent <- factor(gff$Parent, levels=parents)
	return(data.frame(gff))
}


# Main execution
main <- function() {
  # Check required files exist
  check_required_files()
  
  # Read and process input files
  message("Reading GFF file...")
  Hcor <- readGFF(
    file.path(data_dir, "reference/hcor_annotations.gff3"),
    filter=list(type="CDS")
  )
  
  # Format GFF
  gff.exons <- format_gff(Hcor)
  
  # Extract features
  message("Extracting gene features...")
  gene.summary <- extract_gene_features(gff.exons)
  
  # Label genes
  message("Labeling genes...")
  bicycle_genes <- read.table(
    file.path(data_dir, "reference/hcor_gene_cluster2.names"),
    stringsAsFactors=TRUE
  )$V1
  
  gene.summary <- label_genes(gene.summary, bicycle_genes)
  
  # Prepare training data
  message("Preparing training data...")
  fit.summary <- gene.summary %>% 
    filter(label != "unannotated") %>%
    mutate(response = if_else(label == "bicycle", 1, 0))
  
  # Train model
  message("Training GLM model...")
  model <- train_glm_model(fit.summary)
  
  # Find optimal cutoff
  message("Finding optimal cutoff...")
  optimal_cutoff <- find_optimal_cutoff(model, fit.summary)
  
  # Save model and results
  message("Saving results...")
  save(model, file=file.path(output_dir, "bicycle_classifier.rda"))
  write.table(
    data.frame(cutoff=optimal_cutoff),
    file=file.path(output_dir, "optimal_cutoff.txt"),
    row.names=FALSE
  )
  
  # Generate plots
  message("Generating diagnostic plots...")
  generate_diagnostic_plots(model, fit.summary, optimal_cutoff)
  
  message("Done! Model and results saved to ", output_dir)
}

if (!interactive()) {
  main()
}

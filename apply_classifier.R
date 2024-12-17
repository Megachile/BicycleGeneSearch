#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rtracklayer)
})

# Configuration
INPUT_GTF <- "/export/martinsons/adam/bicycle_classifier/data/raw/bkinseyi_annotations.gtf"
INPUT_MODEL <- "/export/martinsons/adam/bicycle_classifier/output/bicycle_classifier.rda"
OUTPUT_DIR <- "/export/martinsons/adam/bicycle_classifier/output/bkin_results"
OUTPUT_PREFIX <- "bkin"
CUTOFF <- 0.68  # Set directly from optimal_cutoff.txt

# Function to predict bicycle genes
predict_bicycle <- function(gtf.file, model.file, cutoff, output_prefix, output_dir) {
  # Read GTF with rtracklayer, filtering for CDS features
  message("Reading GTF file...")
  gff.cds <- import(gtf.file)
  gff.cds <- gff.cds[gff.cds$type == "CDS"]
  
  # Load model
  message("Loading model...")
  e <- new.env()
  load(model.file, envir=e)
  model.name <- ls(e)[1]  # Get the first object name from the loaded environment
  glm.model <- get(model.name, envir=e)
  
  message("Processing GTF data...")
  # Convert to data frame and ensure proper Parent/transcript_id handling
  gff.cds <- as.data.frame(gff.cds)
  if ("transcript_id" %in% colnames(gff.cds)) {
    gff.cds$Parent <- gff.cds$transcript_id
  }
  
  # Format Parent IDs
  parents <- unique(gff.cds$Parent)
  gff.cds$Parent <- factor(gff.cds$Parent, levels = parents)
    
  # Calculate gene statistics
  message("Calculating gene statistics...")
  exons.summary <- gff.cds %>%
    mutate(size = end - start + 1) %>% 
    group_by(Parent) %>%
    summarise(num_total_exons = n(), total_exon_length = sum(size))
  
  # Gene length calculations  
  gene.length.summary <- gff.cds %>%
    group_by(Parent) %>%
    summarise(
      gene_start = min(start),
      gene_end = max(end),
      gene_length = abs(max(end) - min(start)) + 1
    ) %>%
    select(Parent, gene_length)
    
  # Process strand-specific features
  message("Processing exon features...")
  pos_strand <- gff.cds %>% filter(strand == "+")
  neg_strand <- gff.cds %>% filter(strand == "-")
  
  # Calculate first/last exon lengths for positive strand
  pos.first.exons <- pos_strand %>% 
    group_by(Parent) %>% 
    arrange(start) %>%
    filter(row_number() == 1) %>%
    mutate(size = end - start + 1) %>%
    summarise(first_exon_length = size)
    
  pos.last.exons <- pos_strand %>% 
    group_by(Parent) %>%
    arrange(start) %>%
    filter(row_number() != 1 & row_number() == n()) %>%
    mutate(size = end - start + 1) %>%
    summarise(last_exon_length = size)
    
  # Calculate first/last exon lengths for negative strand  
  neg.first.exons <- neg_strand %>% 
    group_by(Parent) %>%
    arrange(desc(start)) %>%
    filter(row_number() == 1) %>%
    mutate(size = end - start + 1) %>%
    summarise(first_exon_length = size)
    
  neg.last.exons <- neg_strand %>% 
    group_by(Parent) %>%
    arrange(desc(start)) %>%
    filter(row_number() != 1 & row_number() == n()) %>%
    mutate(size = end - start + 1) %>%
    summarise(last_exon_length = size)
  
  # Combine strand-specific calculations
  first.exons.summary <- rbind(pos.first.exons, neg.first.exons)
  last.exons.summary <- rbind(pos.last.exons, neg.last.exons)
    
  # Process internal exons
  message("Processing internal exons...")
  internal.exons <- gff.cds %>% 
    group_by(Parent) %>%
    arrange(start) %>%
    filter(row_number() != 1 & row_number() != n())
  
  internal.exons.summary <- internal.exons %>%
    mutate(size = end - start + 1) %>%
    group_by(Parent) %>%
    summarise(
      num_internal_exons = n(),
      exon_mean_length = mean(size),
      exon_var = var(size),
      mode0 = sum(phase == 0),
      mode1 = sum(phase == 1),
      mode2 = sum(phase == 2)
    )
  
  # Combine all statistics
  message("Combining statistics...")
  gene.summary <- merge(
    merge(
      merge(
        merge(
          exons.summary,
          gene.length.summary,
          by = "Parent",
          all = TRUE
        ),
        first.exons.summary,
        by = "Parent",
        all = TRUE
      ),
      last.exons.summary,
      by = "Parent",
      all = TRUE
    ),
    internal.exons.summary,
    by = "Parent",
    all = TRUE
  )
  
  # Remove genes with 3 or fewer exons
  gene.summary <- gene.summary[complete.cases(gene.summary), ]
  
  # Run classifier
  message("Running classifier...")
  predictions <- predict(glm.model, gene.summary, type="response")
  results <- data.frame(
    Parent = gene.summary$Parent,
    response = predictions
  )

  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Saving results...")
  # Save all results
  write.table(
    results,
    file.path(output_dir, paste0(output_prefix, "_classifier_all_transcripts_response.txt")),
    quote = FALSE,
    row.names = FALSE,
    sep = '\t'
  )
    
  # Filter and save bicycle genes
  bicycle.genes <- results[results$response > cutoff, ]
  bicycle.gene.names <- unique(gsub("\\.[^.]*$", "", bicycle.genes$Parent))
  write.table(
    bicycle.gene.names,
    file.path(output_dir, paste0(output_prefix, "_classifier_bicycle_genes.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # Generate response distribution plot
  message("Generating plot...")
  pdf(
    file.path(output_dir, paste0(output_prefix, "_classifier_response_histogram.pdf")),
    width = 4,
    height = 4
  )
  print(ggplot() +
    geom_histogram(aes(results$response), bins = 100) +
    scale_y_log10() +
    geom_vline(aes(xintercept = cutoff), color = "red") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18)
    ) +
    xlab("Response") +
    ylab("Transcript count")
  )
  dev.off()
  
  message("Analysis complete!")
  message(sprintf("Found %d potential bicycle genes", length(bicycle.gene.names)))
}

# Main execution
message("Starting bicycle gene classification...")
predict_bicycle(INPUT_GTF, INPUT_MODEL, CUTOFF, OUTPUT_PREFIX, OUTPUT_DIR)

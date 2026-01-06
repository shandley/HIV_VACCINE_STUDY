# LEfSe Minimal Functions for Microbiome Biomarker Discovery
# Author: Claude
# Date: 2025-04-18

library(microbiomeMarker)
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

#' Run LEfSe analysis and return results
#'
#' @param ps Phyloseq object
#' @param group_var Column name in sample_data defining groups
#' @param taxa_rank Taxonomic rank to aggregate at (default: "Genus")
#' @param kw_cutoff Kruskal-Wallis test p-value cutoff (default: 0.05)
#' @param wilcoxon_cutoff Wilcoxon test p-value cutoff (default: 0.05)
#' @param lda_cutoff LDA score threshold (default: 2.0)
#' @return LEfSe results object
run_lefse_analysis <- function(ps, group_var, 
                              taxa_rank = "Genus",
                              kw_cutoff = 0.05,
                              wilcoxon_cutoff = 0.05,
                              lda_cutoff = 2.0) {
  
  # Run LEfSe analysis
  lefse_results <- run_lefse(
    ps,
    group = group_var,
    taxa_rank = taxa_rank,
    kw_cutoff = kw_cutoff,
    wilcoxon_cutoff = wilcoxon_cutoff,
    lda_cutoff = lda_cutoff
  )
  
  return(lefse_results)
}

#' Extract significant taxa from LEfSe results as a dataframe
#'
#' @param lefse_results LEfSe results object
#' @param lda_cutoff Minimum LDA score to consider significant (default: 2.0)
#' @param output_file Optional CSV output filename
#' @return Dataframe with significant taxa
extract_significant_taxa <- function(lefse_results, lda_cutoff = 2.0, output_file = NULL) {
  
  # Extract data frame of marker taxa
  marker_table <- data.frame(
    feature = lefse_results@marker_table$feature,
    enrich_group = lefse_results@marker_table$enrich_group,
    ef_lda = lefse_results@marker_table$ef_lda,
    pvalue = lefse_results@marker_table$pvalue,
    stringsAsFactors = FALSE
  )
  
  # Filter for significant taxa
  significant_taxa <- marker_table %>%
    filter(ef_lda >= lda_cutoff) %>%
    rename(
      Genus = feature,
      LDA_score = ef_lda,
      Enriched_in = enrich_group
    )
  
  # Calculate adjusted p-values
  significant_taxa$padj <- p.adjust(significant_taxa$pvalue, method = "BH")
  
  # Reorder columns
  significant_taxa <- significant_taxa %>%
    select(Genus, LDA_score, Enriched_in, pvalue, padj)
  
  # Save to CSV if output file is specified
  if (!is.null(output_file)) {
    write.csv(significant_taxa, output_file, row.names = FALSE)
    cat(sprintf("Results saved to %s\n", output_file))
  }
  
  return(significant_taxa)
}

#' Create LEfSe plot with LDA scores only and store in variable
#'
#' @param lefse_results LEfSe results object
#' @param ps Phyloseq object
#' @param group_var Column name in sample_data defining groups
#' @param colors Named vector of colors for groups
#' @param output_file Optional filename to save the plot
#' @param width Width of saved plot (default: 8)
#' @param height Height of saved plot (default: 10)
#' @param variable_name Base name for the plot variable (will be prefixed with "lefse_plot_")
#' @param save_as_tiff Logical, whether to save as TIFF (default: FALSE)
#' @param dpi Resolution for TIFF file (default: 300)
#' @param compression Compression method for TIFF ("none", "lzw", "jpeg", etc.; default: "lzw")
#' @return LDA score plot object
create_lefse_plot <- function(lefse_results, ps, group_var, 
                             colors = NULL, 
                             output_file = NULL,
                             width = 8, height = 10,
                             variable_name = "default",
                             save_as_tiff = FALSE,
                             dpi = 300,
                             compression = "lzw") {
  
  # Set default colors if not provided
  if (is.null(colors)) {
    # Get unique groups
    groups <- unique(sample_data(ps)[[group_var]])
    
    # Create default colors
    default_colors <- c("#4DAF4A", "#E41A1C", "#377EB8", "#984EA3", "#FF7F00")
    colors <- default_colors[1:length(groups)]
    names(colors) <- groups
  }
  
  # Create the LDA score bar plot
  p_lda <- plot_ef_bar(lefse_results) +
    theme_bw() +
    theme(axis.text.y = element_text(face = "italic")) +  # Italicize genus names
    scale_fill_manual(values = colors) +
    labs(title = "LEfSe Results")
  
  # Generate the full variable name
  full_var_name <- paste0("lefse_plot_", variable_name)
  
  # Save the plot if output file is specified
  if (!is.null(output_file)) {
    if (save_as_tiff) {
      # If filename doesn't end with .tiff or .tif, append .tiff
      if (!grepl("\\.tiff?$", output_file, ignore.case = TRUE)) {
        output_file <- paste0(output_file, ".tiff")
      }
      
      # Save as TIFF with compression
      ggsave(
        output_file, 
        p_lda, 
        width = width, 
        height = height, 
        dpi = dpi,
        device = "tiff",
        compression = compression
      )
      cat(sprintf("Plot saved as compressed TIFF to %s (dpi: %d, compression: %s)\n", 
                 output_file, dpi, compression))
    } else {
      # Standard save
      ggsave(output_file, p_lda, width = width, height = height)
      cat(sprintf("Plot saved to %s\n", output_file))
    }
  }
  
  # Print the variable name for reference
  cat(sprintf("Plot stored in variable: %s\n", full_var_name))
  
  # Return the plot
  return(p_lda)
}

#' Extract abundance data from phyloseq for significant LEfSe taxa
#'
#' @param lefse_results LEfSe results object
#' @param ps Phyloseq object
#' @param group_var Column name in sample_data defining groups
#' @param output_file Optional CSV output filename
#' @return Dataframe with abundance data for significant taxa
extract_abundance_data <- function(lefse_results, ps, group_var, output_file = NULL) {
  
  # Extract marker information and sort by LDA score
  markers <- marker_table(lefse_results)
  genus_markers <- markers[order(-abs(markers$ef_lda)), ]  # Sort by LDA score
  genus_order <- genus_markers$feature  # Get ordered genus list
  
  # Prepare abundance data
  ps_genus <- tax_glom(ps, taxrank = "Genus")
  ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x*100/sum(x))
  ps_data <- psmelt(ps_genus_rel)
  
  # Filter to only include significant genera
  ps_data_sig <- ps_data[ps_data$Genus %in% genus_order, ]
  
  # Calculate summary statistics
  abundance_summary <- ps_data_sig %>%
    group_by(Genus, !!sym(group_var)) %>%
    summarize(
      mean_abundance = mean(Abundance, na.rm = TRUE),
      median_abundance = median(Abundance, na.rm = TRUE),
      min_abundance = min(Abundance, na.rm = TRUE),
      max_abundance = max(Abundance, na.rm = TRUE),
      sd_abundance = sd(Abundance, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_abundance))
  
  # Save to CSV if output file is specified
  if (!is.null(output_file)) {
    write.csv(abundance_summary, output_file, row.names = FALSE)
    cat(sprintf("Abundance data saved to %s\n", output_file))
  }
  
  return(abundance_summary)
}
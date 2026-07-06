make_custom_plot <- function(gene_name, subset_cells = NULL, show_x_axis = FALSE) {
  
  # Safe runtime validation check: Fall back to parent object if barcodes are empty
  if (is.null(subset_cells)) {
    obj_subset <- refmap.atac
  } else {
    # Inline subsetting isolates text/expression limits to selected cell barcodes
    obj_subset <- subset(refmap.atac, cells = subset_cells)
  }
  
  # 1. Main Coverage Signal Track
  track_cov <- CoveragePlot(
    object = obj_subset, 
    region = gene_name,
    group.by = "class", 
    ymax = "q95"
  ) & scale_fill_manual(values = c("aparasitemic" = "#58a4cb", "parasitemic" = "#e41a1c")) &
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  # 2. Gene Expression (Violin/Boxplot) Track 
  # Now running on the subset object to avoid duplicate genomic timeline bugs
  track_expr <- ExpressionPlot(
    object = obj_subset, 
    features = gene_name,
    group.by = "class", 
    assay = "RNA"
  ) + scale_fill_manual(values = c("aparasitemic" = "#58a4cb", "parasitemic" = "#e41a1c")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  # 3. Genomic Links Loop Track
  track_links <- LinkPlot(object = obj_subset, region = gene_name) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  # 4. Gene Annotation Structure Track
  track_annot <- AnnotationPlot(object = obj_subset, region = gene_name)
  
  # Handle upper vs lower track x-axis timeline coordinates explicitly
  if (!show_x_axis) {
    track_annot <- track_annot + 
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    track_annot <- track_annot + 
      theme(axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 7))
  }
  
  # Compile track elements into localized patchwork layout blocks
  compiled_layout <- (track_cov | track_expr) / track_links / track_annot + 
    plot_layout(heights = c(8, 1.5, 1), widths = c(4, 1))
  return(compiled_layout)
}

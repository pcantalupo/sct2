#' SignacInfo
#' @description
#' Show information about the Signac ChromatinAssays within a Seurat object
#' @param seurat A Seurat object
#' @param assay Assay to get information for (default is all ChromAssays)
#' @param show_annotation Show genome annotation information (default FALSE)
#' @export
#' @importFrom Signac Fragments Motifs Links Annotation granges seqinfo
#' @importFrom SeuratObject Assays
#' @importFrom GenomicRanges width mcols
#'
#' @examples
#' SignacInfo(multiome_small)
SignacInfo = function(seurat, assay = NULL, show_annotation = FALSE) {

  # Validate assay ------------------------------------------------
  if (is.null(assay)) {
    assays = SeuratObject::Assays(seurat)
  } else {
    if (!assay %in% names(seurat))
      stop("Assay '", assay, "' not present in the Seurat object.")
    assays = assay
  }

  # Filter for ChromatinAssays
  assays = Filter(  # base function
    function(a) inherits(seurat[[a]], "ChromatinAssay"),
    assays
  )
  if (!length(assays))
    stop("No Signac ChromatinAssay found for the supplied request.")


  # Function to summarize one ChromAssay --------------------------
  ChromatinAssayInfo = function(chrom_assay, assay) {

    # Core numbers -------------------------------------------------------
    n_cells    = ncol(chrom_assay)
    n_features = nrow(chrom_assay)
    peak_widths = GenomicRanges::width(Signac::granges(chrom_assay))
    width_summary = summary(peak_widths)
    width_min = width_summary[1]
    width_med = width_summary[3]
    width_max = width_summary[6]
    width_range = paste0(width_min, "-", width_max)


    # Annotation summary ---------------------------------------------------
    anno = Signac::Annotation(chrom_assay)
    genome = "unknown"; n_chrom = NA    # defaults
    anno_summary = NULL
    if (!is.null(anno)) {
      seqinfo = Signac::seqinfo(anno)
      n_chrom = length(seqinfo@seqnames)
      genome = paste0(unique(seqinfo@genome), collapse = ", ")

      md = GenomicRanges::mcols(anno)
      anno_summary = list(
        ranges         = length(anno),
        type_counts    = if ("type" %in% colnames(md))  sort(table(md$type),  decreasing = TRUE) else "column 'type' absent",
        gene_biotypes  = if ("gene_biotype" %in% colnames(md)) sort(table(md$gene_biotype), decreasing = TRUE) else "column 'gene_biotype' absent",
        strand_counts  = table(anno@strand),
        chrom_counts   = sort(table(anno@seqnames), decreasing = TRUE)[1:10],
        seqinfo        = seqinfo
      )
    }


    # Fragments ------------------------------------------------------------
    frags = Signac::Fragments(chrom_assay)
    n_frags = length(frags)
    frag_info = lapply(frags, function(frag) {
      ncells = length(frag@cells)
      path = frag@path
      list(ncells = ncells, path = path)
    })
    frag_numcells = sapply(frag_info, function(frag) {frag$ncells})
    frag_paths = sapply(frag_info, function(frag) {frag$path})


    # Motifs / Links / RegionStats / Positional enrichment -----------------
    motifs  = Signac::Motifs(chrom_assay)
    motif_msg = if (is.null(motifs)) { "None" } else {
      paste0("Motif object with ", length(motifs), " motifs")
    }

    links = Signac::Links(chrom_assay)
    n_links = if (is.null(links)) 0 else length(links)
    max_link_score = NA; min_link_score = NA
    if (n_links > 0) {
      max_link_score = max(links$score, na.rm = TRUE)
      min_link_score = min(links$score, na.rm = TRUE)
    }
    links_msg = paste0(n_links, " links (max: ", round(max_link_score,2),
                       "; min: ", round(min_link_score, 2), ")")

    regionstats = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                     "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT",
                     "GC.percent", "sequence.length")
    region_msg = if (all(regionstats %in% colnames(chrom_assay@meta.features))) {
      "RegionStats has been run already"
    } else { "RegionStats has NOT been run yet" }

    pos_enrich = chrom_assay@positionEnrichment
    pos_terms  = if (length(pos_enrich) == 0) { "None" } else {
      paste(names(pos_enrich), collapse = ", ")
    }


    # Print results ------------------------------------------------------
    cat("\n==============================\n",
        "ChromatinAssay summary (", assay, ")\n",
        "==============================\n", sep = "")
    cat("Cells:            ", n_cells, "\n",
        "Features:         ", n_features, "\n",
        "Width range:      ", width_range, "\n",
        "Width median:     ", as.numeric(width_med), "\n",
        "Genome:           ", genome, "\n",
        "Chromosomes:      ", n_chrom, "\n", sep = "")

    cat("\n----- Fragment objects (", n_frags, ") -----\n", sep = "")
    if (n_frags) {
      for (i in 1:length(frag_info)) {
        cat(paste0(frag_info[[i]]$path, " (", frag_info[[i]]$ncells, " cells)", "\n"))
      }
    } else { cat("None\n") }

    cat("\n----- Motifs -----\n", motif_msg, "\n", sep = "")
    cat("\n----- Links  -----\n", links_msg, "\n", sep = "")
    cat("\n----- Region stats  -----\n", region_msg, "\n", sep = "")
    cat("\n----- Positional enrichment -----\n", pos_terms, "\n", sep = "")

    if (show_annotation) {
      cat("\n----- Annotation -----\n")
      if (is.null(anno_summary)) {
        cat("No annotation stored\n")
      } else {
        cat("Ranges: ", anno_summary$ranges, "\n", sep = "")
        cat("\nGene 'type' counts:"); print(anno_summary$type_counts)
        cat("\nGene biotypes:"); print(anno_summary$gene_biotypes)

        #cat("\nStrand distribution:"); print(anno_summary$strand_counts)
        #cat("\nMost frequent chromosomes:"); print(anno_summary$chrom_counts)
      }
    }


    # Create list for return -------------------------------------------------
    list(
      Assay                  = assay,
      Cells                  = n_cells,
      Features               = n_features,
      Width_Range            = width_range,
      Width_Median           = as.numeric(width_med),
      Genome                 = genome,
      Chromosomes            = n_chrom,
      Fragments = list(
        n_fragment_objects = n_frags,
        fragment_num_cells = frag_numcells,
        fragment_paths     = frag_paths
      ),
      Motifs                 = motif_msg,
      Links = list(
        n_links   = n_links,
        max_score = max_link_score,
        min_score = min_link_score
      ),
      Region_Stats           = region_msg,
      Positional_Enrichment  = pos_terms,
      Annotation_Summary     = anno_summary
    )
  }

  # Iterate over each ChromatinAssay and collect results
  results = setNames(vector("list", length(assays)), assays)
  for (assay in assays) {
    results[[assay]] = ChromatinAssayInfo(seurat[[assay]], assay)
  }

  invisible(results)
}



#' @importFrom tibble rownames_to_column
#' @importFrom readr read_tsv write_tsv
#' @importFrom dplyr arrange case_when filter group_by left_join mutate desc
#' @importFrom dplyr n_distinct select summarise %>%
#' @importFrom ggplot2 ggplot facet_wrap stat_summary ylab xlab theme theme_bw
#' @importFrom ggplot2 scale_x_continuous ggtitle ggsave aes geom_point waiver
#' @importFrom ggplot2 scale_y_continuous scale_y_log10 scale_colour_manual
#' @importFrom ggplot2 element_blank scale_shape_manual geom_function units
#' @importFrom patchwork wrap_plots
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq plotPCA vst results counts
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette dev.off pdf svg
#' @importFrom stats complete.cases cor dnorm kmeans time
#' @importFrom utils head
#' @importFrom grid grid.grabExpr gpar
#' @importFrom ggrepel geom_label_repel
#' @importFrom ComplexHeatmap Heatmap anno_mark draw
#' @importFrom purrr map map_lgl
#' @importFrom graphics plot.new


#' @title Run a full RNASeqhelper analysis
#' @description Single function to call all other functions in the
#'     analysis and perform a full analysis from a few starting
#'     parameters
#' @param tab a matrix of samples (columns) and informative genes
#'     (rows), ideally subset using the function `high_quality_genes'.
#' @param phenotype_data a table with samples (rows) with extra
#'     columns for annotation groups. One of these groups must be
#'     Condition and will be used to normalize the data.
#' @param keep_params A list of parameters to override the default
#'     `high_quality_genes' function. If \code{NULL}, then use the
#'     default.
#' @param heat_params A list of parameters to override the default
#'     `pairwise_heatmap_volcano' function. If \code{NULL}, then use
#'     the default.
#' @param gene list of parameters to override the default
#'     `gene_clusters_by_score' function. If \code{NULL}, then use the
#'     default.
rnaseqhelper <- function(tab, phenotype,
                         keep_params = NULL, heat_params = NULL,
                         gcbs_params = NULL) {

  ## out_dir="test/1_genes"
  keep_genes <- do.call(high_quality_genes, keep_params)
  
  res <- run_deseq(tab, keep_genes, phenotype_data)

  ## out_dir="test/1_genes"
  pca_and_matrices(res, out_dir = "test/1_matrices_and_deseq")

  ## res$ddsObj
  ## out_dirprefix="test/outputs"
  heat_params_defaults <- formals(pairwise_hmap_volcano)
  heat_params <- modifyList(heat_params_defaults, heat_params) ## override
  heat_params$transformed_counts <- res$vFalse
  heat_params$out_dirprefix <- "test/2_heatmaps"
  phv <- do.call(pairwise_hmap_volcano, heat_params)
}


#' @title Run DESeq with sensible defaults
#' @description Runs DESeq with filtering thresholds and shows PCA for
#'   variance-stabilize-transformed data.
#' @param tab a matrix of samples (columns) and informative genes
#'   (rows), ideally subset using the function `high_quality_genes'.
#' @param keep_genes a vector of genes to subset.
#' @param phenotype_data a table with samples (rows) with extra
#'   columns for annotation groups. One of these groups must be
#'   Condition and will be used to normalize the data.
#' @return list of tables
#' @export
run_deseq <- function(tab, keep_genes, phenotype_data) {

  if (!("Condition" %in% colnames(phenotype_data))) {
    stop("Could not find `Condition' column in phenotype data")
  }

  sub_as <- tab[keep_genes, ]
  
  ddsObj <- DESeqDataSetFromMatrix(
    countData = sub_as,
    colData = phenotype_data,
    design = ~Condition
  )
  ddsObj <- DESeq(ddsObj)

  vsd1 <- vst(ddsObj, blind=FALSE)
  p1 <- plotPCA(vsd1, intgroup = c("Condition"), returnData = F )
  vsd2 <- vst(ddsObj, blind=TRUE)
  p2 <- plotPCA(vsd2, intgroup = c("Condition"), returnData = F )
  return(list(tab = tab, phenotype = phenotype_data,
              bFalse = p1, bTrue = p2,
              vFalse = vsd1, vTrue = vsd2,
              ddsObj = ddsObj, sub = sub_as))
}

#' @title Print and Store PCA and Matrices
#' @description Takes the output of a DESeq2 run and generates
#'   normalization tables and PCAs
#' @param res Output of `run_deseq`.
#' @param out_dir String depicting the output directory to place
#'   tables and plots.
#' @return None. Plots and tables and placed into output directory.
pca_and_matrices <- function(res, out_dir = "deseq2") {
  ## Initial input matrices
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  res$tab %>% as.data.frame %>%
    rownames_to_column("gene") %>%
    write_tsv(file.path(out_dir, "input_matrix.tsv"))

  res$phenotype %>% write_tsv(file.path(out_dir, "phenotype_data.tsv"))

  saveRDS(res, file.path(out_dir, "deseq2obj.rds"))

  res$sub %>% as.data.frame %>%
    rownames_to_column("gene") %>%
    write_tsv(file.path(out_dir, "input_matrix.filt.tsv"))

  counts(res$ddsObj, normalized=TRUE) %>% as.data.frame %>%
    rownames_to_column("gene") %>%
    write_tsv(file.path(out_dir, "input_matrix.filt.normalized.tsv"))

  assay(res$vFalse) %>% as.data.frame %>% rownames_to_column("gene") %>%
    write_tsv(
      file.path(out_dir, "input_matrix.filt.normalized.vst_corrected.tsv"))

  pdf(file.path(out_dir, "input_matrix.filt.normalized.vst_corrected_PCA.pdf"),
      width = 7, height = 7)
  print(res$bFalse)
  print(res$bTrue)
  dev.off()
}

#' @title Pairwise heatmaps and volcano plots
#' @description For a given contrast and DESeq object produce volcano
#'     plots of the DE genes, and generate a global heatmap with those
#'     DE genes and a pairwise one with just the numerator and
#'     denominator samples.
#' @param ddsObj A DESeq object,
#' @param transformed_counts REALLY UNSURE WHAT'S GOING ON HERE
#' @param numer A string prefix to select sample columns that will be
#'     part of the numerator part of the contrast.
#' @param denom A string prefix to select sample columns that will
#'     part of the denominator part of the contrast.
#' @param top_ngenes_tocluster A positive integer. How many of the top
#'     log10 padj values to take. Default is 2000.
#' @param top_ngenes_tohighlight A positive integer. How many of the
#'     top log10 padj values to label in the plots. Default is 50.
#' @param lFC_zoom A number to depict the Log2FC threshold of the
#'     zoomed plot.
#' @param pAdj_zoom A number to depict the adjusted p-value threshold
#'     of the zoomed plot.
#' @param kmeans A list? of positive integers to do k-means
#' @param out_dirprefix A character prefix outlining the directory and
#'     basename in which plots and tables will be deployed.
#' @return Void. Plots are deposited to the output directory.
#' @export
pairwise_hmap_volcano <- function(ddsObj,
                                  transformed_counts = NULL,
                                  numer = "FRT", denom = "TAF2",
                                  top_ngenes_tocluster = 2000,
                                  top_ngenes_tohighlight = 50,
                                  lFC_zoom = 1.5, pAdj_zoom = 20,
                                  kmeans = 2, out_dirprefix = ".") {

  message("Started Analysis: ", date())
  
  plot_title <- paste0(numer, " vs ", denom)
  outdir <- file.path(out_dirprefix,
                      tolower(gsub("[^A-Za-z0-9]", "_", plot_title)))

  dir.create(outdir, recursive=T, showWarnings=FALSE)

  ## 1. Perform RNA-seq contrast between numerator and denominator
  dsqres <- results(ddsObj, contrast = c("Condition", numer, denom),
                    cooksCutoff = Inf,
                    independentFiltering = FALSE) %>%
    as.data.frame %>% rownames_to_column("gene") %>%
    mutate(mLog10Padj = -log10(padj))

  norm_counts <- counts(ddsObj, normalized = TRUE)

  top_genes_tocluster <- (dsqres %>% arrange(desc(mLog10Padj)) %>%
                          head(top_ngenes_tocluster))$gene
  top_genes_tohighlight <- (dsqres %>% arrange(desc(mLog10Padj)) %>%
                            head(top_ngenes_tohighlight))$gene

  write_tsv(data.frame(top_genes_tocluster=top_genes_tocluster),
            file.path(outdir, paste0("clustered_genes.top",
                                     top_ngenes_tocluster, ".tsv")))
  write_tsv(data.frame(top_genes_tohighlight = top_genes_tohighlight),
            file.path(outdir, paste0("volcano_genes.tophighlight",
                                     top_ngenes_tohighlight, ".tsv")))
  
  sample_columns <- c(grep(paste0("^", numer),
                           colnames(norm_counts), value = TRUE),
                      grep(paste0("^", denom),
                           colnames(norm_counts), value = TRUE))
  ## Volcano Plots
  p1 <- volcano_plot(dsqres,
                     top_genes_tohighlight,
                     plot_title,
                     curve = list(sd = 0.3, sc = 60, offset = 10),
                     curve_show = F)
  ## Volcano Plots zoomed in
  p2 <- volcano_plot(dsqres %>% filter(abs(log2FoldChange) < lFC_zoom &
                                       mLog10Padj < pAdj_zoom),
                     top_genes_tohighlight,
                     paste0(plot_title, " (zoomed)"),
                     curve = list(sd = 0.25, sc = 5, offset = 8),
                     ylim = c(0,22),
                     curve_show = F)

  hmaps <- wrap_plots(list(p1, p2), ncol = 1, guides = "collect") &
    theme_bw() +
    theme(legend.position = "none")

  volcano_svg <- file.path(outdir, "volcano_pairwise.svg")
  deseq2_out <- file.path(outdir, "deseq2_results.tsv")

  svg(volcano_svg, width = 8, height = 9)
  print(hmaps)
  dev.off()
  message("Saved Volcano: ", volcano_svg)

  dsqres %>%
    mutate(isTopN.gene = gene %in% top_genes_tohighlight) %>%
    write_tsv(deseq2_out)
  message("Saved DESeq2 Results: ", deseq2_out)

  for (kmk in kmeans) {
    message("[Running k=", kmk, "]")

    ## Pairwise Heatmap of Normalised and Transformed Counts
    nice_kmeans_heatmap_norm_and_trans(
      norm_counts, transformed_counts,
      genes_tocluster = top_genes_tocluster,
      sample_columns = sample_columns,
      genes_tohighlight = top_genes_tohighlight,
      dsqres, kmk,
      out_dir = file.path(outdir, paste0("kmeans", kmk)),
      heatprefix = "heatmap_pairwise",
      plot_title = "Heatmap Pairwise")

    ## This below will produce different scaled matrices than the one
    ## above. We can't really pass these back. All we can do is
    ## incorporate the gene plotting stuff into the tail end of this
    ## pipeline, meaning there will be gene plots for each kmeans heatmap.
    
    ## Global Heatmap of Normalised and Transformed Counts
    nice_kmeans_heatmap_norm_and_trans(
      norm_counts, transformed_counts,
      genes_tocluster = top_genes_tocluster,
      sample_columns = NULL,
      genes_tohighlight = top_genes_tohighlight,
      dsqres, kmk,
      out_dir = file.path(outdir, paste0("kmeans", kmk)),
      heatprefix = "heatmap_all",
      plot_title = "Heatmap All")
  }
  NULL
  message("Finished Analysis: ", date())
}


#' @title Nice K-means Heatmap
#' @description Produce a nice K-means clustered heatmap for both
#'   normalized and transformed (DESeq2 VST corrected) matrices
#' @param norms a matrix of normalized counts
#' @param trans a matrix of transformed counts, VST or other produced
#'   by DESeq2.
#' @param genes_to_cluster a vector of gene names. If NULL (the
#'   default), then cluster all genes.
#' @param sample_columns a vector of sample names. If NULL (the
#'   default), then use all samples
#' @param dsqres the output of the DESeq2 function \code{results},
#'   usually called after performing a contrast.
#' @param kmeans A list? of positive integers to do k-means.
#' @param out_dir String depicting the output directory to place
#'   tables and plots.
#' @param heatprefix String to prefix heatmap plot filenames
#' @param plot_title String depicting title to embed into plot,
nice_kmeans_heatmap_norm_and_trans <- function(norms, trans,
                                               genes_tocluster = NULL,
                                               sample_columns = NULL,
                                               genes_tohighlight = NULL,
                                               dsqres,
                                               kmeans, out_dir,
                                               heatprefix, plot_title) {

  if (is.null(genes_tocluster)) {
    message(" - Using all genes in normalised matrix for clustering")
    genes_tocluster <- rownames(norms)
  }
  if (is.null(sample_columns)) {
    message(" - Heatmap all samples in normalised matrix for clustering")
    sample_columns <- colnames(norms)
  } else {
    message(" - Heatmap: ", paste0(sample_columns, collapse = ","))
  }

  ## Heatmaps
  message("   - Using normalized counts")
  res_dsqres <- heatmap_with_geneplots(
    norms[genes_tocluster, sample_columns],
    k = kmeans,
    out_dir = out_dir,
    heatprefix = heatprefix,
    prefix_title = paste0(plot_title, " :"),
    highlight_genes = genes_tohighlight)

  dsq_dsq <- left_join(dsqres,
                       res_dsqres$clusters,
                       by = c("gene" = "gene")) %>%
    mutate(norm_cluster = case_when(
             is.na(cluster) ~ "not in norm heatmap",
             TRUE ~ cluster)) %>% select(-cluster)

  ## Heatmaps using Corrected Normalized Counts
  if (!is.null(trans)) {
    message("   - Using transformed counts too")

    ## Heatmaps
    res_dsqres_trans <- heatmap_with_geneplots(
      trans[genes_tocluster, sample_columns],
      k = kmeans,
      out_dir = out_dir,
      heatprefix = paste0(heatprefix, ".vst_corrected"),
      prefix_title = paste0(plot_title, " (vst corrected) :"),
      highlight_genes = genes_tohighlight
    )
    ## Merge Trans cluster
    dsq_dsq <- left_join(dsq_dsq,
                         res_dsqres_trans$clusters,
                         by = c("gene" = "gene")) %>%
      mutate(trans_cluster = case_when(
               is.na(cluster) ~ "not in trans heatmap",
               TRUE ~ cluster)) %>% select(-cluster)
  }

  save_cluster <- file.path(out_dir,
                            paste0("deseq2.results.cluster.k", kmeans, ".tsv"))
  write_tsv(dsq_dsq, save_cluster)
  message("   - Storing Results k", kmeans, ":", save_cluster)
}


#' @title Heatmap with Gene plots
#' @description Generate a clustered heatmaps for a specific k means
#'   value.
#' @param norm_counts A matrix containing normalised values.
#' @param k An integer for the k-value for kmeans.
#' @param out_dir A character sequence depicting the directory
#' @param heatprefix A character sequence depicting the prefix for
#'   heatmaps
#' @param prefix_title A string to prefix the title of heatmap.
#' @param highlight_hmap_genes A vector of genes to highlight in the
#'   heatmap. 
#' @param gene_plot_genes A list of genes to plot. If NULL, use
#'   highlight_hmap_genes. If FALSE, do not plot.
#' @param width_in A positive integer for the number of inches in the plot width
#' @param height_in A positive integer for the number of inches in the plot height.
#' @return A list of two components; clustered tables, and scaled matrix.
heatmap_with_geneplots <- function(norm_counts, k, out_dir = "heatmaps_k",
                                   heatprefix = "heatmap", prefix_title = "",
                                   highlight_hmap_genes = NULL,
                                   gene_plot_genes = NULL,
                                   width_in = 6, height_in = 6) {
  ## normalised counts should be in the correct sample order already  
  options(repr.plot.height = height_in, repr.plot.width = width_in)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  if (is.null(highlight_hmap_genes)) {
    ## If no genes given, show top N
    top_genes <- head(names(sort(rowMeans(norm_counts), decreasing = T)), 30)
    top_title <- paste0(nrow(norm_counts), " DE genes, top ", 
                        length(top_genes), " highlighted")
  } else {
    top_genes <- highlight_hmap_genes
    top_title <- paste0(nrow(norm_counts), " DE genes, ",
                        length(top_genes), " highlighted")
  }
  if (is.null(gene_plot_genes)) {
    gene_plot_genes <- list(topgenes = top_genes)
  }

  scale_mat <- t(scale(t(norm_counts)))
  scale_mat <- scale_mat[complete.cases(scale_mat), ] ## Remove non-zero

  cluster_assignments <- single_kmeans_heatmap(scale_mat, k,
                                               top_genes, prefix_title, top_title,
                                               out_dir, heatprefix, width_in, height_in)
  
  write_tsv(as.data.frame(norm_counts) %>%
            rownames_to_column("gene"), save_norm)
  message("     - Saved Norm: ", save_norm)
  write_tsv(as.data.frame(scale_mat) %>% rownames_to_column("gene"),
            save_scale)
  message("     - Saved Scale: ", save_norm)
  
  if (gene_plot_genes != FALSE) {
    ## Conversion to long table needs to happen here
    long_norm <- norm_counts %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "Sample", values_to = "value")

    long_scale <- scale_mat %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "Sample", values_to = "value")

    gene_clusters_by_score(long_norm, out_dir = file.path(out_dir, "gene_cluster_norm"))
    message("     - Plotting genes Norm: ", save_norm)
    gene_clusters_by_score(long_scale, out_dir = file.path(out_dir, "gene_cluster_scale"))
    message("     - Plotting genes Scaled : ", save_norm)
    
    gene_plots_by_gene(norm_long, scale_long, gene_plot_genes,
                       outprefix = "gene.lists",
                       out_dir = file.path(out_dir, "gene_lists"))
    message("     - Plotting genes list: ", save_norm)
  }
  return(list(clusters = cluster_assignments, scaled = scale_mat))
}


#' @title A single K-means Heatmap
#' @description Plot a clustered heatmap
#' @param scale_mat A matrix of scaled values.
#' @param k A positive integer  for the k value of kmeans
#' @param top_genes A list of genes to plot
#' @param prefix_title A string to prefix the title of heatmap.
#' @param top_title A character string for the title of the heatmap.
#' @param out_dir A character sequence depicting the directory
#' @param heatprefix A character sequence depicting the prefix for
#'   heatmaps
#' @param width_in A positive integer for the number of inches in the plot width
#' @param height_in A positive integer for the number of inches in the plot height.
#' @return A list of two components; clustered tables, and scaled matrix.
single_kmeans_heatmap <- function(scale_mat, k, top_genes,
                                  prefix_title, top_title,
                                  out_dir, heatprefix,
                                  width_in, height_in) {
  ha <- rowAnnotation(
    foo = anno_mark(
      at = which(rownames(scale_mat) %in% top_genes),
      labels = top_genes))

  if (k > 5) {
    rtitle <- "%s"
  } else {
    rtitle <- "Clust %s"
  }

  hm_now <- Heatmap(scale_mat,
                    column_title = paste0(prefix_title, "  ", top_title),
                    row_km = k,
                    cluster_row_slices = FALSE, ## Arrange the K? NO
                    cluster_rows = TRUE, ## Prevents annotation bunching if TRUE
                    show_row_dend = FALSE,
                    row_gap = unit(3, "mm"),
                    name = "mat",
                    row_title = rtitle,
                    col = colorRampPalette(c("black", "#bb0000"))(100),
                    cluster_columns = FALSE,
                    show_row_names = FALSE,
                    right_annotation = ha,
                    column_names_rot = 45,
                    row_names_gp = gpar(fontsize = 4))

  hm_now_drawn <- draw(hm_now)

  ## Get cluster assignments
  cluster_assignments <- (function() {
    cluster_list <- lapply(row_order(hm_now_drawn),
                           function(x) rownames(scale_mat)[x])
    cluster_table <- do.call(
      rbind, lapply(names(cluster_list),
                    function(n) data.frame(gene = cluster_list[[n]],
                                           cluster = n))
    )
    return(cluster_table)
  })()

  ## Unite cluster assignments with Norm, Scaled, and Pvalue
  output_prefix <- file.path(out_dir, paste0(heatprefix, ".k", k, "."))

  save_svg <- paste0(output_prefix, "svg")
  ## save_svgold <- paste0(output_prefix, "_old.svg")
  save_pdf <- paste0(output_prefix, "pdf")
  ## save_png <- paste0(output_prefix, "png")

  save_scale <- paste0(output_prefix, "scale.tsv")
  save_norm <- paste0(output_prefix, "norm.tsv")
  ## save_cluster <- paste0(output_prefix, "clusters.tsv")

  ## Hack to remove the weird grey lines
  ## TODO: Fix, scaling issue...
  ##bhm <- better_pheatmap(hm_now)

  ##svg(save_svgold, width = width_in, height = height_in)
  ##print(hm_now)
  ##dev.off()

  svg(save_svg, width = width_in, height = height_in)
  ##grid.draw(bhm);
  print(hm_now)
  dev.off()

  pdf(save_pdf, width = width_in, height = height_in)
  ##grid.draw(bhm); grid.newpage()
  print(hm_now)
  dev.off()
  return(cluster_assignments)
}



#' @title Better Pheatmaps
#' @description Remove the greylines that plagues these plots by
#'   modifying the raw ggplot grob.
#' @param ph a pheatmap object
#' @return a modified pheatmap object with lines filled with fill colours.
better_pheatmap <- function(ph) {
  lwdd <- 2
  message("-[better_pheatmap]-")
  ## taken from
  ##https://stackoverflow.com/questions/44318690/no-border-color-in-pheatmap-in-r

  if (class(ph) == "Heatmap") {
    ## ComplexHeatmap
    ph <- grid.grabExpr(draw(ph))
    hmaps <- ph$children
    for (hm in names(hmaps)) {
      ## Get all Rects
      ##rects = grep("rect", names(ph$children[[hm]]$children), value=T)
      rects <- names(ph$children[[hm]]$children)

      if (length(rects) > 0) {
        ##rects = names(ph$children[[hm]]$children)
        for (rr in rects) {
          ## Check for rects with tabular data
          hasdim <- dim(ph$children[[hm]]$children[[rr]]$gp$fill)
          if (!is.null(hasdim))
          {
            message("--Here5")
            ph$children[[hm]]$children[[rr]]$gp$col <-
              ph$children[[hm]]$children[[rr]]$gp$fill
            ph$children[[hm]]$children[[rr]]$gp$lwd <- lwdd
          }
        }
      } else {
        if ("gp" %in% names(ph$children[[hm]])) {
          ## Check for rects with tabular data
          hasdim <- dim(ph$children[[hm]]$gp$fill)
          if (!is.null(hasdim))
          {
            message("Setting Cols for: ", hm)
            ph$children[[hm]]$gp$col <- ph$children[[hm]]$gp$fill
            ph$children[[hm]]$gp$lwd <- lwdd
          }
        }
      }
    }
    return(ph)
  } else {
    ## PHEATMAP
    ## Extract the right grob
    grob_classes <- map(ph$gtable$grobs, class)
    idx_grob <- which(
      map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
    grob_names <- names(ph$gtable$grobs[[idx_grob]]$children)
    idx_rect <- grob_names[grep('rect', grob_names)][1]

    ## Remove borders around cells
    ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <-
      ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$fill
    ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 3
    return(ph)
  }
}

#' @title Get a high quality genes list
#' @description Filter a wide-format matrix for high quality genes
#'   (rows) based on detectability thresholds and smallest occurrence.
#' @param sam_mat A numeric matrix. Genes are rows and Samples as
#'   columns. Typically this is the average sample matrix, with
#'   replicate values summarized into mean expression.
#' @param min_occur A positive integer. A threshold (greater than or
#'   equal to) the number of times a gene is detected in a sample,
#'   across all samples (e.g. Gene X is only kept if it appears in 10
#'   samples with a value greater than the `min_detect' threshold.
#'   Default value is 3.
#' @param min_detect A positive integer. A threshold (greater than or
#'   equal to) for the minimum expression a gene can have for it to be
#'   a valid occurrence in a sample. Default values is 10.
#' @param out_dir A string depicting the output directory to place tables.
#' @return A vector of strings depicting a list of high quality genes.
#'   A table of rowSums is written to a file matching
#'   "smallestGroup*-detected*-keep_genes".
#' @export
high_quality_genes <- function(sam_mat,
                               min_occur = 3,
                               min_detect = 10,
                               out_dir) {
  res <- rowSums(sam_mat >= min_detect) >= min_occur
  keep_genes <- res[res == TRUE]
  drop_genes <- res[res == FALSE]

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  write_tsv(
    data.frame(Gene = names(res), Keep = res),
    file.path(out_dir, paste0(
                         "smallestGroup",
                         min_occur, "-detected",
                         min_detect, "-keep_genes"
                       ))
  )
  message(paste0(
    "Dropped: ", length(drop_genes),
    " genes (",
    as.integer(100 * length(drop_genes) / length(res)), "%)"
  ))
  return(keep_genes)
}


#' @title Produce a Volcano Plot
#' @description Generates a volcano plot from the contrast of a DESeq2
#'   analysis with added curves for graphic effect.
#' @param dsqres A table of DESeq2 results
#' @param degenes A vector of genes.
#' @param title A string depicting the title of the plot
#' @param curve A list of three components: sd, sc, offset.
#' @param curve_scale A numeric curve scaling factor.
#' @param curve_show A boolean on whether or not to show the curve.
#' @param ylim A vector of two components to limit the Y-axis.
#' @return A ggplot2 object of the volcano plot.
volcano_plot <- function(dsqres, degenes, title,
                         curve = list(sd = 0.15, sc = 10, offset = 1),
                         curve_scale = 1, curve_show = FALSE, ylim = NULL) {

  options(repr.plot.height = 12, repr.plot.width = 12)

  ## Extract relevant info from DESeq results
  ana <- dsqres %>%
    select(c(gene, log2FoldChange, padj)) %>%
    mutate(mLog10Padj = -log10(padj)) %>%
    arrange(desc(mLog10Padj))

  ## The main function -- do not edit?
  cust_fun <- function(x, sd = 0.15, sc = 10, offset = 1) {
    offset + dnorm(x, mean = 0, sd = sd) * sc
  }

  ## The slope function -- highlight these genes
  cfun <- function(x) {
    return(cust_fun(x,
                    sd = curve$sd * curve_scale,
                    sc = curve$sc * curve_scale,
                    offset = curve$offset))
  }

  ## We highlight genes in the zoomed zone fitting the curve, but the
  ## main DE genes are shown as shapes and highlighted
  red <- ana %>%
    mutate(isTopN = gene %in% degenes) %>%
    mutate(highlight = isTopN | (mLog10Padj > cust_fun(log2FoldChange)))

  max_x <- max(abs(ana$log2FoldChange)) + 0.05 ##symmetry

  plot1 <- red %>%
    ggplot(aes(x = log2FoldChange,
               y = mLog10Padj,
               colour = highlight,
               shape = isTopN,
               label = gene)) +
    geom_point() +
    scale_colour_manual(values = c("TRUE"  = "red", "FALSE" = "grey")) +
    scale_shape_manual(values =  c("TRUE"  =  5, "FALSE" = 19)) +
    scale_x_continuous(limits = c(-max_x,max_x), breaks = waiver(),
                       n.breaks = 10) +
    geom_label_repel(
      data = red %>% filter(highlight == TRUE) %>% head(15),
      box.padding = 0.5,
      max.overlaps = 30,
      colour = "black") +
    ggtitle(title)

  if (curve_show) {
    plot1 <- plot1 + geom_function(fun = cfun, n = 100, colour = "blue")
  }
  if (!is.null(ylim)) {
    plot1 <- plot1 + scale_y_continuous(limits=ylim)
  }
  return(plot1)
}

#' @title Gene clusters by score
#' @description Generate cluster gene plots for various score
#'     thresholds
#' @param tab A dataframe with long plotting data containing at least
#'     the columns `cluster', `condition', `value', and `time'. It can
#'     take normalised or scaled values as input.
#' @param score_thresh Vector of numerics depicting cluster score
#'     thresholds to filter for high quality genes in each cluster
#'     before plotting. Default is \code{c(0, 0.5, 0.9, 0.99)}
#' @param out_dir String representing the directory to store gene
#'     plots. Default value is "gene_cluster"
#' @return A single patchwork plot containing cluster gene plots for
#'     different score thresholds. This plot is also saved to the
#'     `out_dir' directory under the pattern
#'     "gene_plots-k*_montage.svg"
#' @export
gene_clusters_by_score <- function(tab,
                                   score_thresh = c(0, 0.5, 0.9, 0.99),
                                   out_dir = "gene_cluster") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  res <- lapply(score_thresh, function(x) {
    return(cluster_gene_plots(tab, score_thresh = x, out_dir = out_dir))
  })

  newoutprefix <- file.path(
    out_dir,
    paste0("gene_plots-k", num_clusters(tab))
  )

  rplot <- wrap_plots(res, ncol = 2, guides = "collect") &
    theme(plot.title = element_blank(), axis.title = element_blank())

  ggsave(plot = rplot, filename = paste0(newoutprefix, "_montage.svg"),
         dpi = 800, width = 15, height = 15, units = "inches")

  return(rplot)
}

#' @title Cluster Gene Plots
#' @description Generate ribbon and line plots of gene expression data
#'     split by cluster.
#' @param tab A dataframe with long plotting data containing at least
#'     the columns `cluster', `condition', `value', and `time'. It can
#'     take normalised or scaled values as input.
#' @param score_thresh numeric threshold to filter out genes with
#'     correlation scores not matching these values
#' @param out_dir a string depicting the directory of where the plots
#'     are deposited, with the prefix "gene_plots-k".
#' @return ggplot object ribbon and line plots facetted by cluster.
cluster_gene_plots <- function(tab, score_thresh = 0,
                               out_dir = "gene_cluster") {
  tabn <- tab %>% filter(score >= score_thresh)
  ## Keep all clusters even if they're empty
  tabn$cluster <- factor(tabn$cluster, levels = unique(sort(tab$cluster)))
  ## yes
  dat_labs <- generate_labelling_table(tabn)

  time_breaks <- sort(unique(sort(tab$time)))

  ## labelling function for facet
  lab_fun <- function(s) {
    (dat_labs %>% filter(cluster == s))$ftext
  }

  p1 <- tabn %>% ggplot(aes(
                   x = time, y = value, fill = condition,
                   colour = condition, group = condition)) +
    stat_summary(fun.data = "mean_sdl", geom = "ribbon",
                 alpha = 0.1, colour = NA) +
    stat_summary(fun = mean, geom = "line") +
    facet_wrap("cluster", labeller = labeller(cluster = lab_fun),
               drop = FALSE) +
    ylab("Scaled Expression") +
    ggtitle("Gene Trends by cluster",
            subtitle = paste0(
              "Genes (", length(unique(tabn$gene)),
              ") with cluster affinity scores >= ",
              score_thresh
            )) +
    scale_x_continuous(breaks = time_breaks) +
    theme_bw()

  newoutprefix <- file.path(out_dir,
                            paste0("gene_plots-k", num_clusters(tab),
                                   "-score", score_thresh))
  ggsave(plot = p1, filename = paste0(newoutprefix, ".svg"),
         dpi = 800, width = 10, height = 10, units = "inches")
  tabn %>% write_tsv(paste0(newoutprefix, ".tsv"))
  return(p1)
}


#' @title Gene plots by gene
#' @description Generate plots for genes of interest showing their
#'   normalized and z-score scaled expression.
#' @param norm_long A long dataframe of normalized values. Must
#'   contain columns of: gene, condition, time, replicate, value.
#' @param scale_long A long dataframe of scaled values. Same format as
#'   `norm_long'.
#' @param gois_list A list of vectors. A list of genes which are
#'   referenced by specific names. Unique plots will be generated for
#'   list. e.g.
#'       list(mygeneset1 = c("Msgn1", "Osr1", "Rspo3", "Fgf8", "Wnt3a"),
#'            mygeneset2 = c("Mesp1", "Foxa2", "Sox17", "Lhx1", "Cer1"))
#' @param outprefix A string depicting the output prefix for the
#'   plots. These will be appended with the title of each gene list
#'   and for the normalized and scaled.
#' @param out_dir A string denoting the output directory to store
#'   plots. Default is is "gene_lists".
#' @export
gene_plots_by_gene <- function(norm_long, scale_long, gois_list,
                               outprefix = "gene.lists",
                               out_dir = "gene_lists") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  plot_dims <- function(n) {
    ## Calculates plot and width height for a given N plots
    w <- ceiling(sqrt(n))
    h <- ceiling(n / w)
    return(list(w = w, h = h))
  }

  time_breaks <- sort(unique(sort(norm_long$time)))

  pgene_list <- lapply(names(gois_list), function(glist_name) {
    glist <- gois_list[[glist_name]]
    genes_found <- unique((norm_long %>% filter(gene %in% glist))$gene)

    if (length(genes_found) < 1) {
      message("no genes found for: ", glist_name)
      return(NULL)
    }

    pgene_norm <- norm_long %>%
      filter(gene %in% genes_found) %>%
      ggplot(aes(
        x = time, y = value,
        fill = condition, colour = condition, group = condition
      )) +
      stat_summary(
        fun.data = "mean_sdl", geom = "ribbon",
        alpha = 0.1, colour = NA
      ) +
      stat_summary(fun = mean, geom = "line") +
      geom_jitter(width = 0.1) +
      scale_y_log10() +
      facet_wrap("gene", scales = "free_y") +
      ylab("Log10 Normalised Expression") +
      ggtitle(paste0(
        "Normalised Expression Time Plots ",
        "grouped by Genes of Interest"
      )) +
      scale_x_continuous(breaks = time_breaks) +
      theme_bw()

    pdims <- plot_dims(length(genes_found))
    ggsave(
      plot = pgene_norm,
      filename = file.path(
        out_dir,
        paste0(
          outprefix, ".",
          glist_name, ".normalised.svg"
        )
      ),
      dpi = 800, width = pdims$w * 2, height = pdims$h * 1.5,
      units = "inches"
    )

    pgene_scale <- scale_long %>%
      filter(gene %in% genes_found) %>%
      ggplot(aes(
        x = time, y = value,
        fill = condition, colour = condition, group = condition
      )) +
      stat_summary(
        fun.data = "mean_sdl", geom = "ribbon",
        alpha = 0.1, colour = NA
      ) +
      stat_summary(fun = mean, geom = "line") +
      geom_jitter(width = 0.1) +
      facet_wrap("gene", scales = "free_y") +
      ylab("Scaled Expression") +
      ggtitle("Scaled Expression Time Plots, grouped by Genes of Interest") +
      scale_x_continuous(breaks = time_breaks) +
      theme_bw()

    ggsave(
      plot = pgene_scale,
      filename = file.path(out_dir, paste0(
                                      outprefix, ".",
                                      glist_name, ".scaled.svg"
                                    )),
      dpi = 800, width = pdims$w * 2, height = pdims$h * 1.5, units = "inches"
    )

    return(list(norm = pgene_norm, scale = pgene_scale))
  })
}

#' @title Generate a label table
#' @description Counts how many genes in each cluster for use as
#'   plotting labels in `cluster_gene_plots'
#' @param tabn a gene table of long data with a `cluster' field.
#' @return a table with columns of `cluster', `n' for how many genes
#'   in that cluster, and `ftext' for labelling.
generate_labelling_table <- function(tabn) {
  dat_labs <- tabn %>% group_by(cluster) %>%
    summarise(n = n_distinct(gene)) %>%
    mutate(ftext = paste0("cluster ", cluster, ", ", n, " genes"))

  ## Sanity check that all factors are in summary
  miss <- which(!(levels(tabn$cluster) %in% dat_labs$cluster))
  if (length(miss) > 0) {
    temp <- cbind(cluster = miss, n = 0, ftext = paste0("cluster ",
                                                        miss, ", 0 genes"))
    dat_labs <- rbind(dat_labs, temp)
  }
  dat_labs <- dat_labs %>% arrange(cluster)
  return(dat_labs)
}

#' @title Number of clusters in tabs
#' @description Count the number of unique values in the `cluster'
#'   column of a table
#' @param tab Dataframe. Must contain a `cluster' vector
#' @return Positive integer.
#' @examples
#' RNASeqHelper:::num_clusters(
#'     data.frame(cluster = c(1, 1, 2, 3, 2))
#' )
num_clusters <- function(tab) {
  return(length(unique(tab$cluster)))
}


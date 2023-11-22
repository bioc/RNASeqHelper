
#' @importFrom tibble rownames_to_column
#' @importFrom readr read_tsv write_tsv
#' @importFrom dplyr arrange case_when filter group_by left_join mutate desc
#' @importFrom dplyr n_distinct select summarise %>%
#' @importFrom ggplot2 ggplot facet_wrap stat_summary ylab xlab
#' @importFrom ggplot2 scale_x_continuous ggtitle ggsave theme theme_bw aes
#' @importFrom ggplot2 scale_y_continuous scale_y_log10 geom_point waiver element_blank
#' @importFrom ggplot2 scale_colour_manual scale_shape_manual geom_function units
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
  res_dsqres <- single_kmeans_heatmap(
    norms[genes_tocluster, sample_columns],
    k = kmeans,
    out_dir = out_dir,
    heatprefix = heatprefix,
    prefix_title = paste0(plot_title, " :"),
    highlight_genes = genes_tohighlight
  )

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
    res_dsqres_trans <- single_kmeans_heatmap(
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


single_kmeans_heatmap <- function(norm_counts, k, out_dir = "heatmaps_k",
                                heatprefix = "heatmap", prefix_title = "",
                                highlight_genes = NULL,
                                width_in = 6, height_in = 6) {
  ## normalised counts should be in the correct sample order already
  options(repr.plot.height = height_in, repr.plot.width = width_in)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  if (is.null(highlight_genes)) {
    ## If no genes given, show top N
    top_genes <- head(names(sort(rowMeans(norm_counts),
                                 decreasing = T)), 30)
    top_title <- paste0(nrow(norm_counts), " DE genes, top ",
                        length(top_genes), " highlighted")
  } else {
    top_genes <- highlight_genes
    top_title <- paste0(nrow(norm_counts), " DE genes, ",
                        length(top_genes), " highlighted")
  }

  scaledata <- t(scale(t(norm_counts)))
  scaledata <- scaledata[complete.cases(scaledata), ] ## Remove non-zero

  ha <- rowAnnotation(
    foo = anno_mark(
      at = which(rownames(scaledata) %in% top_genes),
      labels = top_genes))

  if (k > 5) {
    rtitle <- "%s"
  } else {
    rtitle <- "Clust %s"
  }

  hm_now <- Heatmap(scaledata,
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
                           function(x) rownames(scaledata)[x])
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
  ###bhm <- better_pheatmap(hm_now)

  ##svg(save_svgold, width = width_in, height = height_in)
  ##print(hm_now)
  ##dev.off()

  svg(save_svg, width = width_in, height = height_in)
  ##grid.draw(bhm)
  print(hm_now)
  dev.off()

  pdf(save_pdf, width = width_in, height = height_in)
  ##grid.draw(bhm)
  ##grid.newpage()
  print(hm_now)
  dev.off()

  write_tsv(as.data.frame(norm_counts) %>%
            rownames_to_column("gene"), save_norm)
  message("     - Saved Norm: ", save_norm)

  write_tsv(as.data.frame(scaledata) %>% rownames_to_column("gene"),
            save_scale)
  message("     - Saved Scale: ", save_norm)

  return(list(plot = hm_now,
              clusters = cluster_assignments, scaled = scaledata))
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
      as.data.frame(res),
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

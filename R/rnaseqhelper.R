## Mehmet Tekman, 2023

#' @importFrom tidyverse group_by arrange summarise mutate select write_tsv read_tsv summarize n_distinct left_join case_when rownames_to_column write_tsv %>% 
#' @importFrom ggplot2 ggplot facet_wrap stat_summary ylab xlab scale_x_continuous ggtitle ggsave theme theme_bw
#' @importFrom patchwork wrap_plots
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq plotPCA vst
#' @importFrom pheatmap pheatmap

library(tidyverse)
library(RColorBrewer)
library(patchwork)

#' @title Generate a label table
#' @description Counts how many genes in each cluster for use as
#'   plotting labels in `cluster_gene_plots'
#' @param tabn a gene table of long data with a `cluster' field.
#' @return a table with columns of `cluster', `n' for how many genes
#'   in that cluster, and `ftext' for labelling.
#' @examples
#' TODO
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

#' @title Cluster Gene Plots
#' @description Generate ribbon and line plots of gene expression data
#'   split by cluster.
#' @param tab A dataframe with long plotting data containing at least
#'   the columns `cluster', `condition', `value', and `time'
#' @param score_thresh numeric threshold to filter out genes with
#'   correlation scores not matching these values
#' @param dirout a string depicting the directory of where the plots
#'   are deposited, with the prefix "gene_plots-k".
#' @return ggplot object ribbon and line plots facetted by cluster.
cluster_gene_plots <- function(tab, score_thresh = 0, dirout = "gene_cluster") {
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
  
  newoutprefix <- file.path(dirout,
                            paste0("gene_plots-k", num_clusters(tab),
                                   "-score", score_thresh))
  ggplot2::ggsave(
             plot = p1, filename = paste0(newoutprefix, ".svg"),
             dpi = 800, width = 10, height = 10, unit = "in")
  tabn %>% write_tsv(paste0(newoutprefix, ".tsv"))
  return(p1)
}

#' @title Number of clusters in tabs
#' @description Count the number of unique values in the `cluster'
#'   column of a table
#' @param tab Dataframe. Must contain a `cluster' vector
#' @return Positive integer.
#' @examples
#' RNAseqHelper:::num_clusters(
#'     data.frame(cluster = c(1, 1, 2, 3, 2))
#' )
num_clusters <- function(tab) {
  return(length(unique(tab$cluster)))
}

#' @title Gene clusters by score
#' @description Generate cluster gene plots for various score
#'   thresholds
#' @param tab Dataframe containing TODO
#' @param score_thresh Vector of numerics depicting cluster score
#'   thresholds to filter for high quality genes in each cluster
#'   before plotting. Default is \code{c(0, 0.5, 0.9, 0.99)}
#' @param out_dir String representing the directory to store gene
#'   plots. Default value is "gene_cluster"
#' @return A single patchwork plot containing cluster gene plots for
#'   different score thresholds. This plot is also saved to the
#'   `out_dir' directory under the pattern "gene_plots-k*_montage.svg"
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

  ggplot2::ggsave(
             plot = rplot, filename = paste0(newoutprefix, "_montage.svg"),
             dpi = 800, width = 15, height = 15, unit = "in"
           )

  return(rplot)
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
#' @return A vector of strings depicting a list of high quality genes.
#'   A table of rowSums is written to a file matching
#'   "smallestGroup*-detected*-keep_genes".
high_quality_genes <- function(sam_mat,
                               min_occur = 3,
                               min_detect = 10) {
  res <- rowSums(sam_mat >= min_detect) >= min_occur
  keep_genes <- res[res == TRUE]
  drop_genes <- res[res == FALSE]
  write_tsv(
    as.data.frame(res),
    paste0(
      "smallestGroup",
      min_occur, "-detected",
      min_detect, "-keep_genes"
    )
  )
  message(paste0(
    "Dropped: ", length(drop_genes),
    " genes (",
    as.integer(100 * length(drop_genes) / length(res)), "%)"
  ))
  return(keep_genes)
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
            dpi = 800, width = pdims$w * 2, height = pdims$h * 1.5, unit = "in"
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
            dpi = 800, width = pdims$w * 2, height = pdims$h * 1.5, unit = "in"
        )

        return(list(norm = pgene_norm, scale = pgene_scale))
    })
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
    grob_classes <- purrr::map(ph$gtable$grobs, class)
    idx_grob <- which(
      purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
    grob_names <- names(ph$gtable$grobs[[idx_grob]]$children)
    idx_rect <- grob_names[grep('rect', grob_names)][1]

    ## Remove borders around cells
    ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <-
      ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$fill
    ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 3
    return(ph)
  }
}

#' @title Produce a single heatmap
#' @description Generate a single pheatmap and save to file
#' @param scaledata_k A matrix of scaled gene data, rows as genes and
#'   columns as samples.
#' @param ordered_cols A vectors of sample names.
#' @param gaps_idxs vector of row indices that show where to put gaps
#'   into heatmap.
#' @param main A string depicting the title of the heatmap.
#' @param show_genes A vectors of gene names to highlight in the plot.
#' @param out_dir A string depicting the output directory of the
#'   heatmap.
#' @param outfix A string depicting the prefix of the heatmap filename
#' @return None. A heatmap and a corrected "better_pheatmap" is saved
#'   to the output directory.
single_heatmap <- function(scaledata_k, ordered_cols, gaps_idxs, main,
                           show_genes, out_dir = "", outfix = "") {
  ## TODO, feed in the output_svg_prefix
  
  output_svg_prefix <- file.path(out_dir, paste0("genes.k", k, "."))

  colors_kr <- colorRampPalette(c("black", "#bb0000"))(100)
  svg(tempfile(fileext = ".svg"))
  ph <- pheatmap::pheatmap(
                    scaledata_k[, ordered_cols],
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    show_rownames = TRUE,
                    labels_row = show_genes,
                    cellwidth = 40,
                    col = colors_kr,
                    fontsize_row = 8,
                    border_color = NA,
                    gaps_row = gaps_idxs, # gap after each block
                    main = main,
                    ...
                  )
  dev.off()

  svg(paste0(output_svg_prefix, outfix, ".svg"))
  graphics::plot.new()
  ## print(ph)
  print(better_pheatmap(ph))
  dev.off()
}


#' @title Correlation Gene Cluster Scores
#' @description Calculate correlation scores of genes to the centroid
#'   of their clusters.
#' @param scaledata a matrix of numeric scaled data.
#' @param kmolten a dataframe of k cluster centroids.
#' @param k_clusters a vector of integer clusters, associated to the
#'   row (gene) names of the `scaledata'.
#' @param condition_list a vector of sample names.
#' @return a list of gene scores, centroids, sorted gene rows.
run_kmolten_perclust <- function(scaledata, kmolten,
                                 k_clusters,
                                 condition_list) {
  clusters_unique <- sort(unique(kmolten$cluster))
  cores <- list()
  kmolten_list <- list()
  scores <- list()

  for (i in clusters_unique) {
    core_i <- kmolten[kmolten$cluster == i, ]
    core_i$sample <- factor(core_i$sample, levels = condition_list)
    k_i <- scaledata[k_clusters == i, ]
    corescore_i <- function(x) cor(x, core_i$value)
    score_i <- apply(k_i, 1, corescore_i)
    k_i_molten <- as.data.frame(as.table(k_i))
    colnames(k_i_molten) <- c("gene", "sample", "value")
    k_i_molten <- merge(k_i_molten, score_i, by.x = "gene",
                        by.y = "row.names", all.x = T)
    colnames(k_i_molten) <- c('gene', 'sample', 'value', 'score')
    k_i_molten <- k_i_molten[order(k_i_molten$value), ]
    kmolten_list[[i]] <- k_i_molten
    scores[[i]] <- score_i
    cores[[i]] <- core_i
  }
  ## name collected values
  names(kmolten_list) <- paste0("K", seq_along(clusters_unique), "molten")
  names(cores) <- paste0("core", seq_along(clusters_unique))
  names(scores) <- paste0("score", seq_along(clusters_unique))

  return(list(scores = scores, cores = cores, kmolten_list = kmolten_list))
}


#' @title Produce a K-means clustered Heatmap
#' @description Generates a k-means clustered heatmap of z-score
#'   scaled gene expression values, preserving the sample order.
#' @param data A numeric matrix. Must be normalised but not already
#'   scaled.
#' @param k A positive integer. The number of desired partitions of
#'   the data.
#' @param plot_list A named list of vectors. Each name refers to a
#'   preferred subset/ordering of the samples in the numeric `data
#'   matrix. e.g. list(order1 = c("Sample1A",
#'   "Sample1B, "Sample2A", ...), order2 = c("Sample1A", "Sample2A", Sample1B",
#'   ...))
#' @param out_dir A string depicting the output directory of the
#'   heatmap. Default is "heatmaps_k".
#' @param highlight_genes A vector of genes to highlight in the output
#'   plots. Default is NULL, highlighting none.
#' @param ... All other arguments are fed into `pheatmap::pheatmap'
#' @return A list of numeric matrices, partitioned by K value
#'   contained z-score scaled genes values.
do_kmeans <- function(data, k, plot_list, out_dir = "heatmaps_k",
                      highlight_genes = NULL, ...) {
  out_dir <- paste0(out_dir, k)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  scaledata <- t(scale(t(data)))
  scaledata <- scaledata[complete.cases(scaledata), ] ## Remove non-zero

  k_clust <- kmeans(scaledata, centers = k, nstart = 1000, iter.max = 30)
  k_clusters <- k_clust$cluster

  clust_centroid <- function(i, dat, clusters) {
    ind <- (clusters == i)
    colMeans(dat[ind, ])
  }
  ## is a matrix
  k_clustcentroids <- sapply(levels(factor(k_clusters)),
                             clust_centroid, scaledata, k_clusters)

  kmolten <- as.data.frame(as.table(k_clustcentroids))
  colnames(kmolten) <- c("sample", "cluster", "value")
  ## ensure correct factorizing
  kmolten$sample <- factor(kmolten$sample)

  res_rkmopcl <- run_kmolten_perclust(scaledata, kmolten, k_clusters,
                                      plot_list[[1]])
  scores <- res_rkmopcl$scores
  ##cores <- res_rkmopcl$cores
  ##kmolten_list <- res_rkmopcl$kmolten_list

  scaledata_k <-cbind(scaledata,
                      cluster = k_clust$cluster,
                      score = Reduce(`c`, scores)[rownames(scaledata)])
  scaledata_k <- scaledata_k[order(scaledata_k[, "cluster"],
                                   scaledata_k[, "score"]), ]

  select_df <- function(df, val, col_selector) {
    df[df[, col_selector] == val, ]
  }
  df2dflist <- function(df, col_selector) { # actually it is split()
    message("-[df2dflist]-")
    col_vals <- unique(df[, col_selector])
    dfl <- lapply(seq(col_vals), function(i) {
      select_df(df,
                val = col_vals[i],
                col_selector
                )
    })
    names(dfl) <- col_vals
    dfl
  }
  scaledata_list <- df2dflist(scaledata_k, "cluster")
  names(scaledata_list) <- paste("cluster", 1:k, sep = "")

  ## scaledata_k is one matrix! only column 'cluster'
  gaps_idxs <- cumsum(table(scaledata_k[, "cluster"])) 
  show_genes <- rownames(scaledata)
  if (!is.null(highlight_genes)) {
    ## Blank out any non-listed genes
    show_genes[!(show_genes %in% highlight_genes)] <- ""
  }
  main_title <- paste0(nrow(scaledata), " genes, clustered into k=", 
                 k, " ordered by ", outfix)

  for (name in names(plot_list)) {
    message("-[single_heatmap]- ", name)
    single_heatmap(scaledata, plot_list[[name]], gaps_idxs,
                   main_title, show_genes,
                   name)
  }
  return(scaledata_k)
}

#' @title Set colours to fill values in plots
#' @description An experimental function to remove the horrible gray
#'   lines produced in some SVG/PDF renderers due to the colours of
#'   rectangles having no assigned color or size.
#' @param rh A ggplot grob-like object.
#' @return None. It modifies the object in place.
get_filld <- function(rh) {
  getatts <- function(x) attributes(x)$names

  ## If it has fill and dimension, we modify it
  if ("gp" %in% getatts(rh)) {
    if ("fill" %in% getatts(rh$gp)) {
      hasdim = dim(rh$gp$fill)
      if (!is.null(hasdim)) {
        message("Changed: ", rh, " ", hasdim)
        rh$gp$col <- rh$gp$fill
        rh$gp$lwd <- 3
      }
    }
  }

  if ("children" %in% getatts(rh)) {
    for (child in names(rh$children)) {
      rh$children[[child]] <- get_filld(rh$children[[child]])
    }
  }
}


#' @title Run DESeq with sensible defaults
#' @description Runs DESeq with filtering thresholds and shows PCA for
#'   variance-stabilize-transformed data.
#' @param tab a matrix of samples (columns) and genes (rows)
#' @param phenotype_data a table with samples (rows) with extra
#'   columns for annotation groups. One of these groups must be
#'   Condition and will be used to normalize the data.
#' @param min_occur A positive integer. A threshold (greater than or
#'   equal to) the number of times a gene is detected in a sample,
#'   across all samples (e.g. Gene X is only kept if it appears in 10
#'   samples with a value greater than the `min_detect' threshold.
#'   Default value is 3.
#' @param min_detect A positive integer. A threshold (greater than or
#'   equal to) for the minimum expression a gene can have for it to be
#'   a valid occurrence in a sample. Default values is 10.
#' @return list of tables
run_deseq <- function(tab, phenotype_data, min_detect = 10, min_occur = 3) {
  sub_as <- (tab[rowSums(tab > min_detect) > min_samples,])
  print(dim(sub_as))

  ddsObj <- DESeqDataSetFromMatrix(
    countData = sub_as,
    colData = phenotype_data,
    design = ~Condition
  )
  ddsObj <- DESeq(ddsObj)

  vsd1 <- vst(ddsObj, blind=FALSE)
  p1 <- plotPCA(vsd1, intgroup=c("Condition"),returnData =F )
  vsd2 <- vst(ddsObj, blind=TRUE)
  p2 <- plotPCA(vsd2, intgroup=c("Condition"),returnData =F )
  return(list(tab = tab, phenotype = phenotype_data,
              bFalse = p1, bTrue = p2, 
              vFalse = vsd1, vTrue = vsd2,
              ddsObj=ddsObj, sub=sub_as))
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
    scale_x_continuous(lim = c(-max_x,max_x), breaks = waiver(),
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
    plot1 <- plot1 + scale_y_continuous(lim=ylim)
  }
  return(plot1)
}

#' @title Pairwise heatmaps and volcano plots
#' @description For a given contrast and DESeq object produce volcano
#'   plots of the DE genes, and generate a global heatmap with those
#'   DE genes and a pairwise one with just the numerator and
#'   denominator samples.
#' @param ddsObj A DESeq object,
#' @param transformed_counts
#' @param numer A string prefix to select sample columns that will be
#'   part of the numerator part of the contrast.
#' @param denom A string prefix to select sample columns that will
#'   part of the denominator part of the contrast.
#' @param top_ngenes_tocluster A positive integer. How many of the top
#'   log10 padj values to take. Default is 2000.
#' @param top_ngenes_tohighlight A positive integer. How many of the
#'   top log10 padj values to label in the plots. Default is 50.
#' @param lFC_zoom A number to depict the Log2FC threshold of the
#'   zoomed plot.
#' @param pAdj_zoom A number to depict the adjusted p-value threshold of the
#'   zoomed plot.
#' @param kmeans A list? of positive integers to do k-means
#' @param out_dirprefix
#' @return
pairwise_hmap_volcano <- function(ddsObj,
                                  transformed_counts = NULL,
                                  numer = "FRT", denom = "TAF2",
                                  top_ngenes_tocluster = 2000,
                                  top_ngenes_tohighlight = 50,
                                  lFC_zoom = 1.5, pAdj_zoom = 20,
                                  kmeans = 2, out_dirprefix = ".") {
  ntitle <- paste0(numer, " vs ", denom)
  outdir <- file.path(out_dirprefix,
                      tolower(gsub("[^A-Za-z0-9]", "_", ntitle)))

  dir.create(outdir, recursive=T, showWarnings=FALSE)

  ## 1. Perform RNA-seq contrast between numerator and denominator
  dsqres <- DESeq2::results(
                      ddsObj,
                      contrast = c("Condition", numer, denom),
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
                     ntitle,
                     curve = list(sd = 0.3, sc = 60, offset = 10),
                     curve_show = F)
  ## Volcano Plots zoomed in
  p2 <- volcano_plot(dsqres %>% filter(abs(log2FoldChange) < lFC_zoom &
                                       mLog10Padj < pAdj_zoom),
                     top_genes_tohighlight,
                     paste0(ntitle, " (zoomed)"),
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
    message("Calculating k=", kmk)

    ## Pairwise Heatmap of Normalised and Transformed Counts
    nice_kmeans_heatmap_norm_and_trans(
      norm_counts, transformed_counts,
      genes_tocluster = top_genes_tocluster,
      sample_columns = sample_columns,
      genes_tohighlight = top_genes_tohighlight,
      dsqres, kmk,
      out_dir = file.path(outdir, paste0("kmeans", kmk)),
      heatprefix = "heatmap_pairwise",
      ntitle = "Heatmap Pairwise")

    ## Global Heatmap of Normalised and Transformed Counts
    nice_kmeans_heatmap_norm_and_trans(
      norm_counts, transformed_counts,
      genes_tocluster = top_genes_tocluster,
      sample_columns = NULL,
      genes_tohighlight = top_genes_tohighlight,
      dsqres, kmk,
      out_dir = file.path(outdir, paste0("kmeans", kmk)),
      heatprefix = "heatmap_all",
      ntitle = "Heatmap All")
  }
}


nice_kmeans_heatmap_norm_and_trans <- function(norms, trans,
                                               genes_tocluster = NULL, 
                                               sample_columns = NULL, 
                                               genes_tohighlight = NULL, 
                                               dsqres, 
                                               kmeans, out_dir, 
                                               heatprefix, ntitle) {

  if (is.null(genes_tocluster)) {
    message("using all genes in normalised matrix for clustering")
    genes_tocluster = rownames(norms)
  }
  if (is.null(sample_columns)) {
    message("using all samples in normalised matrix for clustering")
    sample_columns = colnames(norms)
  }
  
  ## Heatmaps
  res_dsqres <- nice_kmeans_heatmap(
    norms[genes_tocluster, sample_columns],
    k = kmeans,
    out_dir = out_dir,
    heatprefix = heatprefix,
    prefix_title = paste0(ntitle, " :"),
    highlight_genes = genes_tohighlight
  )
  ## Merge Norm cluster
  ##message("SAVING")
  ##saveRDS(list(dsqres, res_dsqres), 
  ##        paste0(format(Sys.time(), "%H-%M-%S--"), "BUMM.rds"))
  
  dsq_dsq <- left_join(dsqres, 
                       res_dsqres$clusters,
                       by = c("gene" = "gene")) %>% 
    mutate(norm_cluster = case_when(
             is.na(cluster) ~ "not in norm heatmap",
             TRUE ~ cluster)) %>% select(-cluster)
  
  ## Heatmaps using Corrected Normalized Counts
  if (!is.null(trans)) {
    message("Using transformed counts too")

    ## Heatmaps
    res_dsqres_trans <- nice_kmeans_heatmap(
      trans[genes_tocluster, sample_columns], 
      k = kmeans, 
      out_dir = out_dir,
      heatprefix = paste0(heatprefix, ".vst_corrected"),
      prefix_title = paste0(ntitle, " (vst corrected) :"),
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
  message("Saved DESeq k", kmeans, "fill:rgb(33.72549%,0%,0%);: ", save_cluster)
}


nice_kmeans_heatmap <- function(norm_counts, k, out_dir = "heatmaps_k",
                                heatprefix = "heatmap", prefix_title = "",
                                highlight_genes = NULL,
                                width_in = 6, height_in = 6) {
  ## normalised counts should be in the correct sample order already
  suppressPackageStartupMessages(library(ComplexHeatmap))
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
  message("Saved Norm: ", save_norm)

  write_tsv(as.data.frame(scaledata) %>% rownames_to_column("gene"),
            save_scale)
  message("Saved Scale: ", save_norm)

  return(list(plot = hm_now,
              clusters = cluster_assignments, scaled = scaledata))
}

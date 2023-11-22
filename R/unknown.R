                                        # UNUSED
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
#' @param ... All other arguments are fed into `pheatmap'
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
  for (name in names(plot_list)) {
    message("-[single_heatmap]- ", name)
    outfix <- plot_list[[name]]
    main_title <- paste0(nrow(scaledata), " genes, clustered into k=",
                         k, " ordered by ", outfix)

    single_heatmap(scaledata, outfix, gaps_idxs,
                   main_title, show_genes,
                   name)
  }
  return(scaledata_k)
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
  ph <- pheatmap(
    scaledata_k[, ordered_cols],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    labels_row = show_genes,
    cellwidth = 40,
    color = colors_kr,
    fontsize_row = 8,
    border_color = NA,
    gaps_row = gaps_idxs, # gap after each block
    main = main,
    ...
  )
  dev.off()

  svg(paste0(output_svg_prefix, outfix, ".svg"))
  plot.new()
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

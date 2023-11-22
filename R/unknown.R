
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


#' @title Cluster Gene Plots
#' @description Generate ribbon and line plots of gene expression data
#'   split by cluster.
#' @param tab A dataframe with long plotting data containing at least
#'   the columns `cluster', `condition', `value', and `time'
#' @param score_thresh numeric threshold to filter out genes with
#'   correlation scores not matching these values
#' @param out_dir a string depicting the directory of where the plots
#'   are deposited, with the prefix "gene_plots-k".
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

  ggsave(plot = rplot, filename = paste0(newoutprefix, "_montage.svg"),
         dpi = 800, width = 15, height = 15, units = "inches")

  return(rplot)
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

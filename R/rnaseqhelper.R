#' @importFrom dplyr arrange case_when filter group_by left_join
#'     mutate desc n_distinct select summarise %>% .data
#' @importFrom ggplot2 ggplot facet_wrap stat_summary ylab xlab theme
#'     theme_bw scale_x_continuous ggtitle ggsave aes_string
#'     geom_point waiver scale_y_continuous scale_y_log10
#'     scale_colour_manual element_blank scale_shape_manual
#'     geom_function labeller geom_jitter
#' @importFrom ggrepel geom_label_repel
#' @importFrom graphics plot.new
#' @importFrom grid grid.grabExpr gpar unit
#' @importFrom grDevices colorRampPalette dev.off pdf svg
#' @importFrom patchwork wrap_plots
#' @importFrom pheatmap pheatmap
#' @importFrom purrr map map_lgl
#' @importFrom readr read_tsv write_tsv
#' @importFrom stats complete.cases cor dnorm kmeans time
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer all_of
#' @importFrom utils head modifyList
#' @importFrom ComplexHeatmap Heatmap anno_mark draw row_order
#'     rowAnnotation
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq plotPCA vst results
#'     counts
#' @importFrom SummarizedExperiment assay colData

# MAIN #

#' @title Run a full RNASeqhelper analysis
#' @description Single function to call all other functions in the
#'     analysis and perform a full analysis from a few starting
#'     parameters
#' @param tab a matrix of samples (columns) and genes (rows)
#' @param phenotype_data a table with samples (rows) with extra
#'     columns for annotation groups. One of these groups must be
#'     condition and will be used to normalize the data.
#' @param out_dir String depicting the output directory to place
#'     tables and plots.
#' @param numer A string value to select sample columns from the
#'     \code{condition} column that will be part of the numerator part
#'     of the contrast.
#' @param denom A string value to select sample columns from the
#'     \code{condition} column, that will be the denominator part of
#'     the contrast.
#' @param keep_params A list of two components \code{min_occur} and
#'     \code{min_detect} as defined by the default
#'     \code{high_quality_genes} function.
#' @param heat_params A list of parameters to override the default
#'     \code{pairwise_heatmap_volcano} function. If \code{NULL}, then
#'     use the default.
#' @param volcano_params A list of parameters to override the default
#'     \code{do_volcanos} function. If \code{NULL}, then use the
#'     default.
#' @return None.
#' @examples
#' n <- 1000
#' tab <- matrix(as.integer(rnorm(n**2, 1000, 500)), ncol = n)
#' tab <- tab - min(tab)
#' rownames(tab) <- paste0("G", 1:n)
#' colnames(tab) <- paste0("S", 1:n)
#' phenotype_data <- data.frame(
#'     sample = paste0("S", 1:n),
#'     condition = c(rep("red", n / 2), rep("green", n / 2)),
#'     time = as.integer(rnorm(n, 2, 0.5) + 1) * 5
#' )
#' rnaseqhelper(tab, phenotype_data,
#'     out_dir = "/tmp", "green", "red",
#'     heat_params = list(
#'         score_thresh = c(0.2, 0.5),
#'         kmeans_list = 2
#'     )
#' )
#' @export
rnaseqhelper <- function(tab, phenotype_data, out_dir = "/tmp",
                        numer, denom,
                        keep_params = list(),
                        heat_params = list(), volcano_params = list()) {
    ## out_dir="test/1_genes"
    keep_params_defaults <- as.list(formals(high_quality_genes))
    keep_params <- modifyList(keep_params_defaults, keep_params)
    keep_params$sam_mat <- tab
    keep_params$out_dir <- file.path(out_dir, "0_genes")

    keep_genes <- do.call(high_quality_genes, keep_params)

    res <- run_deseq(tab, keep_genes, phenotype_data,
        out_dir = file.path(out_dir, "1_matrices_and_deseq")
    )
    ## res$ddsObj
    ## out_dirprefix="test/outputs"
    heat_params_defaults <- as.list(formals(pairwise_hmap_volcano))
    heat_params <- modifyList(heat_params_defaults, heat_params)
    heat_params$ddsObj <- res$ddsObj
    heat_params$transformed_counts <- assay(res$vFalse)
    heat_params$out_dirprefix <- file.path(out_dir, "2_heatmaps")
    heat_params$numer <- numer
    heat_params$denom <- denom
    heat_params$phv <- do.call(pairwise_hmap_volcano, heat_params)
}

#' @title Run DESeq with sensible defaults
#' @description Runs DESeq with filtering thresholds and shows PCA for
#'     variance-stabilize-transformed data.
#' @param tab a matrix of samples (columns) and informative genes
#'     (rows), ideally subset using the function
#'     \code{high_quality_genes}
#' @param keep_genes a vector of genes to subset.
#' @param phenotype_data a table with samples (rows) with extra
#'     columns for annotation groups. One of these groups must be
#'     condition and will be used to normalize the data.
#' @param out_dir String depicting the output directory to place
#'     tables and plots.
#' @return list of tables
#' @examples
#' n <- 1000
#' tab <- matrix(as.integer(rnorm(n**2, 1000, 500)), ncol = n)
#' tab <- tab - min(tab)
#' rownames(tab) <- paste0("G", 1:n)
#' colnames(tab) <- paste0("S", 1:n)
#' phenotype_data <- data.frame(
#'     sample = paste0("S", 1:n),
#'     condition = c(rep("red", n / 2), rep("green", n / 2)),
#'     time = as.integer(rnorm(n, 2, 0.5) + 1) * 5
#' )
#' keep_genes <- paste0("G", 1:n)
#' res <- run_deseq(tab, keep_genes, phenotype_data, "/tmp")
#' @export
run_deseq <- function(tab, keep_genes, phenotype_data, out_dir) {

    if (!("condition" %in% colnames(phenotype_data))) {
        stop("Could not find `condition' column in phenotype data")
    }
    sub_as <- tab[keep_genes, ]

    ddsObj <- DESeqDataSetFromMatrix(
        countData = sub_as,
        colData = phenotype_data,
        design = ~condition
    )
    ddsObj <- DESeq(ddsObj)

    vsd1 <- vst(ddsObj, blind = FALSE)
    p1 <- plotPCA(vsd1, intgroup = c("condition"),
                returnData = FALSE)
    vsd2 <- vst(ddsObj, blind = TRUE)
    p2 <- plotPCA(vsd2, intgroup = c("condition"),
                returnData = FALSE)
    res <- list(
        tab = tab, phenotype = phenotype_data,
        bFalse = p1, bTrue = p2,
        vFalse = vsd1, vTrue = vsd2,
        ddsObj = ddsObj, sub = sub_as
    )    
    pca_and_matrices(res, out_dir)    
    return(res)
}

#' @title Print and Store PCA and Matrices
#' @description Takes the output of a DESeq2 run and generates
#'     normalization tables and PCAs
#' @param res Output of \code{run_deseq}.
#' @param out_dir String depicting the output directory to place
#'     tables and plots.
#' @return None. Plots and tables and placed into output directory.
pca_and_matrices <- function(res, out_dir = "deseq2") {
    ## Initial input matrices
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    res$tab %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        write_tsv(file.path(out_dir, "input_matrix.tsv"))

    res$phenotype %>% write_tsv(file.path(out_dir, "phenotype_data.tsv"))

    saveRDS(res, file.path(out_dir, "deseq2obj.rds"))

    res$sub %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        write_tsv(file.path(out_dir, "input_matrix.filt.tsv"))

    counts(res$ddsObj, normalized = TRUE) %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        write_tsv(file.path(out_dir, "input_matrix.filt.normalized.tsv"))

    assay(res$vFalse) %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        write_tsv(
            file.path(
                out_dir,
                "input_matrix.filt.normalized.vst_corrected.tsv"
            )
        )

    pdf(
        file.path(
            out_dir,
            "input_matrix.filt.normalized.vst_corrected_PCA.pdf"
        ),
        width = 7, height = 7
    )
    plot(res$bFalse)
    plot(res$bTrue)
    dev.off()
}

#' @title Top N Genes
#' @description From a DESeq2 Result extract the top N genes as
#'     determined by descending order of the minus Log10 Adjusted
#'     P-value. Write the result to file as well as return it.
#' @param dsqres the output of the DESeq2 function \code{results},
#'     usually called after performing a contrast.
#' @param top_ng A positive integer for the N value of genes to
#'     extract.
#' @param out_dir String depicting the output directory to place
#'     tables and plots.
#' @param prefix String to act as the prefix filename.
#' @return A vector of N gene names
#' @examples
#' dsqres <- data.frame(mLog10Padj = 1:10 / 10, gene = paste0("G", 1:10))
#' RNASeqHelper:::top_n_genes(dsqres, 3, out_dir = "/tmp", prefix = "test")
top_n_genes <- function(dsqres, top_ng, out_dir, prefix) {
    tgenes <- (dsqres %>% arrange(desc(.data[["mLog10Padj"]])) %>%
                head(top_ng))$gene
    write_tsv(data.frame(tgenes = tgenes),
                file.path(out_dir, paste0(prefix, top_ng, ".tsv")))
    return(tgenes)
}

#' @title Pairwise heatmaps and volcano plots
#' @description For a given contrast and DESeq object produce volcano
#'     plots of the DE genes, and generate a global heatmap with those
#'     DE genes and a pairwise one with just the numerator and
#'     denominator samples.
#' @param ddsObj A DESeq object,
#' @param transformed_counts REALLY UNSURE WHAT'S GOING ON HERE
#' @param numer A string value to select sample columns from the
#'     \code{condition} column that will be part of the numerator part
#'     of the contrast.
#' @param denom A string value to select sample columns from the
#'     \code{condition} column, that will be the denominator part of
#'     the contrast.
#' @param top_ngenes_tocluster A positive integer. How many of the top
#'     log10 padj values to take. Default is 2000.
#' @param top_ngenes_tohighlight A positive integer. How many of the
#'     top log10 padj values to label in the plots. Default is 50.
#' @param score_thresh Vector of numerics depicting cluster score
#'     thresholds to filter for high quality genes in each cluster
#'     before plotting. Default is \code{c(0, 0.5, 0.9, 0.99)}
#' @param genes_of_interest A named list of gene groups to plot in
#'     ribbot plots. If NULL
#' @param volcano_params A list of two components: "global" and
#'     "zoomed" each with "curve" and "curve_show" components
#'     describing the curve to fit to the volcano plots by \code{sd}
#'     (standard deviation from center mean), \code{sc} vertical scale
#'     factor, \code{offset} vertical offset. Extra params for the
#'     \code{zoomed} component are \code{ylim}, \code{padj} an upper
#'     limit on the adjusted p-value, \code{lfc} and an upper-limit on
#'     the log2 fold change. If \code{curve_show} is set to FALSE in
#'     either \code{global} or \code{zoomed} then no curve is shown.
#' @param kmeans_list A list of positive integers to do k-means
#' @param out_dirprefix A character prefix outlining the directory and
#'     basename in which plots and tables will be deployed.
#' @return Void. Plots are deposited to the output directory.
pairwise_hmap_volcano <- function(ddsObj, transformed_counts = NULL,
                                    numer, denom,
                                    top_ngenes_tocluster = 2000,
                                    top_ngenes_tohighlight = 50,
                                    score_thresh = c(0.5, 0.9),
                                    genes_of_interest = NULL,
                                    volcano_params = NULL,
                                    kmeans_list = c(2,5,8,16),
                                    out_dirprefix = ".") {
    message("Started Analysis: ", date())
    plot_title <- paste0(numer, " vs ", denom)
    outdir <- file.path(out_dirprefix, tolower(gsub("[^A-Za-z0-9]", "_",
                                                    plot_title)))
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    ## 1. Perform RNA-seq contrast between numerator and denominator
    dsqres <- results(ddsObj, contrast = c("condition", numer, denom),
                        cooksCutoff = Inf, independentFiltering = FALSE) %>%
        as.data.frame %>% rownames_to_column("gene") %>%
        mutate(mLog10Padj = -log10(.data[["padj"]]))

    norm_counts <- counts(ddsObj, normalized = TRUE)
    top_genes_tocluster <- top_n_genes(dsqres, top_ngenes_tocluster,
                                    outdir, "clustered_genes.top")
    top_genes_tohighlight <- top_n_genes(dsqres, top_ngenes_tohighlight,
                                    outdir, "volcano_genes.tophighlight")
    do_volcanos(dsqres, top_genes_tohighlight, plot_title, outdir,
                volcano_params)
    ## Sample columns are selected from numer and denom in the condition
    pheno <- colData(ddsObj)
    sample_columns <- c(
        pheno[pheno$condition == numer, ]$sample,
        pheno[pheno$condition == denom, ]$sample
    )
    genes <- list(
        cluster = top_genes_tocluster, highlight = top_genes_tohighlight,
        interest = genes_of_interest,  score_thresh = score_thresh
    )
    zzz <- lapply(kmeans_list, function(k) {
        message("[Running k=", k, "]")
        do_kmeans(norm_counts, transformed_counts, colData(ddsObj),
                genes, sample_columns, dsqres, k, outdir)
    })
    NULL
    message("Finished Analysis: ", date())
}

#' @title Perform K-means Pairwise and Global
#' @description Produce two K-means plots where one contains pairwise
#'     comparisons between samples of interest, and one contains the
#'     same pairwise comparisons but extended over all samples.
#' @param norm_counts a matrix of normalized counts
#' @param transformed_counts a matrix of transformed counts, VST or
#'     other produced by DESeq2
#' @param phenotype_data a data.frame of phenotype data which has the
#'     same rows as the number of samples.
#' @param genes A list with at least 4 components: \code{cluster}, a vector
#'     of genes to use for clustering, \code{highlight}, a vector of genes
#'     to highlight in the heatmap; \code{interest}, a list of genes to
#'     plot, where if NULL, use \code{highlight}, and if FALSE do not use;
#'     \code{score_thresh} a vector of numerics depicting cluster score
#'     thresholds to filter for high quality genes in each cluster
#'     before plotting. Default is \code{c(0, 0.5, 0.9, 0.99)}
#' @param sample_columns a vector of sample names. If NULL (the
#'     default), then use all samples
#' @param dsqres the output of the DESeq2 function \code{results},
#'     usually called after performing a contrast.
#' @param k A positive integer to do k-means.
#' @param out_dir String depicting the output directory to place
#'     tables and plots.
#' @return None. Produces only tables and plots in the output
#'     directory.
do_kmeans <- function(norm_counts, transformed_counts, phenotype_data,
                    genes, sample_columns, dsqres, k, out_dir) {
    ## Pairwise Heatmap of Normalised and Transformed Counts
    p1<- kmeans_heatmaps(
        norm_counts, transformed_counts, phenotype_data,
        genes = genes,
        sample_columns = sample_columns,
        dsqres, k,
        out_dir = file.path(out_dir, paste0("kmeans", k)),
        heatprefix = "heatmap_pairwise",
        plot_title = "Heatmap Pairwise")

    ## This below will produce different scaled matrices than the
    ## one above. We can't really pass these back. All we can do
    ## is incorporate the gene plotting stuff into the tail end
    ## of this pipeline, meaning there will be gene plots for
    ## each kmeans heatmap.

    ## Global Heatmap of Normalised and Transformed Counts
    p2 <- kmeans_heatmaps(
        norm_counts, transformed_counts, phenotype_data,
        genes = genes,
        sample_columns = NULL,
        dsqres, k,
        out_dir = file.path(out_dir, paste0("kmeans", k)),
        heatprefix = "heatmap_all",
        plot_title = "Heatmap All")
}

#' @title Nice K-means Heatmap
#' @description Produce a nice K-means clustered heatmap for both
#'     normalized and transformed (DESeq2 VST corrected) matrices
#' @param norms a matrix of normalized counts
#' @param trans a matrix of transformed counts, VST or other produced
#'     by DESeq2.
#' @param pheno a data.frame of phenotype data
#' @param sample_columns a vector of sample names. If NULL (the
#'     default), then use all samples
#' @param genes A list with at least 4 components: \code{cluster}, a vector
#'     of genes to use for clustering, \code{highlight}, a vector of genes
#'     to highlight in the heatmap; \code{interest}, a list of genes to
#'     plot, where if NULL, use \code{highlight}, and if FALSE do not use;
#'     \code{score_thresh} a vector of numerics depicting cluster score
#'     thresholds to filter for high quality genes in each cluster
#'     before plotting. Default is \code{c(0, 0.5, 0.9, 0.99)}
#' @param dsqres the output of the DESeq2 function \code{results},
#'     usually called after performing a contrast.
#' @param k A positive integer to do k-means.
#' @param out_dir String depicting the output directory to place
#'     tables and plots.
#' @param heatprefix String to prefix heatmap plot filenames
#' @param plot_title String depicting title to embed into plot
kmeans_heatmaps <- function(norms, trans, pheno,
                            genes, sample_columns = NULL,
                            dsqres, k, out_dir, heatprefix, plot_title) {
    if (is.null(genes$cluster)) {
        message(" - Using all genes in normalised matrix for clustering")
        genes$cluster <- rownames(norms)
    }
    if (is.null(sample_columns)) {
        message(" - Heatmap all samples in normalised matrix for clustering")
        sample_columns <- colnames(norms)
    } else {
        if (length(sample_columns) <= 10) {
            message(" - Heatmap on samples: ",
                    paste0(sample_columns, collapse = ","))
        } else {
            message(" - Heatmap on ", length(sample_columns), " samples")
        }
    }
    message("   - Using normalized counts")
    res_dsqres <- heatmap_with_geneplots(
        norms[genes$cluster, sample_columns], pheno, k = k, genes = genes,
        out_dir = out_dir,
        heatprefix = heatprefix, prefix_title = paste0(plot_title, " :"),
    )
    dsq_dsq <- left_join(dsqres, res_dsqres$clusters,
                        by = c("gene" = "gene")) %>%
        mutate(norm_cluster = case_when(
                is.na(.data[["cluster"]]) ~ "not in norm heatmap",
                TRUE ~ .data[["cluster"]])) %>% select(-.data[["cluster"]])
    if (!is.null(trans)) { ## Heatmaps using Corrected Normalized Counts
        message("   - Using transformed counts too")
        res_dsqres_trans <- heatmap_with_geneplots(
            trans[genes$cluster, sample_columns], pheno, k = k,
            genes = genes, out_dir = out_dir,
            heatprefix = paste0(heatprefix, ".vst_corrected"),
            prefix_title = paste0(plot_title, " (vst corrected) :")
        )
        dsq_dsq <- left_join(dsq_dsq, res_dsqres_trans$clusters,
                                by = c("gene" = "gene")) %>%
            mutate(trans_cluster = case_when(
                is.na(.data[["cluster"]]) ~ "not in trans heatmap",
                TRUE ~ .data[["cluster"]])) %>% select(-.data[["cluster"]])
    }
    save_cluster <- file.path(out_dir,
                            paste0("deseq2.results.cluster.k", k, ".tsv"))
    write_tsv(dsq_dsq, save_cluster)
    message("   - Storing Results k", k, ":", save_cluster)
}

#' @title Heatmap with Gene plots
#' @description Generate a clustered heatmaps for a specific k means
#'     value.
#' @param norm_counts A matrix containing normalised values.
#' @param pheno_data A data.frame of phenotype data.
#' @param k An integer for the k-value for kmeans.
#' @param genes A list with at least 3 components: \code{highlight}, a
#'     vector of genes to highlight in the heatmap; \code{interest}, a
#'     list of genes to plot, where if NULL, use \code{highlight}, and
#'     if FALSE do not use; \code{score_thresh} a vector of numerics
#'     depicting cluster score thresholds to filter for high quality
#'     genes in each cluster before plotting. Default is \code{c(0,
#'     0.5, 0.9, 0.99)}
#' @param out_dir A character sequence depicting the directory
#' @param heatprefix A character sequence depicting the prefix for
#'     heatmaps
#' @param prefix_title A string to prefix the title of heatmap.
#' @param width_in A positive integer for the number of inches in the
#'     plot width.
#' @param height_in A positive integer for the number of inches in the
#'     plot height.
#' @return A list of two components; clustered tables, and scaled
#'     matrix.
#' @examples
#' n <- 100
#' norm_counts <- matrix(rnorm(n**2, mean = 5), nrow = n)
#' pheno <- data.frame(sample=paste0("S", 1:n),
#'                     condition = c(rep("red", n / 2), rep("green", n / 2)),
#'                     time = as.integer(rnorm(n, 2, 0.5) + 1) * 5)
#' rownames(norm_counts) <- paste0("G", 1:n)
#' colnames(norm_counts) <- paste0("S", 1:n)
#' genes <- list(highlight = paste0("G", 1:10),
#'              interest = paste0("G", c(2,3)),
#'              score_thresh = c(0.2, 0.3))
#' res <- heatmap_with_geneplots(norm_counts, pheno, 2,
#'                               genes, out_dir = "/tmp")
#' @export
heatmap_with_geneplots <- function(norm_counts, pheno_data, k,
                                    genes, out_dir = "heatmaps_k",
                                    heatprefix = "heatmap",
                                    prefix_title = "",
                                    width_in = 6, height_in = 6) {
    ## normalised counts should be in the correct sample order already
    options(repr.plot.height = height_in, repr.plot.width = width_in)
    if (!dir.exists(out_dir)) { dir.create(out_dir) }
    if (is.null(genes$highlight)) {
        top_genes <- head(names(sort(rowMeans(norm_counts),   ## If no genes
                                    decreasing = TRUE)), 30)   ## given show
        top_title <- paste0(nrow(norm_counts), " DE genes, top ", ## top N
                            length(top_genes), " highlighted")
    } else {
        top_genes <- genes$highlight
        top_title <- paste0(nrow(norm_counts), " DE genes, ",
                            length(top_genes), " highlighted")
    }
    if (is.null(genes$interest)) {
        genes$interest <- list(topgenes = top_genes)
    }
    scale_mat <- t(scale(t(norm_counts)))
    scale_mat <- scale_mat[complete.cases(scale_mat), ] ## Remove 0s

    cluster_assignments <- single_kmeans_heatmap(
        scale_mat, k, top_genes, prefix_title, top_title,
        out_dir, heatprefix, width_in, height_in)
    cluster_assignments_scores <- calculate_cluster_corr(
        cluster_assignments, scale_mat, out_dir, heatprefix)

    save_norm <- file.path(out_dir, paste0("k", k, "_norm.tsv"))
    save_scale <- file.path(out_dir, paste0("k", k, "_scale.tsv"))
    write_tsv(as.data.frame(norm_counts) %>%
                rownames_to_column("gene"), save_norm)
    message("     - Saved Norm: ", save_norm)
    write_tsv(as.data.frame(scale_mat) %>% rownames_to_column("gene"),
                save_scale)
    message("     - Saved Scale: ", save_norm)
    if (!is.null(genes$interest)) {
        do_gene_plots(norm_counts, scale_mat, pheno_data,
                        cluster_assignments_scores,
                        score_thresh = genes$score_thresh,
                        genes_of_interest = genes$interest, out_dir)
    }
    return(list(clusters = cluster_assignments_scores, scaled = scale_mat))
}

#' @title Calculate Cluster Correlation for a single cluster
#' @description Calculate cluster correlation scores for a single
#'     cluster of genes
#' @param clust_assign A two-column table of 'gene' and 'cluster'.
#' @param scale_mat A matrix with genes as rows and samples as
#'     columns.
#' @param i A positive integer representing a cluster number.
#' @return A three-column table of 'gene', 'cluster' and 'score'.
#' @examples
#' n <- 100
#' scale_mat <- matrix(rnorm(n**2, mean = 0), nrow = n)
#' rownames(scale_mat) <- paste0("G", 1:n)
#' colnames(scale_mat) <- paste0("S", 1:n)
#' clust_assign <- data.frame(
#'     gene = rownames(scale_mat),
#'     cluster = c(rep(1, n / 50), rep(2, n / 50)),
#'     value = rnorm(n, 100, 10)
#' )
#' cac <- calculate_cluster_corr_i(clust_assign, scale_mat, 2)
#' @export
calculate_cluster_corr_i <- function(clust_assign, scale_mat, i) {
    clust_sub <- clust_assign[clust_assign$cluster == i, ]
    genes_in_i <- clust_sub$gene
    matrix_in_i <- scale_mat[genes_in_i, ]
    mean_of_samples_in_i <- colMeans(matrix_in_i)

    sample_centroids_in_i <- data.frame(centroid = mean_of_samples_in_i) %>%
        rownames_to_column("sample") %>%
        mutate(cluster = i) %>%
        select(c(.data[["sample"]], .data[["cluster"]], .data[["centroid"]]))

    corr_basis <- left_join(
        x = sample_centroids_in_i,
        y = (as.data.frame(t(matrix_in_i)) %>% rownames_to_column("sample")),
        by = "sample"
    )
    ## Two column table of "gene" and "score"
    corr_scores <- data.frame(score = vapply(genes_in_i, function(gene) {
        cor(corr_basis[["centroid"]], corr_basis[[gene]])
    }, 0)) %>% rownames_to_column("gene")

    return(left_join(clust_sub, corr_scores, by = "gene"))
}

#' @title Calculate Cluster Correlation
#' @description Calculate mean expression of each sample for each
#'     cluster to have gene centroids, and use these to calculate
#'     correlation scores for each gene against the centroid.
#' @param clust_assign A two-column table of 'gene' and 'cluster'.
#' @param scale_mat A matrix with genes as rows and samples as
#'     columns.
#' @param out_dir A string denoting the output directory to store
#'     plots and tables.
#' @param prefix_str A string prefix for the filename
#' @return A three-column table of 'gene', 'cluster' and 'score'.
#' @examples
#' n <- 100
#' scale_mat <- matrix(rnorm(n**2, mean = 0), nrow = n)
#' rownames(scale_mat) <- paste0("G", 1:n)
#' colnames(scale_mat) <- paste0("S", 1:n)
#' clust_assign <- data.frame(
#'     gene = rownames(scale_mat),
#'     cluster = c(rep(1, n / 50), rep(2, n / 50)),
#'     value = rnorm(n, 100, 10)
#' )
#' ca <- calculate_cluster_corr(clust_assign, scale_mat, "/tmp", "red")
#' @export
calculate_cluster_corr <- function(clust_assign, scale_mat,
                                    out_dir, prefix_str) {
    all_clusters <- unique(sort(clust_assign$cluster))
    clust_assignments <- do.call(rbind, lapply(all_clusters, function(cl) {
        calculate_cluster_corr_i(clust_assign, scale_mat, cl)
    }))
    write_tsv(clust_assignments,
        file = file.path(
            out_dir,
            paste0(prefix_str, ".cluster_assignments.tsv")
        )
    )
    return(clust_assignments)
}


#' @title Perform Gene Plots
#' @description Perform cluster gene plots for normalised and scaled
#'     data as well as plotting genes of interest for both.
#' @param norm_counts a matrix of normalized count data, with genes as
#'     rows and samples as columns.
#' @param scale_mat a matrix of scaled count data, with genes as rows
#'     and samples as columns.
#' @param pheno_data a data.frame of phenotype data.
#' @param gene_cluster_scores A data frame with three columns: 'gene',
#'     'cluster', and 'scores'.
#' @param score_thresh Vector of numerics depicting cluster score
#'     thresholds to filter for high quality genes in each cluster
#'     before plotting. Default is \code{c(0, 0.5, 0.9, 0.99)}
#' @param genes_of_interest A named list of gene vectors. A list of
#'     genes which are referenced by specific names. Unique plots will
#'     be generated for list. e.g. list(mygeneset1 = c("Msgn1",
#'     "Osr1", "Rspo3", "Fgf8", "Wnt3a"), mygeneset2 = c("Mesp1",
#'     "Foxa2", "Sox17", "Lhx1", "Cer1"))
#' @param out_dir A string denoting the output directory to store
#'     plots and tables.
#' @return None.
#' @examples
#' n <- 100
#' norm_counts <- matrix(rnorm(n**2, mean = 5), nrow = n)
#' rownames(norm_counts) <- paste0("G", 1:n)
#' colnames(norm_counts) <- paste0("S", 1:n)
#' scale_mat <- matrix(rnorm(n**2, mean = 0), nrow = n)
#' rownames(scale_mat) <- paste0("G", 1:n)
#' colnames(scale_mat) <- paste0("S", 1:n)
#' pheno <- data.frame(
#'     sample = paste0("S", 1:n),
#'     condition = c(rep("red", n / 2), rep("green", n / 2)),
#'     time = as.integer(rnorm(n, 2, 0.5) + 1) * 5
#' )
#' gene_cluster_scores <- data.frame(
#'     gene = rownames(scale_mat),
#'     cluster = c(rep(1, n / 50), rep(2, n / 50)),
#'     score = rnorm(n, 0.5, 0.25)
#' )
#' score_thresh <- c(0.2, 0.3)
#' genes_of_interest <- paste0("G", 5:50)
#' res <- do_gene_plots(norm_counts, scale_mat, pheno,
#'                      gene_cluster_scores, score_thresh,
#'                      genes_of_interest,
#'                      out_dir = "/tmp")
#' @export
do_gene_plots <- function(norm_counts, scale_mat, pheno_data,
                        gene_cluster_scores, score_thresh,
                        genes_of_interest, out_dir) {

    ## bind the phenotype data to the matrix data
    long_norm <- left_join(
        left_join(
            as.data.frame(t(norm_counts)) %>% rownames_to_column("sample"),
            as.data.frame(pheno_data), by = "sample" ) %>%
        pivot_longer(all_of(rownames(norm_counts)), names_to = "gene"),
        gene_cluster_scores, by = "gene" )

    long_scale <- left_join(
        left_join(
            as.data.frame(t(scale_mat)) %>% rownames_to_column("sample"),
            as.data.frame(pheno_data), by = "sample" ) %>%
        pivot_longer(all_of(rownames(scale_mat)), names_to = "gene"),
        gene_cluster_scores, by = "gene" )

    gene_clusters_by_score(long_norm,
                            score_thresh = score_thresh,
                            out_dir = file.path(out_dir, "gene_cluster_norm"))
    message("     - Plotting genes Norm")
    gene_clusters_by_score(
        long_scale,
        score_thresh = score_thresh,
        out_dir = file.path(out_dir, "gene_cluster_scale"))
    message("     - Plotting genes Scaled")

    gene_plots_by_gene(long_norm, long_scale, genes_of_interest,
                    outprefix = "gene.lists",
                    out_dir = file.path(out_dir, "gene_lists"))
    message("     - Plotting genes List")
}

#' @title Cluster Assignments
#' @description Get the cluster assignment for all genes
#' @param hm_now_drawn A ComplexHeatmap object invoked with
#'     \code{draw}.
#' @param matobj A matrix with genes as rownames.
#' @return A two column data frame of Gene and Cluster in ascending
#'     order of cluster number.
cluster_assignments <- function(hm_now_drawn, matobj) {
    cluster_list <- lapply(
        row_order(hm_now_drawn),
        function(x) rownames(matobj)[x]
    )
    cluster_table <- do.call(
        rbind, lapply(
            names(cluster_list),
            function(n) {
                data.frame(
                    gene = cluster_list[[n]],
                    cluster = n
                )
            }
        )
    )
    return(cluster_table)
}

#' @title A single K-means Heatmap
#' @description Plot a clustered heatmap
#' @param scale_mat A matrix of scaled values.
#' @param k A positive integer for the k value of kmeans
#' @param top_genes A list of genes to plot
#' @param prefix_title A string to prefix the title of heatmap.
#' @param top_title A character string for the title of the heatmap.
#' @param out_dir A character sequence depicting the directory
#' @param heatprefix A character sequence depicting the prefix for
#'     heatmaps
#' @param width_in A positive integer for the number of inches in the
#'     plot width.
#' @param height_in A positive integer for the number of inches in the
#'     plot height.
#' @return A list of two components; clustered tables, and scaled
#'     matrix.
#' @examples
#' n <- 100
#' scale_mat <- matrix(rnorm(n**2, mean = 5), nrow = n)
#' rownames(scale_mat) <- paste0("G", 1:n)
#' colnames(scale_mat) <- paste0("S", 1:n)
#' top_genes <- paste0("G", n - 50:n - 20)
#' res <- single_kmeans_heatmap(scale_mat, 2, top_genes, "test", "test",
#'                              out_dir = "/tmp", "test", 7, 7)
#' @export
single_kmeans_heatmap <- function(scale_mat, k, top_genes,
                                    prefix_title, top_title,
                                    out_dir, heatprefix,
                                    width_in, height_in) {
    ha <- rowAnnotation(
        foo = anno_mark(
            at = which(rownames(scale_mat) %in% top_genes),
            labels = top_genes))
    if (k > 5) { rtitle <- "%s" }
    else { rtitle <- "Clust %s" }
    hm_now <- Heatmap(
        scale_mat, column_title = paste0(prefix_title, "  ", top_title),
        row_km = k, cluster_row_slices = FALSE, ## Arrange the K? NO
        cluster_rows = TRUE,  ## Prevents annotation bunching if TRUE
        show_row_dend = FALSE, row_gap = unit(3, "mm"), name = "mat",
        row_title = rtitle, right_annotation = ha,
        col = colorRampPalette(c("black", "#bb0000"))(100),
        cluster_columns = FALSE, show_row_names = FALSE,
        column_names_rot = 45, row_names_gp = gpar(fontsize = 4))

    hm_now_drawn <- draw(hm_now)
    cluster_assignments <- cluster_assignments(hm_now_drawn, scale_mat)
    ## Unite cluster assignments with Norm, Scaled, and Pvalue
    output_prefix <- file.path(out_dir,
                                paste0(heatprefix, ".k", k, "."))
    save_svg <- paste0(output_prefix, "svg")
    ## save_svgold <- paste0(output_prefix, "_old.svg")
    save_pdf <- paste0(output_prefix, "pdf")
    ## save_png <- paste0(output_prefix, "png")
    save_scale <- paste0(output_prefix, "scale.tsv")
    save_norm <- paste0(output_prefix, "norm.tsv")
    ## save_cluster <- paste0(output_prefix, "clusters.tsv")
    ##bhm <- better_pheatmap(hm_now)
    ##svg(save_svgold, width = width_in, height = height_in)
    ##print(hm_now)
    ##dev.off()

    svg(save_svg, width = width_in, height = height_in)
    plot(hm_now)     ##grid.draw(bhm);
    dev.off()

    pdf(save_pdf, width = width_in, height = height_in)
    plot(hm_now)     ##grid.draw(bhm); grid.newpage()
    dev.off()
    return(cluster_assignments)
}

#' @title Better Pheatmaps
#' @description Remove the greylines that plagues these plots by
#'     modifying the raw ggplot grob.
#' @param ph a pheatmap object
#' @return a modified pheatmap object with lines filled with fill
#'     colours.
better_pheatmap <- function(ph) {
    lwdd <- 2
    message("-[better_pheatmap]-")
    if (inherits(ph) == "Heatmap") { ## ComplexHeatmap
        ph <- grid.grabExpr(draw(ph))
        hmaps <- ph$children
        for (hm in names(hmaps)) { ## Get all Rects
            rects <- names(ph$children[[hm]]$children)
            if (length(rects) > 0) {
                ## rects = names(ph$children[[hm]]$children)
                for (rr in rects) { ## Check for rects with tabular data
                    hasdim <- dim(ph$children[[hm]]$
                        children[[rr]]$gp$fill)
                    if (!is.null(hasdim)) {
                        ph$children[[hm]]$children[[rr]]$gp$col <-
                            ph$children[[hm]]$children[[rr]]$gp$fill
                        ph$children[[hm]]$children[[rr]]$gp$lwd <- lwdd
                    }
                }
            } else {
                if ("gp" %in% names(ph$children[[hm]])) {
                    ## Check for rects with tabular data
                    hasdim <- dim(ph$children[[hm]]$gp$fill)
                    if (!is.null(hasdim)) {
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
        grob_classes <- map(ph$gtable$grobs, class) ## Extract the right grob
        idx_grob <- which(
            map_lgl(grob_classes, function(cl) "gTree" %in% cl)
        )[1]
        grob_names <- names(ph$gtable$grobs[[idx_grob]]$children)
        idx_rect <- grob_names[grep("rect", grob_names)][1]
        ## Remove borders around cells
        ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <-
            ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$fill
        ph$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 3
        return(ph)
    }
}

#' @title Get a high quality genes list
#' @description Filter a wide-format matrix for high quality genes
#'     (rows) based on detectability thresholds and smallest
#'     occurrence.
#' @param sam_mat A numeric matrix. Genes are rows and Samples as
#'     columns. Typically this is the average sample matrix, with
#'     replicate values summarized into mean expression.
#' @param min_occur A positive integer. A threshold (greater than or
#'     equal to) the number of times a gene is detected in a sample,
#'     across all samples (e.g. Gene X is only kept if it appears in
#'     10 samples with a value greater than the \code{min_detect}
#'     threshold. Default value is 3.
#' @param min_detect A positive integer. A threshold (greater than or
#'     equal to) for the minimum expression a gene can have for it to
#'     be a valid occurrence in a sample. Default values is 10.
#' @param out_dir A string depicting the output directory to place
#'     tables.
#' @return A vector of strings depicting a list of high quality genes.
#'     A table of rowSums is written to a file matching
#'     "smallestGroup*-detected*-keep_genes".
#' @examples
#' sam_mat <- matrix(1:100, nrow = 10, dimnames = list(
#'     paste0("G", 1:10),
#'     paste0("C", 1:10)
#' ))
#' keep <- high_quality_genes(sam_mat, 10, 7, "/tmp")
#' names(keep) == paste0("G", 7:10)
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
    message(
        "Dropped: ", length(drop_genes),
        " genes (",
        as.integer(100 * length(drop_genes) / length(res)), "%)"
    )
    return(keep_genes)
}


#' @title Do volcano plots
#' @description Generate General and Zoomed volcano plots
#' @param dsqres A DESeq2 Result object
#' @param top_genes_tohighlight A vector of genes to highlight in the
#'     volcano plots
#' @param plot_title A string depicting the title of the gplots
#' @param outdir An output directory path fro the plots
#' @param volcano_params A list of two components: "global" and
#'     "zoomed" each with "curve" and "curve_show" components
#'     describing the curve to fit to the volcano plots by \code{sd}
#'     (standard deviation from center mean), \code{sc} vertical scale
#'     factor, \code{offset} vertical offset. Extra params for the
#'     \code{zoomed} component are \code{ylim}, \code{padj} an upper
#'     limit on the adjusted p-value, \code{lfc} and an upper-limit on
#'     the log2 fold change. If \code{curve_show} is set to FALSE in
#'     either \code{global} or \code{zoomed} then no curve is shown.
#' @return None.
#' @examples
#' n <- 100
#' dsqres <- data.frame(
#'     gene = paste0("G", 1:n),
#'     data = rnorm(n, 5, 2),
#'     log2FoldChange = 1e-5,
#'     mLog10Padj = rnorm(n, 1000, 2),
#'     padj = 1/rnorm(n, 1000, 2)
#' )
#' top_genes_tohighlight <- paste0("G", (n-50):(n-20))
#' res <- RNASeqHelper:::do_volcanos(dsqres, top_genes_tohighlight,
#'                                   "my plot", "/tmp")
do_volcanos <- function(dsqres, top_genes_tohighlight, plot_title, outdir,
                        volcano_params = NULL) {
    if (is.null(volcano_params)) {
        volcano_params <- list(
            global = list(
                curve = list(sd = 0.3, sc = 60, offset = 10),
                curve_show = FALSE
            ),
            zoomed = list(
                curve = list(sd = 0.25, sc = 5, offset = 8), ylim = c(0, 22),
                curve_show = FALSE, padj = 20, lfc = 1.5
            )
        )
    }
    ## Volcano Plots
    p1 <- volcano_plot(dsqres, top_genes_tohighlight, plot_title,
        curve = volcano_params$global$curve,
        curve_show = volcano_params$global$curve_show
    )
    ## Volcano Plots zoomed in
    p2 <- volcano_plot(
        dsqres %>%
        filter(abs(.data[["log2FoldChange"]]) <
                volcano_params$zoomed$lfc &
                .data[["mLog10Padj"]] < volcano_params$zoomed$padj),
        top_genes_tohighlight, paste0(plot_title, " (zoomed)"),
        curve = volcano_params$zoomed$curve,
        ylim = volcano_params$zoomed$ylim,
        curve_show = volcano_params$zoomed$curve_show
    )

    hmaps <- wrap_plots(list(p1, p2), ncol = 1, guides = "collect") &
        theme_bw() + theme(legend.position = "none")

    volcano_svg <- file.path(outdir, "volcano_pairwise.svg")
    deseq2_out <- file.path(outdir, "deseq2_results.tsv")

    svg(volcano_svg, width = 8, height = 9)
    plot(hmaps)
    dev.off()
    message("Saved Volcano: ", volcano_svg)

    dsqres %>%
        mutate(isTopN.gene = .data[["gene"]] %in%
            top_genes_tohighlight) %>%
        write_tsv(deseq2_out)
    message("Saved DESeq2 Results: ", deseq2_out)
}

#' @title Produce a Volcano Plot
#' @description Generates a volcano plot from the contrast of a DESeq2
#'     analysis with added curves for graphic effect.
#' @param dsqres A table of DESeq2 results
#' @param degenes A vector of genes.
#' @param title A string depicting the title of the plot
#' @param curve A list of three components: sd, sc, offset.
#' @param curve_scale A numeric curve scaling factor.
#' @param curve_show A boolean on whether or not to show the curve.
#' @param ylim A vector of two components to limit the Y-axis.
#' @return A ggplot2 object of the volcano plot.
#' @examples
#' n <- 100
#' dsqres <- data.frame(gene = paste0("G", 1:n),
#'                      data = rnorm(n, 5, 2),
#'                      log2FoldChange = 1e-5,
#'                      padj = 1 / rnorm(n, 5, 2))
#' degenes <- paste0("G", (n-50):(n-20))
#' res <- volcano_plot(dsqres, degenes, "my plot")
#' @export
volcano_plot <- function(dsqres, degenes, title,
                        curve = list(sd = 0.15, sc = 10, offset = 1),
                        curve_scale = 1, curve_show = FALSE, ylim = NULL) {
    options(repr.plot.height = 12, repr.plot.width = 12)
    ## Extract relevant info from DESeq results
    ana <- dsqres %>% select(c(.data[["gene"]], .data[["log2FoldChange"]],
                                .data[["padj"]])) %>%
        mutate(mLog10Padj = -log10(.data[["padj"]])) %>%
        arrange(desc(.data[["mLog10Padj"]]))

    cust_fun <- function(x, sd = 0.15, sc = 10, offset = 1) {
        offset + dnorm(x, mean = 0, sd = sd) * sc
    }
    cfun <- function(x) { ## The slope function -- highlight these genes
        return(cust_fun(x,
            sd = curve$sd * curve_scale,
            sc = curve$sc * curve_scale, offset = curve$offset
        ))
    }
    ## We highlight genes in the zoomed zone fitting the curve, but
    ## the main DE genes are shown as shapes and highlighted
    red <- ana %>%
        mutate(isTopN = .data[["gene"]] %in% degenes) %>%
        mutate(highlight = .data[["isTopN"]] |
            (.data[["mLog10Padj"]] >
                cust_fun(.data[["log2FoldChange"]])))

    max_x <- max(abs(ana$log2FoldChange)) + 0.05 ## symmetry
    plot1 <- red %>% ggplot(aes_string(
                        x = "log2FoldChange", y = "mLog10Padj",
                        colour = "highlight", shape = "isTopN",
                        label = "gene")) + geom_point() +
        scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
        scale_shape_manual(values = c("TRUE" = 5, "FALSE" = 19)) +
        scale_x_continuous(limits = c(-max_x, max_x), breaks = waiver(),
                        n.breaks = 10) +
        geom_label_repel(data = red %>%
                            filter(.data[["highlight"]] == TRUE) %>%
                            head(15), box.padding = 0.5, max.overlaps = 30,
                        colour = "black") + ggtitle(title)

    if (curve_show) {
        plot1 <- plot1 + geom_function(fun = cfun, n = 100, colour = "blue")
    }
    if (!is.null(ylim)) {
        plot1 <- plot1 + scale_y_continuous(limits = ylim)
    }
    return(plot1)
}

#' @title Gene clusters by score
#' @description Generate cluster gene plots for various score
#'     thresholds
#' @param tab A dataframe with long plotting data containing at least
#'     the columns \code{gene}, \code{cluster}, \code{condition},
#'     \code{value}, \code{time}, and \code{score}. It can take
#'     normalised or scaled values as input.
#' @param score_thresh Vector of numerics depicting cluster score
#'     thresholds to filter for high quality genes in each cluster
#'     before plotting. Default is \code{c(0, 0.5, 0.9, 0.99)}
#' @param out_dir String representing the directory to store gene
#'     plots. Default value is "gene_cluster"
#' @return A single patchwork plot containing cluster gene plots for
#'     different score thresholds. This plot is also saved to the
#'     \code{out_dir} directory under the pattern
#'     "gene_plots-k*_montage.svg"
#' @examples
#' n <- 100
#' tab <- data.frame(
#'     gene = paste0("G", 1:n),
#'     cluster = as.integer(rnorm(n, 5, 2) + 1),
#'     condition = c(rep("red", n/2), rep("green", n/2)),
#'     value = rnorm(n, 10, 2) + 1,
#'     time = as.integer(rnorm(n, 2, 0.5) + 1) * 5,
#'     score = as.integer(rnorm(n, 2, 0.5) + 1) / 4
#' )
#' plot <- gene_clusters_by_score(tab, c(0.2, 0.7), out_dir = "/tmp")
#' @export
gene_clusters_by_score <- function(tab,
                                score_thresh = c(0, 0.5, 0.9, 0.99),
                                out_dir = "gene_cluster") {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
    }
    res <- lapply(score_thresh, function(x) {
        return(cluster_gene_plots(tab, score_thresh = x,
                                out_dir = out_dir))
    })
    newoutprefix <- file.path(out_dir,
                            paste0("gene_plots-k", num_clusters(tab)))

    rplot <- wrap_plots(res, ncol = 2, guides = "collect") &
        theme(plot.title = element_blank(),
            axis.title = element_blank())

    ggsave(plot = rplot, filename = paste0(newoutprefix,
                                            "_montage.svg"),
                dpi = 800, width = 15, height = 15, units = "in")

    return(rplot)
}

#' @title Cluster Gene Plots
#' @description Generate ribbon and line plots of gene expression data
#'     split by cluster.
#' @param tab A dataframe with long plotting data containing at least
#'     the columns \code{cluster}, \code{condition}, \code{value},
#'     \code{time}, and \code{score}. It can take normalised or scaled
#'     values as input.
#' @param score_thresh numeric threshold to filter out genes with
#'     correlation scores not matching these values
#' @param out_dir a string depicting the directory of where the plots
#'     are deposited, with the prefix "gene_plots-k".
#' @return ggplot object ribbon and line plots facetted by cluster.
cluster_gene_plots <- function(tab, score_thresh = 0,
                                out_dir = "gene_cluster") {
    tabn <- tab %>% filter(.data[["score"]] >= score_thresh)
    ## Keep all clusters even if they're empty
    tabn$cluster <- factor(tabn$cluster, levels = unique(sort(tab$cluster)))
    dat_labs <- generate_labelling_table(tabn)
    time_breaks <- sort(unique(sort(tab$time)))
    ## labelling function for facet
    lab_fun <- function(s) {
        (dat_labs %>% filter(.data[["cluster"]] == s))$ftext
    }

    p1 <- tabn %>% ggplot(aes_string(
                        x = "time", y = "value", fill = "condition",
                        colour = "condition", group = "condition")) +
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
                                paste0("gene_plots-k",
                                    num_clusters(tab),
                                    "-score", score_thresh))
    ggsave(plot = p1, filename = paste0(newoutprefix, ".svg"),
            dpi = 800, width = 10, height = 10, units = "in")
    tabn %>% write_tsv(paste0(newoutprefix, ".tsv"))
    return(p1)
}

#' @title Gene plots by gene
#' @description Generate plots for genes of interest showing their
#'     normalized and z-score scaled expression.
#' @param norm_long A long dataframe of normalized values. Must
#'     contain columns of: gene, condition, time, replicate, value.
#' @param scale_long A long dataframe of scaled values. Same format as
#'     \code{norm_long}.
#' @param genes_of_interest A named list of gene vectors. A list of
#'     genes which are referenced by specific names. Unique plots will
#'     be generated for list. e.g. list(mygeneset1 = c("Msgn1",
#'     "Osr1", "Rspo3", "Fgf8", "Wnt3a"), mygeneset2 = c("Mesp1",
#'     "Foxa2", "Sox17", "Lhx1", "Cer1"))
#' @param outprefix A string depicting the output prefix for the
#'     plots. These will be appended with the title of each gene list
#'     and for the normalized and scaled.
#' @param out_dir A string denoting the output directory to store
#'     plots. Default is is "gene_lists".
gene_plots_by_gene <- function(norm_long, scale_long, genes_of_interest,
                                outprefix = "gene.lists",
                                out_dir = "gene_lists") {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
    }
    pgene_list <- lapply(
        names(genes_of_interest), function(glist_name) {
            glist <- genes_of_interest[[glist_name]]
            genes_found <- unique((norm_long %>%
                                    filter(.data[["gene"]] %in% glist))$gene)
            if (length(genes_found) < 1) {
                message("no genes found for: ", glist_name)
                return(NULL)
            }
            pgene_norm <- single_gene_plot(
                norm_long, genes_found, glist_name,
                "Log10 Normalised Expression",
                "Normalised Expression Plots grouped by Genes of Interest",
                out_dir, outprefix, ".normalised.svg", TRUE
            )
            pgene_scale <- single_gene_plot(
                scale_long, genes_found, glist_name,
                "Scaled Expression",
                "Scaled Expression Plots grouped by Genes of Interest",
                out_dir, outprefix, ".scaled.svg", FALSE
            )
            return(list(norm = pgene_norm, scale = pgene_scale))
        }
    )
}


#' @title Single Gene Plot
#' @description Perform a plot of all genes for a single instance of
#'     long table data and save to file, as well as return the plot.
#' @param long_data A tibble containing columns: gene, time, condition, value
#' @param genes_found A vector of gene names
#' @param glist_name A string for the group name for said vector of genes
#' @param ylab_text A string for Y-axis label text
#' @param title_text A string for the title text
#' @param out_dir A string denoting the output directory to store
#'     plots. Default is is "gene_lists".
#' @param outprefix A string for the filename prefix
#' @param filesuffix A string depicting the filename suffix and extension.
#' @param scaley10 A boolean on whether to scale the Y-axis by log10
#' @return A ggplot2 object
single_gene_plot <- function(long_data, genes_found, glist_name,
                            ylab_text, title_text,
                            out_dir, outprefix, filesuffix, scaley10) {
    time_breaks <- sort(unique(sort(long_data$time)))

    pgene <- long_data %>%
        filter(.data[["gene"]] %in% genes_found) %>%
        ggplot(aes_string(
            x = "time", y = "value", fill = "condition",
            colour = "condition", group = "condition"
        )) +
        stat_summary(
            fun.data = "mean_sdl", geom = "ribbon",
            alpha = 0.1, colour = NA
        ) +
        stat_summary(fun = mean, geom = "line") +
        geom_jitter(width = 0.1) +
        facet_wrap("gene", scales = "free_y") +
        ylab(ylab_text) +
        ggtitle(title_text) +
        scale_x_continuous(breaks = time_breaks) +
        theme_bw()

    if (scaley10) {
        pgene <- pgene + scale_y_log10()
    }
    plot_dims <- function(n) { ## Calculates plot dims for a given N plots
        w <- ceiling(sqrt(n))
        return(list(w = w, h = ceiling(n / w)))
    }
    pdims <- plot_dims(length(genes_found))
    ggsave(
        plot = pgene,
        filename = file.path(
            out_dir, paste0(outprefix, ".", glist_name, filesuffix)
        ),
        dpi = 800, width = pdims$w * 2, height = pdims$h * 1.5,
        units = "in"
    )
    return(pgene)
}


#' @title Generate a label table
#' @description Counts how many genes in each cluster for use as
#'     plotting labels in \code{cluster_gene_plots}
#' @param tabn a gene table of long data with a \code{cluster},
#'     \code{gene} and \code{value} field.
#' @return a table with columns of \code{cluster}, \code{n} for how
#'     many genes in that cluster, and \code{ftext} for labelling.
#' @examples
#' tabn <- data.frame(
#'     cluster = c(1, 1, 2, 2),
#'     gene = paste0("G", 1:4),
#'     value = 1:4
#' )
#' res <- RNASeqHelper:::generate_labelling_table(tabn)
generate_labelling_table <- function(tabn) {
    dat_labs <- tabn %>%
        group_by(.data[["cluster"]]) %>%
        summarise(n = n_distinct(.data[["gene"]])) %>%
        mutate(ftext = paste0("cluster ", .data[["cluster"]], ", ",
                            .data[["n"]], " genes"))

    ## Sanity check that all factors are in summary
    miss <- which(!(levels(tabn$cluster) %in% dat_labs$cluster))
    if (length(miss) > 0) {
        temp <- cbind(cluster = miss, n = 0,
                        ftext = paste0("cluster ", miss, ", 0 genes"))
        dat_labs <- rbind(dat_labs, temp)
    }
    dat_labs <- dat_labs %>% arrange(.data[["cluster"]])
    return(dat_labs)
}

#' @title Number of clusters in tabs
#' @description Count the number of unique values in the
#'     \code{cluster} column of a table
#' @param tab Dataframe. Must contain a \code{cluster} vector
#' @return Positive integer.
#' @examples
#' RNASeqHelper:::num_clusters(
#'     data.frame(cluster = c(1, 1, 2, 3, 2))
#' )
num_clusters <- function(tab) {
    return(length(unique(tab$cluster)))
}

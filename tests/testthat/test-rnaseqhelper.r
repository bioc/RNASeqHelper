
make_mat_and_phenotype <- function(ngenes=1000, nsamples=10,
                                   mean=1000, var=500) {
    tab <- matrix(as.integer(rnorm(ngenes*nsamples, mean, var)),
                  ncol = nsamples)
    tab <- tab - min(tab)
    rownames(tab) <- paste0("G", 1:ngenes)
    colnames(tab) <- paste0("S", 1:nsamples)
    phenotype_data <- data.frame(
        sample = paste0("S", 1:nsamples),
        condition = c(rep("red", nsamples / 2),
                      rep("green", nsamples / 2)),
        time = as.integer(rnorm(nsamples, 2, 0.5) + 1) * 5
    )
    return(list(mat = tab, pheno = phenotype_data))
}

cleanup_outdir <- function(out_dir, numer="red", denom="green") {
    ## Unused. Output directory artefacts come from examples,
    ## not tests.
    out_dir2 <- file.path(out_dir, paste0(numer, "_vs_", denom))
    if (dir.exists(out_dir2)) {
        unlink(out_dir2, recursive = TRUE) ## rm -rf equivalent
        message("Removed: ", out_dir2)
    }
}

test_that("rnaseqhelper", {
    tmp <- make_mat_and_phenotype(ngenes = 1000, nsamples = 10)
    tab <- tmp[["mat"]]
    phenotype_data <- tmp[["pheno"]]

    rnaseqhelper(tab, phenotype_data,
                out_dir = tempdir(), "green", "red",
                heat_params = list(
                    score_thresh = c(0.2, 0.5),
                    kmeans_list = 2))
})

test_that("run_deseq", {
    tmp <- make_mat_and_phenotype(ngenes = 2000, nsamples = 10)
    tab <- tmp[["mat"]]
    phenotype_data <- tmp[["pheno"]]
    keep_genes <- paste0("G", 1:1000)
    res <- run_deseq(tab, keep_genes, phenotype_data, tempdir())
})

test_that("run_deseq_notime", {
    tmp <- make_mat_and_phenotype(ngenes = 2000, nsamples = 10)
    tab <- tmp[["mat"]]
    phenotype_data <- tmp[["pheno"]]
    phenotype_data[["time"]] <- NULL ## Blank out time component
    keep_genes <- paste0("G", 1:1000)
    res <- run_deseq(tab, keep_genes, phenotype_data, tempdir())
})

test_that("top_n_genes", {
    dsqres <- data.frame(mLog10Padj = 1:10 / 10, gene = paste0("G", 1:10))
    RNASeqHelper:::top_n_genes(dsqres, 3, out_dir = tempdir(), prefix = "test")
})

test_that("heatmap_with_geneplots", {
    tmp <- make_mat_and_phenotype(ngenes = 100, nsamples = 10, mean=5, var=1)
    ## convert to decimal
    norm_counts <- apply(tmp[["mat"]], 1:2, as.numeric) * 1.5
    phenotype_data <- tmp[["pheno"]]
    genes <- list(
        highlight = paste0("G", 1:10),
        interest = paste0("G", c(2, 3)),
        score_thresh = c(0.2, 0.3)
    )
    res <- heatmap_with_geneplots(norm_counts, phenotype_data, 2,
                                genes, out_dir = tempdir())
})

test_that("calculate_cluster_corr", {
    ngenes <- 100
    nsamples <- 10
    scale_mat <- matrix(rnorm(ngenes*nsamples, mean = 0), nrow = ngenes)
    rownames(scale_mat) <- paste0("G", seq(ngenes))
    colnames(scale_mat) <- paste0("S", seq(nsamples))
    clust_assign <- data.frame(
        gene = rownames(scale_mat),
        cluster = c(rep(1, ngenes / 50), rep(2, ngenes / 50)),
        value = rnorm(ngenes, 100, 10)
    )
    ca <- calculate_cluster_corr(clust_assign, scale_mat, tempdir(), "red")
})

test_that("do_gene_plots", {
    ngenes <- 100
    tmp <- make_mat_and_phenotype(ngenes=ngenes, nsamples=10, mean = 5, var=1)
    norm_counts <- apply(tmp[["mat"]], 1:2, as.numeric) * 1.5  ## decimal
    scale_mat <- t(scale(t(norm_counts)))
    pheno <- tmp[["pheno"]]

    gene_cluster_scores <- data.frame(
        gene = rownames(scale_mat),
        cluster = c(rep(1, ngenes / 50), rep(2, ngenes / 50)),
        score = rnorm(ngenes, 0.5, 0.25)
    )
    score_thresh <- c(0.2, 0.3)
    genes_of_interest <- paste0("G", 5:50)
    res <- do_gene_plots(norm_counts, scale_mat, pheno,
                        gene_cluster_scores, score_thresh,
                        genes_of_interest,
                        out_dir = tempdir())
})

test_that("single_kmeans_heatmap", {
    n <- 100
    scale_mat <- matrix(rnorm(n*10, mean = 5), nrow = n)
    rownames(scale_mat) <- paste0("G", 1:n)
    colnames(scale_mat) <- paste0("S", 1:10)
    top_genes <- paste0("G", n - 50:n - 20)
    res <- single_kmeans_heatmap(scale_mat, 2, top_genes, "test", "test",
                                out_dir = tempdir(), "test", 7, 7)
})

test_that("high_quality_genes", {
    sam_mat <- matrix(1:100, nrow = 10,
                    dimnames = list(
                        paste0("G", seq(10)),
                        paste0("C", seq(10))
                    ))
    keep <- high_quality_genes(sam_mat, 10, 7, tempdir())
    names(keep) == paste0("G", 7:10)
})

test_that("do_volcanos", {
    n <- 100
    dsqres <- data.frame(
        gene = paste0("G", 1:n),
        data = rnorm(n, 5, 2),
        log2FoldChange = 1e-5,
        mLog10Padj = rnorm(n, 1000, 2),
        padj = 1/rnorm(n, 1000, 2)
    )
    top_genes_tohighlight <- paste0("G", (n-50):(n-20))
    res <- RNASeqHelper:::do_volcanos(dsqres, top_genes_tohighlight,
                                    "my plot", tempdir())
})

test_that("volcano_plot", {
    n <- 100
    dsqres <- data.frame(gene = paste0("G", 1:n),
                        data = rnorm(n, 5, 2),
                        log2FoldChange = 1e-5,
                        padj = 1 / rnorm(n, 5, 2))
    degenes <- paste0("G", (n-50):(n-20))
    res <- volcano_plot(dsqres, degenes, "my plot")
})

test_that("gene_clusters_by_score", {
    n <- 100
    tab <- data.frame(
        gene = paste0("G", 1:n),
        cluster = as.integer(rnorm(n, 5, 2) + 1),
        condition = c(rep("red", n/2), rep("green", n/2)),
        value = rnorm(n, 10, 2) + 1,
        time = as.integer(rnorm(n, 2, 0.5) + 1) * 5,
        score = as.integer(rnorm(n, 2, 0.5) + 1) / 4
    )
    plot <- gene_clusters_by_score(tab, c(0.2, 0.7),
                                    out_dir = tempdir())
})

test_that("generate_labelling_table", {
    tabn <- data.frame(
        cluster = c(1, 1, 2, 2),
        gene = paste0("G", 1:4),
        value = 1:4
    )
    res <- generate_labelling_table(tabn)
})

test_that("num_clusters", {
    expect_identical(
        num_clusters(data.frame(
            cluster = c(1, 1, 2, 3, 2)
        )),
        3L
    )
})
